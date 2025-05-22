__version__ = "2025020401"
import configparser
import copy
import json
import shutil
from collections import defaultdict
from pathlib import Path

from accel import Atoms, Box, System, Systems
from accel.util import Execmd, Log
from accel.util.constants import Elements
from accel.util.matrix import Matrix
from acceltools.plot import PlotBox
from acceltools.series import SeriesBox

Execmd.add("sdconvert", r"/app/schrodinger/utilities/sdconvert")

pdir = Path().cwd()
tdir = pdir / "template"

config = configparser.ConfigParser()
config.read(pdir / "099_config" / "config.ini")
cfg = config["USER"]

Log.console(show=cfg.getboolean("console_show"))
Log.file(pdir.joinpath("_" + Path(__file__).stem).with_suffix(".out"))
Log.write(f"{Path(__file__).absolute()}: version {__version__}")
Log.write("Starting ACCeL NEB Step1")

if not cfg.getboolean("keep_calc"):
    Log.write("ACCeL NEB Step1 is deactivated")
    exit()

energy_limit_to_empopt = cfg.getfloat("energy_limit_to_empopt")
energy_limit_to_neb = cfg.getfloat("energy_limit_to_neb")
max_limit_to_neb = cfg.getint("max_limit_to_neb")
energy_limit_to_reapro_dft = cfg.getfloat("energy_limit_to_reapro_dft")
energy_limit_to_ts_dft = cfg.getfloat("energy_limit_to_ts_dft")
max_limit_to_dftopt = cfg.getint("max_limit_to_dftopt")
max_limit_for_local_reapro_in_total = cfg.getint("max_limit_for_local_reapro_in_total")

additional_valid_range = cfg.getint("additional_valid_range")
acceptable_invalid_atoms = cfg.getint("acceptable_invalid_atoms")

cov_scaling = cfg.getfloat("cov_scaling")

calc_symm_for_all_ts = cfg.getboolean("calc_symm_for_all_ts")

fix_isomeric_ez = cfg.getboolean("fix_isomeric_ez")

USE_CI_XYZ = cfg.getboolean("use_ci_xyz")

NEB_KEY = "NEB-TS"
if USE_CI_XYZ:
    NEB_KEY = "NEB-CI"
if cfg.getboolean("zoom_neb"):
    NEB_KEY = "ZOOM-" + NEB_KEY
XTB_KEY = "XTB2"
XTB_VER = 2
if cfg.getboolean("use_xtb1"):
    XTB_KEY = "XTB1"
    XTB_VER = 1

lewis = defaultdict(dict)
for line in json.loads(cfg["lewis"]):
    lewis[line[0]][(min(int(line[1]), int(line[2])), max(int(line[1]), int(line[2])))] = float(line[3])
lewis: dict[str, dict[tuple[int], float]] = dict(lewis)

lewis_f = cfg.getfloat("lewis_force_constant")

ignored_atoms_in_bond_check: dict[str, list[int]] = json.loads(cfg["ignored_atoms_in_bond_check"])
mm_const: dict[str, list[list]] = json.loads(cfg["mm_constraints"])
mm_const_str: dict[str, list[str]] = {}
for key, ls in mm_const.items():
    constr = ""
    for t in ls:
        if len(t) != 9:
            raise ValueError
        constr += (
            f"\n{t[0].replace(' ', ''):^6} {int(t[1]):>6} {int(t[2]):>6} {int(t[3]):>6} {int(t[4]):>6}"
            + f" {float(t[5]):>10.4f} {float(t[6]):>10.4f} {float(t[7]):>10.4f} {float(t[8]):>10.4f}"
        )
    mm_const_str[key] = constr

ignored_bonds_in_bond_check: dict[str, list[list[int]]] = json.loads(cfg["ignored_bonds_in_bond_check"])


def modify_bonds(box: Box, ref_box: Box):
    for label, cs in box.get().labels.items():
        lewis_dict = lewis.get(label.split("_")[1])
        if lewis_dict is not None:
            for c in cs:
                if c.data.get("diagram_role") == "ts":
                    continue
                for b in lewis_dict.keys():
                    c.atoms.bonds[b[0], b[1]] = 0
        ignored_atoms_list = ignored_atoms_in_bond_check.get(label.split("_")[1])
        if ignored_atoms_list is not None and ref_box.get().labels.get(label) is not None:
            ref_atoms = ref_box.get().labels[label].get().atoms
            for c in cs:
                if c.data.get("diagram_role") == "ts":
                    continue
                for a_num in ignored_atoms_list:
                    for not_bonding_atom in c.atoms.get(a_num).bonds:
                        c.atoms.bonds[a_num, not_bonding_atom.number] = 0
                    for bonding_atom in ref_atoms.get(a_num).bonds:
                        c.atoms.bonds[a_num, bonding_atom.number] = ref_atoms.bonds[a_num, bonding_atom.number]
        ignored_bonds_list = ignored_bonds_in_bond_check.get(label.split("_")[1])
        if ignored_bonds_list is not None and ref_box.get().labels.get(label) is not None:
            ref_atoms = ref_box.get().labels[label].get().atoms
            for c in cs:
                if c.data.get("diagram_role") == "ts":
                    continue
                for a_nums in ignored_bonds_list:
                    c.atoms.bonds[a_nums[0], a_nums[1]] = ref_atoms.bonds[a_nums[0], a_nums[1]]
    return None


# direction of pair generation: "both", "reactant_to_product", "product_to_reactant", None, "r2p", "p2r"
path_flag = None
if cfg["path_flag"].lower() in ("auto", "none"):
    path_flag = None
else:
    if cfg["path_flag"].lower() in ("both", "reactant_to_product", "product_to_reactant"):
        path_flag = cfg["path_flag"].lower()
    elif cfg["path_flag"].lower() == "r2p":
        path_flag = "reactant_to_product"
    elif cfg["path_flag"].lower() == "p2r":
        path_flag = "product_to_reactant"
    else:
        Log.write(f"could not detect path_flag: {cfg['path_flag']}")
        Log.write("path_flag will be determined automatically")


input_box = Box()
for role_dir in ["011_reactant", "012_product"]:
    for p in (pdir / role_dir).glob(pattern="*_*.mae"):
        c = System(file_path=p)
        c.data["role"] = role_dir.split("_")[1]
        input_box.add(c)

stem_name = input_box.get().get().name.split("_")[0]

input_box.labeling().read_atoms().calc_symm(calc_all=True)

input_box.calc_stereo()

rxn_bonds: dict[tuple[str], set[tuple[int]]] = {}
rxn_atoms: dict[tuple[str], set[int]] = {}
for _rea_name in [_c.name for _c in input_box.get().has_data("role", "reactant")]:
    for _pro_name in [_c.name for _c in input_box.get().has_data("role", "product")]:
        _bs = input_box.get().get(_rea_name).atoms.bonds.to_set() ^ input_box.get().get(_pro_name).atoms.bonds.to_set()
        rxn_bonds[(_rea_name, _pro_name)] = _bs
        rxn_atoms[(_rea_name, _pro_name)] = set()
        for b in _bs:
            rxn_atoms[(_rea_name, _pro_name)].add(b[0])
            rxn_atoms[(_rea_name, _pro_name)].add(b[1])
Log.write(f"reacting bonds: {rxn_bonds}")
Log.write(f"reacting atoms: {rxn_atoms}")

rxn_bonds_list = []
for rp_label in rxn_bonds:
    for b in rxn_bonds[rp_label]:
        rxn_bonds_list.append([rp_label[0], rp_label[1], b[0], b[1]])
config["DEFAULT"]["rxn_bonds"] = json.dumps(rxn_bonds_list)

rxn_atoms_list = []
for rp_label in rxn_atoms:
    for a in rxn_atoms[rp_label]:
        rxn_atoms_list.append([rp_label[0], rp_label[1], a])
config["DEFAULT"]["rxn_atoms"] = json.dumps(rxn_atoms_list)

with pdir.joinpath("030_solvent").joinpath("solvent_key").open() as f:
    solvent_key = f.readlines()[0]

ignored_symbols_in_bond_check: dict[str, list[str]] = json.loads(cfg["ignored_symbols_in_bond_check"])
for label, cs in input_box.get().labels.items():
    symbols = ignored_symbols_in_bond_check.get(label.split("_")[1])
    if symbols is None:
        continue
    num_list: list[int] = []
    for symb in symbols:
        symb = Elements.canonicalize(symb)
        for a in cs.get().atoms:
            if a.symbol == symb:
                num_list.append(a.number)
    Log.write(f"ignored atoms added by symbol:{label}: {num_list}")
    ignored_atoms_in_bond_check[label.split("_")[1]] = (
        ignored_atoms_in_bond_check.get(label.split("_")[1], []) + num_list
    )


input_box.copy_files("100_ff_reapro_cs", True)

input_box.set_data("LEWIS", "")
for c in input_box.get():
    if lewis.get(c.label.split("_")[1]) is not None:
        for t, fc in lewis.get(c.label.split("_")[1]).items():
            c.data[
                "LEWIS"
            ] += f"\n FXDI {t[0]:>7}{t[1]:>7}      0      0 {lewis_f:>10.4f} {fc:>10.4f}     0.0000     0.0000"
    if mm_const_str.get(c.label.split("_")[1]) is not None:
        c.data["LEWIS"] += mm_const_str.get(c.label.split("_")[1])

auop_str = ""
auop_int = cfg.getint("auop")
if auop_int != 0:
    auop_str = f"\n AUOP       0      0      0      0 {auop_int:>5}.0000     0.0000     0.0000     0.0000"

input_box.write_input(
    tdir / "mm_cs.com", "100_ff_reapro_cs", arg={"MCMM": f"{cfg.getint('cs_steps'):>6}", "AUOP": auop_str}
)

input_box.write_input(tdir / "mm_cs.sbc", "100_ff_reapro_cs")
input_box.write_input(tdir / "mm.qsh", "100_ff_reapro_cs")
cs_box = Box("100_ff_reapro_cs/*.qsh").run()

for c in cs_box.get():
    c.name += "-out"
rp_box = Box(cs_box.search(suffix=".maegz").plugin.mae.get_unzip())
rp_box.plugin.mae.calc_energy()
rp_box.labeling()
for c in rp_box.get():
    c.data["role"] = input_box.get().labels[c.label].get().data["role"]

rp_box.energy_limit(energy_limit_to_empopt)
rp_box.map_numbers(input_box)
rp_box.calc_bonds(cov_scaling)
rp_box.write_xyz(pdir / "102_ff_reapro_xyz")

if cfg.get("charge") is not None:
    for c in rp_box.get():
        c.charge = cfg.getint("charge")

UHF_KEY = ""
if cfg.get("multiplicity") is not None:
    UHF_KEY = "--uhf " + str(cfg.getint("multiplicity") - 1)
    for c in rp_box.get():
        c.multiplicity = cfg.getint("multiplicity")


if cfg.getboolean("prep_opt"):
    for label, cs in rp_box.get().labels.items():
        const_atom_set: set[int] = set()
        for rp_label in rxn_atoms:
            if label in rp_label:
                const_atom_set = const_atom_set | rxn_atoms[rp_label]
        prep_const = "Constraints\n"
        for anum in const_atom_set:
            prep_const += "{C " + str(anum - 1) + " C}\n"
        prep_const += "end\n"
        Box(cs).set_data("PREP_CONST", prep_const)
    rp_box.write_input(tdir / "orca_opt_prep.inp", pdir / "109_emp_reapro_prep_opt", arg={"XTB_KEY": XTB_KEY})
    rp_box.write_input(tdir / "orca_opt.qsh", pdir / "109_emp_reapro_prep_opt")
    rp_box.run()
    rp_box.search(suffix=".out").check_end().read_atoms()

if cfg.getboolean("xtb_opt"):
    rp_box.write_xyz(pdir / "110_emp_reapro_opt", centering=False)
    rp_box.write_input(
        tdir / "xtb_opt.qsh", pdir / "110_emp_reapro_opt", arg={"XTB_VER": str(XTB_VER), "UHF_KEY": UHF_KEY}
    )
else:
    rp_box.write_input(tdir / "orca_opt.inp", pdir / "110_emp_reapro_opt", arg={"XTB_KEY": XTB_KEY})
    rp_box.write_input(tdir / "orca_opt.qsh", pdir / "110_emp_reapro_opt")
rp_box.run()

rp_box.search(suffix=".out").check_end().read_energy()
if cfg.getboolean("xtb_opt"):
    for c in rp_box.get():
        c.path = c.path.with_suffix(".xtbopt.xyz")
rp_box.read_atoms().calc_bonds(cov_scaling)

rot_mats_dict: dict[str, list[Matrix]] = {c.label: c.data["rotamer"] for c in input_box.get()}
num_mats_dict: dict[str, list[Matrix]] = {c.label: c.data["numisomer"] for c in input_box.get()}
for c in rp_box.get():
    atoms_list = c.atoms.to_list()
    c.data["rotamer"] = [Matrix(atoms_list).bind(m._matrix) for m in rot_mats_dict[c.label]]
    c.data["numisomer"] = [Matrix(atoms_list).bind(m._matrix) for m in num_mats_dict[c.label]]
    c.data["has_symm"] = True

rp_box.rmsd_limit(cfg.getfloat("rmsd_limit"))

modify_bonds(rp_box, input_box)

for c in rp_box.get():
    c_bonds_set = c.atoms.bonds.to_set()
    ref_bonds_set = input_box.get().labels[c.label].get().atoms.bonds.to_set()
    if c_bonds_set != input_box.get().labels[c.label].get().atoms.bonds.to_set():
        c.deactivate("bonds_check")
        Log.write(f"invalid bonds of {c.name}: {c_bonds_set}")
        Log.write(f"reference bonds: {ref_bonds_set}")
        Log.write(f"extra: {c_bonds_set - ref_bonds_set}: lack: {ref_bonds_set - c_bonds_set}")

alive_rp = rp_box.get()

rp_box.energy_limit(energy_limit_to_neb)
if path_flag is None:
    rea_lens = [len(_b) for _b in rp_box.get().has_data("role", "reactant").labels.values()]
    rea_ave = sum(rea_lens) / len(rea_lens)
    pro_lens = [len(_b) for _b in rp_box.get().has_data("role", "product").labels.values()]
    pro_ave = sum(pro_lens) / len(pro_lens)
    if (abs(rea_ave - pro_ave) / (rea_ave + pro_ave)) < 0.25:
        path_flag = "both"
    elif rea_ave < pro_ave:
        path_flag = "reactant_to_product"
    else:
        path_flag = "product_to_reactant"
Log.write(f"path_flag: {path_flag}")

rp_box.energy_limit(energy_limit_to_neb, max_limit=max_limit_to_neb)

if len(rp_box.get()) > max_limit_for_local_reapro_in_total:
    config["DEFAULT"]["need_hpc"] = "True"

g16_qsh_file = "g16.qsh"
rqbs_port = cfg.getint("rqbs_port")
if 49152 <= rqbs_port and rqbs_port <= 65535:
    g16_qsh_file = "g16_rqbs.qsh"
    config["DEFAULT"]["need_hpc"] = "False"
if cfg.getboolean("need_hpc"):
    g16_qsh_file = "g16_hpc.qsh"

if cfg.getboolean("no_chk"):
    g16_qsh_file = g16_qsh_file.split(".")[0] + "_nochk.qsh"

g16_opt_file = "g16_opt.gjf"
g16_tsopt_file = "g16_tsopt.gjf"
if cfg.getboolean("qst_opt"):
    g16_tsopt_file = "g16_tsopt_qst.gjf"
if cfg.getboolean("sp_only"):
    g16_opt_file = "g16_sp.gjf"
    g16_tsopt_file = "g16_sp.gjf"

Box(rp_box.get().has_data("role", "reactant")).write_xyz("111_emp_rea_xyz")
Box(rp_box.get().has_data("role", "product")).write_xyz("112_emp_pro_xyz")

for c in rp_box.get():
    if c.data["role"] == "reactant":
        if path_flag == "both" or path_flag == "reactant_to_product":
            c.data["rxn_mates"] = input_box.get().has_data("role", "product")
        else:
            c.deactivate("neb_path_detection")
    elif c.data["role"] == "product":
        if path_flag == "both" or path_flag == "product_to_reactant":
            c.data["rxn_mates"] = input_box.get().has_data("role", "reactant")
        else:
            c.deactivate("neb_path_detection")
    else:
        c.deactivate("neb_path_detection")

ts_box = Box()
for c in rp_box.get():
    for input_c in c.data["rxn_mates"]:
        input_c: System = input_c
        ts_c = c.duplicate()
        ts_c.name = c.name + "_" + input_c.name
        ts_c.label = input_c.name
        ts_c.atoms.bonds.bind(copy.deepcopy(input_c.atoms.bonds._matrix))
        ts_c.atoms._copy_charges(input_c.atoms)
        ts_c.data["rxn_mates"] = Systems().bind([c])
        ts_c.data["LEWIS"] = input_c.data["LEWIS"]
        atoms_list = ts_c.atoms.to_list()
        ts_c.data["rotamer"] = [Matrix(atoms_list).bind(_m._matrix) for _m in input_c.data["rotamer"]]
        ts_c.data["numisomer"] = [Matrix(atoms_list).bind(_m._matrix) for _m in input_c.data["numisomer"]]
        ts_c.data["has_symm"] = True
        ts_box.add(ts_c)

ts_box.set_data("CONSTRAINT", "")
ts_box.plugin.mae.write_mae("120_ff_pair_mini")
ts_box.write_input(tdir / "mm_mini.com", "120_ff_pair_mini")
ts_box.write_input(tdir / "mm_mini.sbc", "120_ff_pair_mini")
ts_box.write_input(tdir / "mm.qsh", "120_ff_pair_mini")
ts_box.run()

for c in ts_box.get():
    c.data["original_name"] = c.name
    c.name += "-out"
ts_box.search(suffix=".maegz").plugin.mae.read_maegz()
for c in ts_box.get():
    c.name = c.data["original_name"]

ts_box.set_data("db_chirality", False)
for c in ts_box.get():
    c.data["db_failed_dict"] = defaultdict(lambda: 0)


def get_transs(atoms: Atoms) -> list[tuple[int]]:
    transs = []
    for b_form_to, b_type in atoms.bonds.to_dict().items():
        if b_type == 2:
            atom_a = atoms.get(b_form_to[0])
            atom_b = atoms.get(b_form_to[1])
            if len(atom_a.bonds) > len(atom_b.bonds):
                atom_a, atom_b = atom_b, atom_a
            if len(atom_a.bonds) >= 2:
                if (not fix_isomeric_ez) and len(atom_b.bonds) <= 1:
                    continue
                for nb_a in [a for a in atom_a.single if a is not atom_b]:
                    for nb_b in [b for b in atom_b.single if b is not atom_a]:
                        dihed_ang = atoms.get_dihedral(nb_a.number, atom_a.number, atom_b.number, nb_b.number)
                        if dihed_ang > 90 or dihed_ang < -90:
                            transs.append((nb_a.number, atom_a.number, atom_b.number, nb_b.number))
    return transs


for force_index in range(cfg.getint("max_cycle_of_force_mini")):
    for c in ts_box.get().has_data("db_chirality", False):
        calcd_transs = get_transs(c.atoms)
        failed_transs = [t for t in get_transs(input_box.get().get(c.label).atoms) if t not in calcd_transs]
        if len(failed_transs) == 0:
            c.data["db_chirality"] = True
            continue
        for t in failed_transs:
            c.data["db_failed_dict"][t] += 50
        c.data["CONSTRAINT"] = ""
        for t, fc in c.data["db_failed_dict"].items():
            c.data[
                "CONSTRAINT"
            ] += f"\n FXTA {t[0]:>7}{t[1]:>7}{t[2]:>7}{t[3]:>7} {round(fc):>5}.0000   180.0000     0.0000     0.0000"

    recalc_box = Box(ts_box.get().has_data("db_chirality", False))
    if len(recalc_box) == 0:
        Log.write("No invalid double bond chirality")
        break
    else:
        Log.write(f"invalid double bond chirality detected in {len(recalc_box)} structures")

    recalc_box.search("120_ff_pair_mini", suffix=".mae").read_atoms()
    recalc_box.plugin.mae.write_mae(f"120_ff_pair_mini_recalc{force_index}")
    recalc_box.write_input(tdir / "mm_mini.com", f"120_ff_pair_mini_recalc{force_index}")
    recalc_box.write_input(tdir / "mm_mini.sbc", f"120_ff_pair_mini_recalc{force_index}")
    recalc_box.write_input(tdir / "mm.qsh", f"120_ff_pair_mini_recalc{force_index}")
    recalc_box.run()

    for c in recalc_box.get():
        c.name += "-out"
    recalc_box.search(suffix=".maegz").plugin.mae.read_maegz()
    for c in recalc_box.get():
        c.name = c.data["original_name"]


ts_box.write_xyz("121_after_pair_mini", centering=False)

if cfg.getboolean("pair_prep_opt"):
    for label, cs in ts_box.get().labels.items():
        const_atom_set: set[int] = set()
        for rp_label in rxn_atoms:
            if label in rp_label:
                const_atom_set = const_atom_set | rxn_atoms[rp_label]
        prep_const = "Constraints\n"
        for anum in const_atom_set:
            prep_const += "{C " + str(anum - 1) + " C}\n"
        prep_const += "end\n"
        Box(cs).set_data("PREP_CONST", prep_const)
    ts_box.write_input(tdir / "orca_opt_prep.inp", pdir / "121_emp_pair_prep_opt", arg={"XTB_KEY": XTB_KEY})
    ts_box.write_input(tdir / "orca_opt.qsh", pdir / "121_emp_pair_prep_opt")
    ts_box.run()
    ts_box.search(suffix=".out").check_end().read_atoms()

if cfg.getboolean("xtb_opt"):
    ts_box.write_xyz(pdir / "122_emp_pair_opt", centering=False)
    ts_box.write_input(
        tdir / "xtb_opt.qsh", pdir / "122_emp_pair_opt", arg={"XTB_VER": str(XTB_VER), "UHF_KEY": UHF_KEY}
    )
else:
    ts_box.write_input(tdir / "orca_opt.inp", pdir / "122_emp_pair_opt", arg={"XTB_KEY": XTB_KEY})
    ts_box.write_input(tdir / "orca_opt.qsh", pdir / "122_emp_pair_opt")
ts_box.run()

ts_box.search(suffix=".out").check_end()
if cfg.getboolean("xtb_opt"):
    for c in ts_box.get():
        c.path = c.path.with_suffix(".xtbopt.xyz")
ts_box.read_atoms().calc_bonds(cov_scaling)
modify_bonds(ts_box, input_box)

if cfg.getboolean("check_chirality"):
    ts_box.calc_stereo()
    for c in ts_box.get():
        for a, input_atom in zip(c.atoms, input_box.get().get(c.label).atoms):
            all_stereo: set = a.stereo | input_atom.stereo
            if "nR" in all_stereo and "nS" in all_stereo:
                c.deactivate(f"{c.name}: invalid RS chirality")
                continue
            if "nE" in all_stereo and "nZ" in all_stereo:
                c.deactivate(f"{c.name}: invalid EZ chirality")
                continue

for c in ts_box.get():
    rxn_mate_mols: Systems = c.data["rxn_mates"]
    c.data["mate_name"] = rxn_mate_mols.get().name
    ref_bonds_set = input_box.get().labels["_".join(c.name.split("_")[3:5])].get().atoms.bonds.to_set()
    c_bonds_set = c.atoms.bonds.to_set()
    if c_bonds_set != ref_bonds_set:
        c.deactivate("bonds_check")
        Log.write(f"invalid bonds of {c.name}: {c_bonds_set}")
        Log.write(f"reference bonds: {ref_bonds_set}")
        Log.write(f"extra: {c_bonds_set - ref_bonds_set}: lack: {ref_bonds_set - c_bonds_set}")


(pdir / "113_rea_pro_data").mkdir(exist_ok=True)
rp_box.export_data(pdir / "113_rea_pro_data" / stem_name)
rp_box.write_xyz(pdir / "123_emp_ts_neb", centering=False)
ts_box.write_input(
    tdir / "orca_neb.inp",
    pdir / "123_emp_ts_neb",
    arg={"XTB_KEY": XTB_KEY, "NEB_KEY": NEB_KEY, "NIMAGES": str(cfg.getint("nimages"))},
)
ts_box.write_input(tdir / "orca_neb.qsh", pdir / "123_emp_ts_neb")
ts_box.run()

for c in ts_box.get():
    c.data["diagram_role"] = "ts"
    connection = ["_".join(c.name.split("_")[0:2]), "_".join(c.name.split("_")[3:5])]
    if connection[0] not in input_box.get().has_data("role", "reactant").labels.keys():
        connection[0], connection[1] = connection[1], connection[0]
    c.data["diagram_connection"] = connection
    c.label = "_".join(connection)

ts_box.search(suffix=".out").check_end()
if not USE_CI_XYZ:
    ts_box.read_energy()
    ts_box.read_atoms()
else:
    for c in ts_box.get():
        c.name = c.name + "_NEB-CI_converged"
    ts_box.search(suffix=".xyz").read_atoms()
    for c in ts_box.get():
        c.name = "_".join(c.name.split("_")[:5])
    ts_box.search(suffix=".out")
ts_box.calc_bonds(cov_scaling).calc_symm(calc_all=calc_symm_for_all_ts).rmsd_limit(cfg.getfloat("rmsd_limit"))
modify_bonds(ts_box, input_box)

idist = cfg.getint("image_distance")
eidx = 2 if not USE_CI_XYZ else 3
for c in ts_box.get():
    with c.path.open() as f:
        ls = f.readlines()
    out_of_summary = True
    neb_summary = None
    for image_idx, line in enumerate(ls):
        if ("Image     E(Eh)" if not USE_CI_XYZ else "Image Dist.(Ang.)") in line:
            out_of_summary = False
            neb_summary = []
            continue
        if out_of_summary:
            continue
        if not line.startswith(" "):
            out_of_summary = True
            continue
        sp_l = line.split()
        if sp_l[0] == "TS":
            neb_summary.append(("TS", float(sp_l[eidx])))
            ts_energy = float(sp_l[eidx])
            ts_position = len(neb_summary) - 1
        elif USE_CI_XYZ and sp_l[-1] == "CI":
            neb_summary.append((int(sp_l[0]), float(sp_l[eidx])))
            c.energy = float(sp_l[2])
            neb_summary.append(("TS", float(sp_l[eidx])))
            ts_energy = float(sp_l[eidx])
            ts_position = len(neb_summary) - 1
        else:
            neb_summary.append((int(sp_l[0]), float(sp_l[eidx])))
    if neb_summary is None:
        c.deactivate("NEB summary was not found.")
        continue
    Log.write(f"{c.name}: NEB: {neb_summary}")
    before_ts = True
    r_idx = 0
    p_idx = len(neb_summary) - 2
    for image_idx, energy in neb_summary:
        if image_idx == "TS":
            before_ts = False
            continue
        if energy > ts_energy or ((ts_position - idist) <= image_idx and image_idx <= (ts_position + (idist - 1))):
            continue
        if before_ts:
            r_idx = image_idx
        else:
            p_idx = image_idx
            break
    trj_path = c.path.with_name(f"{c.path.stem}_MEP_trj.xyz")
    if not trj_path.exists():
        c.deactivate("MEP_trj file was not found.")
        continue
    with trj_path.open() as f:
        ls = f.readlines()
    ls = [line.replace("\n", "") for line in ls if line.startswith("  ")]
    c.data["RAXYZ"] = "\n".join(ls[(len(c.atoms) * r_idx) : (len(c.atoms) * (r_idx + 1))])
    c.data["PAXYZ"] = "\n".join(ls[(len(c.atoms) * p_idx) : (len(c.atoms) * (p_idx + 1))])


for label, cs in ts_box.get().labels.items():
    rxn_tuple = ("_".join(label.split("_")[0:2]), "_".join(label.split("_")[2:4]))
    for b in rxn_bonds[rxn_tuple]:
        Box().bind(cs).calc_length(b[0], b[1])


ts_box.write_xyz(pdir / "124_emp_nebts_xyz")

ts_box.write_xyz(pdir / "125_emp_ts_freq", centering=False)
ts_box.write_input(tdir / "xtb_freq.qsh", pdir / "125_emp_ts_freq", arg={"XTB_VER": str(XTB_VER), "UHF_KEY": UHF_KEY})
ts_box.run()
ts_box.search(suffix=".out").check_end()

ts_box.search(suffix=".g98.out")

ts_box.copy_files(pdir / "127_emp_ts_freq_analysis", suffix=".g98")
if len(ts_box.get()) != 0:
    freq_box = Box(pdir / "127_emp_ts_freq_analysis" / "*.g98")
    freq_box.plugin.txt.read_text()
    freq_box.plugin.txt.delete_lines_by_key("Raman Activ")
    freq_box.plugin.txt.delete_lines_by_key("Depolar")
    freq_box.plugin.txt.write_text(suffix=".g98")

if USE_CI_XYZ:
    Log.write("skipping check_freq because use_ci_xyz is activated")
else:
    ts_box.plugin.gau.check_freq(im=1)
ts_box.plugin.gau.read_vibration()

irc_box = Box()
if cfg.getboolean("irc_test"):
    ts_box.write_input(
        tdir / "orca_irc.inp", pdir / "125_emp_ts_irc_F", arg={"DIRECTION": "forward", "XTB_KEY": XTB_KEY}, link=False
    )
    ts_box.write_input(tdir / "orca_irc.qsh", pdir / "125_emp_ts_irc_F", link=False)
    ts_box.write_input(
        tdir / "orca_irc.inp", pdir / "125_emp_ts_irc_B", arg={"DIRECTION": "backward", "XTB_KEY": XTB_KEY}, link=False
    )
    ts_box.write_input(tdir / "orca_irc.qsh", pdir / "125_emp_ts_irc_B", link=False)
    irc_box.add(Box(pdir / "125_emp_ts_irc_F" / "*.qsh").set_data("irc", "IRC_F"))
    irc_box.add(Box(pdir / "125_emp_ts_irc_B" / "*.qsh").set_data("irc", "IRC_B"))
    irc_box.run()

    def get_trjxyz(c: System, increment_num: int) -> Systems:
        with c.path.open() as f:
            ls = f.readlines()
        is_title_line = False
        new_systems = Systems()
        current_idx = 0
        for line in ls:
            l_split = line.split()
            if (not is_title_line) and len(l_split) == 1 and l_split[0].isnumeric():
                is_title_line = True
            elif is_title_line:
                new_c = System()
                current_idx += increment_num
                new_c.name = f"{c.name}_{current_idx:03}"
                new_c.data["irc_idx"] = current_idx
                new_c.data["xyz_title"] = line.replace("\n", "")
                new_systems.append(new_c)
                is_title_line = False
            elif len(line.split()) == 4:
                new_c.atoms.append(line.split())
        return new_systems

    irc_trj = Box()
    irc_trj.add(Box(pdir / "125_emp_ts_irc_F" / "*_trj.xyz").set_data("increment_num", 1))
    irc_trj.add(Box(pdir / "125_emp_ts_irc_B" / "*_trj.xyz").set_data("increment_num", -1))
    irc_trj_all = Box(pdir / "124_emp_nebts_xyz" / "*.xyz").set_data("irc_idx", 0).read_atoms()
    for c in irc_trj.get():
        irc_trj_all.add(get_trjxyz(c, c.data["increment_num"]))
    irc_trj_all.labeling(index_list=[0, 1, 2, 3, 4])

    for label, cs in irc_trj_all.get().labels.items():
        if len(cs) == 1:
            for c in cs:
                c.deactivate("IRC_NOT_PERFORMED")

    SeriesBox(irc_trj_all).write_trjxyz(pdir / "128_emp_irc_analysis", "irc_idx")

    for c in irc_box.get():
        c.name += "_" + c.data["irc"]
    irc_box.search(suffix=".xyz").labeling(index_list=[0, 1, 2, 3, 4])
    irc_box.read_atoms()
    for label, cs in irc_box.get().labels.items():
        if len(cs) != 2:
            for c in cs:
                c.deactivate("could not get irc pair")
            continue
        rp_bonds_set: dict[str, set] = {}
        a_label = "_".join(label.split("_")[0:2])
        a_mol = input_box.get().get(a_label)
        rp_bonds_set[a_mol.data["role"]] = a_mol.atoms.bonds.to_set()
        b_label = "_".join(label.split("_")[3:5])
        b_mol = input_box.get().get(b_label)
        rp_bonds_set[b_mol.data["role"]] = b_mol.atoms.bonds.to_set()
        temp_box = Box(cs)
        for p_label, q_label in [(a_label, b_label), (b_label, a_label)]:
            Log.logger.debug("checking bonding supposing the following structures")
            cs[0].label = p_label
            cs[1].label = q_label
            temp_box.calc_bonds(cov_scaling)
            modify_bonds(temp_box, input_box)
            for c in cs:
                bonds_set = c.atoms.bonds.to_set()
                is_r = rp_bonds_set["reactant"] == bonds_set
                is_p = rp_bonds_set["product"] == bonds_set
                if is_r ^ is_p:
                    if is_r:
                        c.data["role"] = "reactant"
                    else:
                        c.data["role"] = "product"
                else:
                    c.data["role"] = "unknown"
            if (cs[0].data.get("role"), cs[1].data.get("role")) in (("reactant", "product"), ("product", "reactant")):
                if not cfg.getboolean("check_irc_chirality"):
                    break
                temp_box.calc_stereo()
                valid_chirality = True
                for c in cs:
                    for a, input_atom in zip(c.atoms, input_box.get().get(c.label).atoms):
                        all_stereo: set = a.stereo | input_atom.stereo
                        if "nR" in all_stereo and "nS" in all_stereo:
                            Log.write(f"{c.name}: invalid RS chirality: {str(a)}")
                            valid_chirality = False
                        if "nE" in all_stereo and "nZ" in all_stereo:
                            Log.write(f"{c.name}: invalid EZ chirality: {str(a)}")
                            valid_chirality = False
                if valid_chirality:
                    break
        else:
            for c in cs:
                c.deactivate("inappropriate terminal bonding")
    irc_box.labeling(index_list=[0, 1, 2, 3, 4])
    active_names = [c.label for c in irc_box.get()]
    for c in ts_box.get():
        if c.name in active_names:
            continue
        c.deactivate("inappropriate irc terminal")
else:
    for label, cs in ts_box.get().labels.items():
        rxn_tuple = ("_".join(label.split("_")[0:2]), "_".join(label.split("_")[2:4]))
        Box().bind(cs).plugin.gau.check_bonding(
            list(rxn_atoms[rxn_tuple]),
            additional_valid_range=additional_valid_range,
            acceptable_invalid_atoms=acceptable_invalid_atoms,
        )

if cfg.getboolean("use_irc_terminal_for_qst") and cfg.getboolean("irc_test"):
    for c in ts_box.get():
        irc_rp = irc_box.get().has_label(c.name)
        if len(irc_rp) != 2:
            Log.write(f"{c.name} :could not get irc terminals: invalid length {len(irc_rp)}")
            continue
        irc_r_c = irc_rp[0]
        irc_p_c = irc_rp[1]
        if (irc_r_c.data.get("role"), irc_p_c.data.get("role")) != ("reactant", "product"):
            irc_r_c, irc_p_c = irc_p_c, irc_r_c
        if (irc_r_c.data.get("role"), irc_p_c.data.get("role")) != ("reactant", "product"):
            Log.write(f"{c.name} :could not get irc terminals: invalid role")
            continue
        c.data["RAXYZ"] = "\n".join([" ".join([str(v) for v in a.axyz]) for a in irc_r_c.atoms])
        c.data["PAXYZ"] = "\n".join([" ".join([str(v) for v in a.axyz]) for a in irc_p_c.atoms])

ts_box.write_xyz(pdir / "126_emp_ts_xyz")
ts_box.energy_limit(energy_limit_to_ts_dft, max_limit=max_limit_to_dftopt)
ts_box.write_input(
    tdir / g16_tsopt_file,
    pdir / "203_dft_ts_opt",
    arg={"SOLVENT_KEY": solvent_key, "MODEL_OPT": cfg["model_opt"], "MODEL_SP": cfg["model_sp"]},
)
ts_box.write_input(tdir / g16_qsh_file, pdir / "203_dft_ts_opt", arg={"PORT": rqbs_port})
ts_box.search(pdir / "125_emp_ts_freq", suffix=".out").plugin.xtb.read_free_energy()

(pdir / "126_emp_ts_data").mkdir(exist_ok=True)
ts_box.export_data(pdir / "126_emp_ts_data" / stem_name)

low_rp_box = Box(alive_rp).set_state(True).energy_limit(energy_limit_to_reapro_dft, max_limit=max_limit_to_dftopt)
low_r_box = Box(low_rp_box.get().has_data("role", "reactant"))
low_r_box.write_input(
    tdir / g16_opt_file,
    pdir / "201_dft_rea_opt",
    arg={"SOLVENT_KEY": solvent_key, "MODEL_OPT": cfg["model_opt"], "MODEL_SP": cfg["model_sp"]},
)
low_r_box.write_input(tdir / g16_qsh_file, pdir / "201_dft_rea_opt", arg={"PORT": rqbs_port})
low_p_box = Box(low_rp_box.get().has_data("role", "product"))
low_p_box.write_input(
    tdir / g16_opt_file,
    pdir / "202_dft_pro_opt",
    arg={"SOLVENT_KEY": solvent_key, "MODEL_OPT": cfg["model_opt"], "MODEL_SP": cfg["model_sp"]},
)
low_p_box.write_input(tdir / g16_qsh_file, pdir / "202_dft_pro_opt", arg={"PORT": rqbs_port})
low_rp_box.write_xyz(pdir / "114_emp_reapro_freq", centering=False)
low_rp_box.write_input(
    tdir / "xtb_freq.qsh", pdir / "114_emp_reapro_freq", arg={"XTB_VER": str(XTB_VER), "UHF_KEY": UHF_KEY}
)
low_rp_box.run()
low_rp_box.search(suffix=".out").check_end().plugin.xtb.read_free_energy()

for c in low_rp_box.get():
    c.data["diagram_role"] = c.data["role"]

ts_labels: dict[str, str] = {}
for c in ts_box.get():
    ts_labels[c.name] = c.label
config["DEFAULT"]["ts_labels"] = json.dumps(ts_labels)

plot_mc = Box(ts_box.get() + low_rp_box.get())
(pdir / "130_energy_diagram_xtb").mkdir(exist_ok=True)
PlotBox(plot_mc).diagram(pdir / "130_energy_diagram_xtb" / stem_name)
plot_mc.only_minimum().export_data(pdir / "130_energy_diagram_xtb" / stem_name)

if cfg.getboolean("need_hpc"):
    shutil.copy(tdir / "atm.py", pdir / "201_dft_rea_opt" / "atm.py")
    shutil.copy(tdir / "qsb", pdir / "201_dft_rea_opt" / "qsb")
    shutil.copy(tdir / "atm.py", pdir / "202_dft_pro_opt" / "atm.py")
    shutil.copy(tdir / "qsb", pdir / "202_dft_pro_opt" / "qsb")
    shutil.copy(tdir / "atm.py", pdir / "203_dft_ts_opt" / "atm.py")
    shutil.copy(tdir / "qsb", pdir / "203_dft_ts_opt" / "qsb")

with (pdir / "099_config" / "config.ini").open("w") as f:
    config.write(f)

Log.write("ACCel NEB Step1 terminated normally")
