__version__ = "2025020401"
import configparser
import json
from collections import defaultdict
from pathlib import Path

from accel import Box, System
from accel.util import Execmd, Log
from accel.util.constants import Elements
from rdkit import Chem
from rdkit.Chem import AllChem

Execmd.add("sdconvert", r"/app/schrodinger/utilities/sdconvert")

config = configparser.ConfigParser()

pdir = Path().cwd()
with (pdir / "template" / "default.ini").open() as f:
    config_ls = f.readlines()
config_ls.append("\n[USER]\n")
for p in pdir.glob("*.ini"):
    with p.open() as f:
        config_ls.extend(f.readlines())
    config_ls.append("\n")
config.read_string("".join(config_ls))
cfg = config["USER"]

max_mw = cfg.getint("max_mw")
lewis_distance_scaling = cfg.getfloat("lewis_distance_scaling")

Log.console(show=cfg.getboolean("console_show"))
Log.file(pdir.joinpath("_" + Path(__file__).stem).with_suffix(".out"))
Log.write(f"{Path(__file__).absolute()}: version {__version__}")
Log.write("Starting ACCeL NEB Step0")

alphabets = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

solvents = {
    "chloroform": ["chloroform", "chcl3", "cdcl3"],
    "water": ["water", "h2o", "d2o"],
    "ethanol": ["ethanol", "etoh", "etod"],
    "methanol": ["methanol", "meoh", "meod"],
    "dimethylsulfoxide": ["dimethylsulfoxide", "dmso"],
    "benzene": ["benzene", "c6h6", "c6d6", "phh"],
    "toluene": ["toluene", "c6h5ch3", "c6d5cd3", "phme"],
    "acetonitrile": ["acetonitrile", "mecn", "ch3cn", "cd3cn"],
    "acetone": ["acetone", "ch3coch3"],
    "tetrahydrofuran": ["tetrahydrofuran", "thf"],
    "diethylether": ["diethylether", "et2o", "ether"],
    "dichloromethane": ["dichloromethane", "dcm", "ch2cl2", "methylenechloride"],
    "dichloroethane ": ["dichloroethane", "ch2clch2cl", "dce"],
}

if (pdir / "Need_HPC.txt").exists() or (pdir / "Need_DFT_on_SuperComputer.txt").exists():
    cfg["need_hpc"] = "True"

for rxn_path in pdir.glob("*.rxn"):
    pdir.joinpath("reactant").mkdir(exist_ok=True)
    pdir.joinpath("product").mkdir(exist_ok=True)
    with rxn_path.open() as f:
        ls = f.readlines()
    rxn_label = ["reactant" for _ in range(int(ls[4][0:3]))]
    rxn_label = rxn_label + ["product" for _ in range(int(ls[4][3:6]))]
    rxn_mols = [[] for _ in range(len(rxn_label))]
    mols_idx = -1
    for line in ls:
        if "$MOL" in line:
            mols_idx += 1
            continue
        if mols_idx == -1:
            continue
        rxn_mols[mols_idx].append(line)
    for idx, rxn_mol in enumerate(rxn_mols):
        out_path = pdir / rxn_label[idx]
        out_path = out_path / (rxn_path.stem.split("_")[0] + "_" + alphabets[idx : idx + 1])
        with out_path.with_suffix(".mol").open("w", newline="\n") as f:
            f.writelines(rxn_mol)

inp_box = Box()
not_mol_box = Box()
for role in ["reactant", "product"]:
    for p in (pdir / role).glob(pattern="*_*.*"):
        c = System(file_path=p)
        c.data["role"] = role
        c.name = "_".join(c.name.split("_")[0:2])
        c.label = c.name
        if p.suffix != ".mol":
            not_mol_box.add(c)
            continue
        Box(c).read_atoms()
        for a in c.atoms:
            a.data["original_number"] = a.number
        with c.path.open() as f:
            ls = f.readlines()
        sdf_bonds_order = {}
        sdf_bonds_chiral = {}
        for line in ls[(4 + int(ls[3][0:3])) : (4 + int(ls[3][0:3]) + int(ls[3][3:6]))]:
            if int(line[6:9]) in [8, 9]:
                sdf_bonds_order[(int(line[0:3]), int(line[3:6]))] = int(line[6:9])
            if int(line[9:12]) != 0:
                sdf_bonds_chiral[(int(line[0:3]), int(line[3:6]))] = int(line[9:12])
        if c.atoms.bonds is None:
            splitted_atoms = c.atoms
        else:
            splitted_atoms = c.modeler.get_splitted()
        if len(splitted_atoms) == 1 and len(sdf_bonds_order) == 0:
            inp_box.add(c)
            c.data["lewis_idxs"] = [[] for _ in c.atoms]
            for idx, bond in enumerate(sdf_bonds_order.keys()):
                for a in c.atoms:
                    if a.number not in bond:
                        continue
                    if sdf_bonds_order[bond] == 9:
                        c.data["lewis_idxs"][a.number - 1].append(idx)
        else:
            for atoms in splitted_atoms:
                atom_map = [a.data["original_number"] for a in atoms]
                sc = System(file_path=p)
                sc.data["role"] = role
                sc.label = c.label
                sc.atoms = atoms
                if len(splitted_atoms) == 1:
                    sc.name = c.name
                else:
                    sc.name = f"{c.name}{splitted_atoms.index(atoms) + 1}"
                Box(sc).write_mol("002_splitted_mols")
                with sc.path.open() as f:
                    ls = f.readlines()
                ls[1] = ls[1][:20] + "2" + ls[1][21:]
                for bond, chiral in sdf_bonds_chiral.items():
                    if bond[0] not in atom_map or bond[1] not in atom_map:
                        continue
                    for idx, line in enumerate(ls):
                        if line.startswith(f"{atom_map.index(bond[0])+1:>3}{atom_map.index(bond[1])+1:>3}"):
                            ls[idx] = f"{line[0:9]}{chiral:>3}{line[12:]}"
                        elif line.startswith(f"{atom_map.index(bond[1])+1:>3}{atom_map.index(bond[0])+1:>3}"):
                            if chiral == 1:
                                ls[idx] = f"{line[0:9]}{6:>3}{line[12:]}"
                            elif chiral == 6:
                                ls[idx] = f"{line[0:9]}{1:>3}{line[12:]}"
                            else:
                                ls[idx] = f"{line[0:9]}{chiral:>3}{line[12:]}"
                with sc.path.open("w") as f:
                    f.writelines(ls)
                sc.data["lewis_idxs"] = [[] for _ in atoms]
                for idx, bond in enumerate(sdf_bonds_order.keys()):
                    for a in atoms:
                        if a.data["original_number"] not in bond:
                            continue
                        if sdf_bonds_order[bond] == 9:
                            sc.data["lewis_idxs"][a.number - 1].append(idx)
                inp_box.add(sc)

for c in inp_box.get():
    rdmol = Chem.MolFromMolFile(str(c.path), removeHs=False)
    rdmol = AllChem.AddHs(rdmol)
    AllChem.EmbedMolecule(rdmol)
    try:
        AllChem.MMFFGetMoleculeForceField(rdmol, AllChem.MMFFGetMoleculeProperties(rdmol)).Minimize(maxIts=200)
    except AttributeError:
        AllChem.UFFGetMoleculeForceField(rdmol).Minimize(maxIts=200)
    c.data["txt"] = Chem.MolToMolBlock(rdmol)
    t_box = Box(c).copy_files("000_2d_mols")
    t_box.plugin.txt.write_text("001_3d_structures", True, suffix=".sdf")

inp_box.read_atoms()
for c in inp_box.get().has_data("lewis_idxs"):
    for idx, a in enumerate(c.atoms):
        try:
            a.data["lewis_idx"] = c.data["lewis_idxs"][idx]
        except IndexError:
            a.data["lewis_idx"] = None

Box(inp_box.get().has_bonds(False)).calc_bonds(cfg.getfloat("cov_scaling"))
for label, cs in inp_box.get().labels.items():
    if len(cs) == 1:
        continue
    csz = cs[0].duplicate()
    cs[0].deactivate()
    for c in cs[1:]:
        c: System = c
        csz.modeler.incorporate(c.atoms)
        c.deactivate()
    csz.name = label
    Box(csz).write_mol("002_splitted_mols")
    inp_box.add(csz)

box = Box(inp_box.get())
for c in not_mol_box.get():
    if c.path.suffix != ".mae":
        Box(c).read_atoms()
        box.add(c)
    else:
        if c.data["role"] == "reactant":
            Box(c).copy_files("011_reactant").read_atoms().write_mol("021_reactant_sdf")
        if c.data["role"] == "product":
            Box(c).copy_files("012_product").read_atoms().write_mol("022_product_sdf")

Box(box.get().has_bonds(False)).calc_bonds(cfg.getfloat("cov_scaling"))
if len(box) != 0:
    target_atoms = box.get().has_data("role", "reactant").get().atoms
    mol_wt_list = [c.atoms.mw for c in box.get()]
    if len(set(mol_wt_list)) != 1:
        Log.write("Check the input structures again.")
        config["DEFAULT"]["keep_calc"] = "False"
        exit()

mapped_mc = Box()
for c in box.get():
    mapped_atoms = c.modeler.get_mapped(target_atoms)
    if len(mapped_atoms) == 0:
        c.deactivate("No_mapping")
        continue
    if c.atoms != target_atoms:
        c.atoms = mapped_atoms[0]

    for idx, one_atoms in enumerate(mapped_atoms):
        new_c = c.duplicate()
        new_c.atoms = one_atoms
        new_c.name += f"_N{(idx + 1):02}"
        mapped_mc.add(new_c)

mapped_mc.write_mol("003_mapped_candidates", link=False)
box.write_mol("004_mapped_inputs", link=False, centering=False)
box.plugin.mae.write_mae("005_mae_inputs")

lewis = {}
for c in box.get():
    lewis_dict: dict[tuple[int], float] = {}
    lewis_idx_dict: dict[int, list[int]] = defaultdict(list)
    for a in c.atoms:
        if "lewis_idx" not in a.data:
            break
        if a.data["lewis_idx"] is None:
            continue
        for i in a.data["lewis_idx"]:
            lewis_idx_dict[i].append(a.number)
    else:
        for nums in lewis_idx_dict.values():
            if len(nums) != 2:
                continue
            vdw_dist = 0.0
            for n in (nums[0], nums[1]):
                t_vdw_dist = Elements.get_element(c.atoms.get(n).symbol)["vdw"]
                if t_vdw_dist is None:
                    Log.write(f"vdw radius of {c.atoms.get(n)} is unknown: 1.50 will be used")
                    vdw_dist += 1.50
                else:
                    vdw_dist += t_vdw_dist
            lewis_dict[(nums[0], nums[1])] = vdw_dist * lewis_distance_scaling
        lewis[c.name.split("_")[1]] = lewis_dict

Log.write(f"Lewis info: {lewis}")
lewis_lines = []
for n in lewis:
    for b, length in lewis[n].items():
        lewis_lines.append([n, b[0], b[1], round(length, 4)])
config["DEFAULT"]["lewis"] = json.dumps(lewis_lines)

inp_box.search("005_mae_inputs", suffix=".mae").read_atoms()

for c in inp_box.get():
    if c.data["role"] == "reactant":
        Box(c).copy_files("011_reactant").write_mol("021_reactant_sdf")
    if c.data["role"] == "product":
        Box(c).copy_files("012_product").write_mol("022_product_sdf")


mol_wt_list = [c.atoms.mw for c in Box([pdir / "011_reactant", pdir / "012_product"]).read_atoms().get()]
if len(set(mol_wt_list)) != 1:
    Log.write("Check the input structures again.")
    config["DEFAULT"]["keep_calc"] = "False"
    exit()

if max_mw < max(mol_wt_list):
    Log.write("Need_HPC")
    config["DEFAULT"]["need_hpc"] = "True"

solvent_files = list(pdir.glob("[Ss][Oo][Ll][Vv][Ee][Nn][Tt]*"))
with (pdir / "template" / "solvent_key.txt").open() as f:
    solvent_template = f.readlines()[0].replace("\n", "")
if len(solvent_files) == 0:
    solvent_keyword = " "
else:
    with solvent_files[0].open() as f:
        ls = f.readlines()
    for sol, names in solvents.items():
        for _n in names:
            if _n not in ls[0].lower():
                continue
            solvent_keyword = solvent_template.replace("#SOLVENT#", sol)
            break
        else:
            continue
        break
    else:
        t_sol_name = ls[0].replace(" ", "").replace("\n", "")
        Log.write("could not recognize solvent: {t_sol_name} will be used")
        solvent_keyword = solvent_template.replace("#SOLVENT#", t_sol_name)
sdir = pdir.joinpath("030_solvent")
sdir.mkdir(exist_ok=True)
with sdir.joinpath("solvent_key").open("w") as f:
    f.write(solvent_keyword)

(pdir / "099_config").mkdir(exist_ok=True)
with (pdir / "099_config" / "config.ini").open("w") as f:
    config.write(f)

Log.write("ACCel NEB Step0 terminated normally")
