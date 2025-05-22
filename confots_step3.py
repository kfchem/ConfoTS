__version__ = "2023051801"
import configparser
import json
import shutil
from collections import defaultdict
from pathlib import Path

from accel import Box
from accel.util import Log, Units
from acceltools.plot import PlotBox

pdir = Path().cwd()
config = configparser.ConfigParser()
config.read(pdir / "099_config" / "config.ini")
cfg = config["USER"]

additional_valid_range = cfg.getint("additional_valid_range")
acceptable_invalid_atoms = cfg.getint("acceptable_invalid_atoms")

Log.console(show=cfg.getboolean("console_show"))
Log.file(pdir.joinpath("_" + Path(__file__).stem).with_suffix(".out"))
Log.write(f"{Path(__file__).absolute()}: version {__version__}")
Log.write("Starting ACCeL NEB Step3")

if not cfg.getboolean("keep_calc"):
    Log.write("ACCeL NEB Step3 is deactivated")
    exit()

rxn_atoms = defaultdict(set)
for line in json.loads(cfg["rxn_atoms"]):
    rxn_atoms[(line[0], line[1])].add(int(line[2]))
rxn_atoms: dict[tuple[str], set[int]] = dict(rxn_atoms)
ts_labels: dict[str, str] = json.loads(cfg["ts_labels"])

rea_box = Box(pdir / "201_dft_rea_opt" / "*.log")
rea_box.labeling().set_data("diagram_role", "reactant")

ts_box = Box(pdir / "203_dft_ts_opt" / "*.log")
ts_box.set_data("diagram_role", "ts")
for c in ts_box.get():
    c.label = ts_labels[c.name]
    c.data["diagram_connection"] = [
        "_".join(c.label.split("_")[0:2]),
        "_".join(c.label.split("_")[2:4]),
    ]

pro_box = Box(pdir / "202_dft_pro_opt" / "*.log")
pro_box.labeling().set_data("diagram_role", "product")

all_box = Box(rea_box.get() + ts_box.get() + pro_box.get())

if len(all_box) == 0:
    Log.write("No DFT output")
    config["DEFAULT"]["keep_calc"] = "False"
    exit()

all_box.check_end().read_atoms().read_energy()

if cfg.getboolean("sp_only"):
    ts_box.search("125_emp_ts_freq", suffix=".out")
    rea_box.search("114_emp_reapro_freq", suffix=".out")
    pro_box.search("114_emp_reapro_freq", suffix=".out")
    all_box.read_thermal()
    all_box.calc_energy(["g16_scf", "xtb_g_rrho_contrib"], unit=Units.hartree)
    ts_box.search("203_dft_ts_opt", suffix=".log")
    rea_box.search("201_dft_rea_opt", suffix=".log")
    pro_box.search("202_dft_pro_opt", suffix=".log")
else:
    rea_box.check_freq()
    ts_box.check_freq(im=1)
    pro_box.check_freq()

    all_box.read_thermal().calc_free_energy()
    ts_box.plugin.gau.read_vibration()

    for label, cs in ts_box.get().labels.items():
        rxn_tuple = ("_".join(label.split("_")[0:2]), "_".join(label.split("_")[2:4]))
        Box().bind(cs).plugin.gau.check_bonding(
            list(rxn_atoms[rxn_tuple]),
            additional_valid_range=additional_valid_range,
            acceptable_invalid_atoms=acceptable_invalid_atoms,
        )

all_box.calc_bonds(cfg.getfloat("cov_scaling")).rmsd_limit(cfg.getfloat("rmsd_limit"))
all_box.calc_distribution()

stem_name = all_box.get().get().name.split("_")[0]

all_box.copy_files(pdir / "300_dft_all_log")
all_box.write_xyz(pdir / "301_dft_all_xyz", link=False)
all_box.export_data(pdir / "301_dft_all_xyz" / stem_name)


(pdir / "310_energy_diagram_dft").mkdir(exist_ok=True)
PlotBox(all_box).diagram(pdir / "310_energy_diagram_dft" / stem_name)


tdir = pdir / "template"
g16_qsh_file = "g16.qsh"
if cfg.getboolean("need_hpc"):
    g16_qsh_file = "g16_hpc.qsh"
with pdir.joinpath("030_solvent").joinpath("solvent_key").open() as f:
    solvent_key = f.readlines()[0]
rqbs_port = cfg.getint("rqbs_port")
if 49152 <= rqbs_port and rqbs_port <= 65535:
    g16_qsh_file = "g16_rqbs.qsh"
if cfg.getboolean("no_chk"):
    g16_qsh_file = g16_qsh_file.split(".")[0] + "_nochk.qsh"


irc_file = "g16_irc.gjf"
irc_qsh_file = g16_qsh_file
reverse_statement = "reverse"
if cfg.getboolean("sp_only"):
    irc_file = "orca_irc.inp"
    irc_qsh_file = "orca_irc.qsh"
    reverse_statement = "backward"

ts_box.write_input(
    tdir / irc_file,
    pdir / "403_irc_forward_all",
    arg={"DIRECTION": "forward", "SOLVENT_KEY": solvent_key, "MODEL_OPT": cfg["model_opt"]},
    link=False,
)
ts_box.write_input(tdir / irc_qsh_file, pdir / "403_irc_forward_all", link=False, arg={"PORT": rqbs_port})
ts_box.write_input(
    tdir / irc_file,
    pdir / "404_irc_reverse_all",
    arg={"DIRECTION": reverse_statement, "SOLVENT_KEY": solvent_key, "MODEL_OPT": cfg["model_opt"]},
    link=False,
)
ts_box.write_input(tdir / irc_qsh_file, pdir / "404_irc_reverse_all", link=False, arg={"PORT": rqbs_port})

all_box.only_minimum()
all_box.export_data(pdir / "310_energy_diagram_dft" / stem_name)
all_box.copy_files(pdir / "311_dft_min_log")
all_box.write_xyz(pdir / "312_dft_min_xyz")

ts_box.write_input(
    tdir / irc_file,
    pdir / "401_irc_forward",
    arg={"DIRECTION": "forward", "SOLVENT_KEY": solvent_key, "MODEL_OPT": cfg["model_opt"]},
)
ts_box.write_input(tdir / irc_qsh_file, pdir / "401_irc_forward", arg={"PORT": rqbs_port})
ts_box.write_input(
    tdir / irc_file,
    pdir / "402_irc_reverse",
    arg={"DIRECTION": reverse_statement, "SOLVENT_KEY": solvent_key, "MODEL_OPT": cfg["model_opt"]},
)
ts_box.write_input(tdir / irc_qsh_file, pdir / "402_irc_reverse", arg={"PORT": rqbs_port})

if cfg.getboolean("need_hpc"):
    shutil.copy(tdir / "qsb", pdir / "401_irc_forward" / "qsb")
    shutil.copy(tdir / "qsb", pdir / "402_irc_reverse" / "qsb")
    shutil.copy(tdir / "qsb", pdir / "403_irc_forward_all" / "qsb")
    shutil.copy(tdir / "qsb", pdir / "404_irc_reverse_all" / "qsb")

with (pdir / "099_config" / "config.ini").open("w") as f:
    config.write(f)

Log.write("ACCel NEB Step3 terminated normally")
