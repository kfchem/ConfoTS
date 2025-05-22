__version__ = "2025052201"
import configparser
from pathlib import Path

from accel import Box
from accel.util import Log
from acceltools.series import SeriesBox

pdir = Path().cwd()
config = configparser.ConfigParser()
config.read(pdir / "099_config" / "config.ini")
cfg = config["USER"]

tdir = pdir / "template"

Log.console(show=cfg.getboolean("console_show"))
Log.file(pdir.joinpath("_" + Path(__file__).stem).with_suffix(".out"))
Log.write(f"{Path(__file__).absolute()}: version {__version__}")
Log.write("Starting ConfoTS Step5")

if not cfg.getboolean("keep_calc"):
    Log.write("ConfoTS Step5 is deactivated")
    exit()

if not cfg.getboolean("run_irc"):
    Log.write("ConfoTS Step5 is deactivated")
    exit()

with pdir.joinpath("030_solvent").joinpath("solvent_key").open() as f:
    solvent_key = f.readlines()[0]

ircbox = Box()
if cfg.getboolean("irc_all"):
    ircbox.add(pdir / "403_irc_forward_all" / "*.log")
    ircbox.add(pdir / "404_irc_reverse_all" / "*.log")
else:
    ircbox.add(pdir / "401_irc_forward" / "*.log")
    ircbox.add(pdir / "402_irc_reverse" / "*.log")

ircbox.check_end()

box = Box(ircbox.plugin.gau.get_irc())
box.calc_bonds(cfg.getfloat("cov_scaling"))

for c in box.get():
    if "R_TS" in c.name:
        c.deactivate("duplicated TS")

for cs in box.get().labels.values():
    min_rxn_coord = min([c.data["rxn_coordinate"] for c in cs])
    max_rxn_coord = max([c.data["rxn_coordinate"] for c in cs])
    terminal = Box(
        [cs.has_data("rxn_coordinate", min_rxn_coord).get(), cs.has_data("rxn_coordinate", max_rxn_coord).get()]
    )
    terminal.write_input(
        tdir / "g16_opt.gjf",
        pdir / "502_terminal_opt",
        link=False,
        arg={"SOLVENT_KEY": solvent_key, "MODEL_OPT": cfg["model_opt"], "MODEL_SP": cfg["model_sp"]},
    )
    terminal.write_input(tdir / "g16.qsh", pdir / "502_terminal_opt", link=False)

SeriesBox(box).write_trjxyz(pdir / "501_irc_trj", order_key="rxn_coordinate")
box.export_data(pdir / "501_irc_trj" / pdir.name)

Log.write("ConfoTS Step5 terminated normally")
