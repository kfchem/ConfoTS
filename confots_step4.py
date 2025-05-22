__version__ = "2022111301"
import configparser
from pathlib import Path

from accel import Box
from accel.util import Log

pdir = Path().cwd()
config = configparser.ConfigParser()
config.read(pdir / "099_config" / "config.ini")
cfg = config["USER"]

Log.console(show=cfg.getboolean("console_show"))
Log.file(pdir.joinpath("_" + Path(__file__).stem).with_suffix(".out"))
Log.write(f"{Path(__file__).absolute()}: version {__version__}")
Log.write("Starting ACCeL NEB Step4")

if not cfg.getboolean("keep_calc"):
    Log.write("ACCeL NEB Step4 is deactivated")
    exit()

if cfg.getboolean("need_hpc"):
    Log.write("ACCeL NEB Step4 is deactivated")
    exit()

if not cfg.getboolean("run_irc"):
    Log.write("ACCeL NEB Step4 is deactivated")
    exit()


dft_box = Box()
if cfg.getboolean("irc_all"):
    dft_box.add(pdir / "403_irc_forward_all" / "*.qsh")
    dft_box.add(pdir / "404_irc_reverse_all" / "*.qsh")
else:
    dft_box.add(pdir / "401_irc_forward" / "*.qsh")
    dft_box.add(pdir / "402_irc_reverse" / "*.qsh")

dft_box.run()

Log.write("ACCel NEB Step4 terminated normally")
