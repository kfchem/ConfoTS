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
Log.write("Starting ACCeL NEB Step2")

if not cfg.getboolean("keep_calc"):
    Log.write("ACCeL NEB Step2 is deactivated")
    exit()

if cfg.getboolean("need_hpc"):
    Log.write("ACCeL NEB Step2 is deactivated")
    exit()

dft_box = Box()
dft_box.add(pdir / "203_dft_ts_opt" / "*.qsh")
if len(dft_box) == 0:
    Log.write("No input for TS calculation: Aborted")
    config["DEFAULT"]["keep_calc"] = "False"
    exit()
dft_box.add(pdir / "201_dft_rea_opt" / "*.qsh")
dft_box.add(pdir / "202_dft_pro_opt" / "*.qsh")
dft_box.run()

for _ in range(3):
    Box(dft_box.plugin.gau.get_resub()).run()

(pdir / "204_dft_data").mkdir(exist_ok=True)
dft_box.export_data(pdir / "204_dft_data" / dft_box.get().get().name.split("_")[0])

with (pdir / "099_config" / "config.ini").open("w") as f:
    config.write(f)

Log.write("ACCel NEB Step2 terminated normally")
