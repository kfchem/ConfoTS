__version__ = "20230500701"
import datetime
import statistics
from pathlib import Path

atom_symbols = [
    "",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Uub",
    "Uut",
    "Uuq",
    "Uup",
    "Uuh",
    "Uus",
    "Uuo",
]

if Path("_batch_log").exists():
    log_path = Path("_batch_log").joinpath("atm.log")
else:
    log_path = Path("atm_log.out")
with log_path.open("a") as f:
    f.write("{} {}\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "atm.py was invoked."))


def log_print(comment, console_print=True):
    with log_path.open("a") as f:
        f.write("{} {}\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), comment))
    if console_print:
        print(comment)


def yes_no(comment, default_value="y"):
    if default_value is True:
        default_value = "y"
    if default_value is False:
        default_value = "n"
    while True:
        _ans = default_value
        _txt = "{} y/n (default = {}): ".format(comment, default_value)
        _your_ans = input(_txt)
        log_print(_txt + _your_ans, console_print=False)
        if _your_ans != "":
            _ans = _your_ans
        if _ans.lower() in ["y", "yes", "ye"]:
            return True
        elif _ans.lower() in ["n", "no"]:
            return False


def gen_gjf_from_final_stdorient(_p: Path):
    with _p.open() as f:
        _ls = f.readlines()
    geom_no = 0
    for i, _l in enumerate(_ls):
        if "Standard orientation:" in _l:
            geom_no = i
    if geom_no == 0:
        log_print("Input File Error: " + _p.name)
        return False
    geom = []
    for _l in _ls[(geom_no + 5) :]:
        if "---------------------------------------------------------------------" in _l:
            break
        geom.append(_l.split())
    geom = ["{1:<2} {0[3]:>11} {0[4]:>11} {0[5]:>11}\n".format(_g, atom_symbols[int(_g[1])]) for _g in geom]

    with _p.with_suffix(".gjf").open() as f:
        _gjf_ls = f.readlines()
    blank_line_index = [i for i, _gjf_l in enumerate(_gjf_ls) if _gjf_l == "\n"]
    for line in _gjf_ls:
        if "qst3" in line.lower():
            before_geom = blank_line_index[5] + 2
            after_geom = blank_line_index[6]
            break
    else:
        before_geom = blank_line_index[1] + 2
        after_geom = blank_line_index[2]
    _out_lines = _gjf_ls[0:before_geom]
    _out_lines.extend(geom)
    _out_lines.extend(_gjf_ls[after_geom:])

    _suffix = "_"
    for _ext in [".gjf", ".log", ".chk"]:
        while _p.with_suffix(_ext + _suffix).exists():
            if not yes_no("{} is found. continue?".format(_p.with_suffix(_ext + _suffix).name), default_value="n"):
                log_print("{} will not be in the que.".format(_p.stem))
                return False
            _suffix = _suffix + "_"

    for _ext in [".gjf", ".log", ".chk"]:
        if _p.with_suffix(_ext).exists():
            _p.with_suffix(_ext).rename(_p.with_suffix(_ext + _suffix))
            log_print("{} was renamed as {}.".format(_p.with_suffix(_ext).name, _p.with_suffix(_ext + _suffix).name))

    with _p.with_suffix(".gjf").open("w") as f:
        f.writelines(_out_lines)
    log_print("{} was generated.".format(_p.with_suffix(".gjf").name))
    return True


no_log_err = {
    "name": "no log file error",
    "statement": None,
    "default_resubmission": True,
    "from_final_geom": False,
    "paths": list(),
}

unknown_err = {
    "name": "unknown error",
    "statement": None,
    "default_resubmission": True,
    "from_final_geom": False,
    "paths": list(),
}

imaginary_freq_err = {
    "name": "inappropriate imaginary frequency",
    "statement": "imaginary frequencies (negative Signs)",
    "default_resubmission": True,
    "from_final_geom": True,
    "paths": list(),
}

errs = dict()
errs["float_exp"] = {
    "name": "floating point exception",
    "statement": "Error: floating point exception, integer divide by zero",
    "default_resubmission": True,
    "from_final_geom": False,
}
errs["alloc_mem"] = {
    "name": "memory allocation error",
    "statement": "could not allocate memory.: Cannot allocate memory",
    "default_resubmission": True,
    "from_final_geom": False,
}
errs["pcmmku_fail"] = {
    "name": "PCMMkU failure",
    "statement": "Inv3 failed in PCMMkU.",
    "default_resubmission": False,
    "from_final_geom": True,
}
errs["scf_conv"] = {
    "name": "SCF convergence failure",
    "statement": "SCF has not converged.",
    "default_resubmission": False,
    "from_final_geom": True,
}
errs["frombx_problem"] = {
    "name": "FormBX problem",
    "statement": "FormBX had a problem.",
    "default_resubmission": True,
    "from_final_geom": True,
}
errs["link9999"] = {
    "name": "link9999 error",
    "statement": "Error termination request processed by link 9999.",
    "default_resubmission": True,
    "from_final_geom": True,
}

for _key in errs.keys():
    errs[_key]["paths"] = list()

for gjf_p in sorted(Path.cwd().glob("*.gjf")):
    log_print("{} is processing.".format(gjf_p.name))
    with gjf_p.open("r") as f:
        gjf_ls = f.readlines()
    tot_jobs = 1
    check_freq = False
    is_ts_freq = False
    for i, gjf_l in enumerate(gjf_ls):
        if "--link1--" in gjf_l.lower():
            tot_jobs += 1
        if "freq" in gjf_l.lower() and "opt" in gjf_l.lower():
            tot_jobs += 1
        if "ts" in gjf_l.lower() and check_freq is False:
            is_ts_freq = True
        if "qst" in gjf_l.lower() and check_freq is False:
            is_ts_freq = True
        if "freq" in gjf_l.lower():
            check_freq = True

    log_p = gjf_p.with_suffix(".log")
    if not log_p.exists():
        no_log_err["paths"].append(gjf_p)
        continue

    norm_term = 0
    err_flag = False

    with log_p.open("r") as f:
        log_ls = f.readlines()

    for log_l in reversed(log_ls):
        if norm_term == 0:
            for _e in errs.values():
                if _e["statement"] in log_l:
                    _e["paths"].append(log_p)
                    err_flag = True
                    log_print("{} was detected.".format(_e["name"]))
                    break
            else:
                if "Normal termination" in log_l:
                    norm_term += 1
                continue
            break
        if "Normal termination" in log_l:
            norm_term += 1

    if not err_flag:
        if norm_term == tot_jobs:
            log_print("{} terminated normally.".format(log_p.name))
            if check_freq:
                imlines = [log_l for log_l in log_ls if imaginary_freq_err["statement"] in log_l]
                if imlines == []:
                    num_freq = 0
                else:
                    num_freq = int(imlines[-1].split()[1])
                if (is_ts_freq and num_freq != 1) or (not is_ts_freq and num_freq != 0):
                    imaginary_freq_err["paths"].append(log_p)
                    err_flag = True
                    log_print("{} was detected: {} imaginary freq.".format(imaginary_freq_err["name"], num_freq))
        else:
            unknown_err["paths"].append(log_p)
            log_print("unknown termination was detected.")

print()

que_qsub_ps = []

errs["no_log"] = no_log_err
errs["imaginary_freq"] = imaginary_freq_err
errs["link9999_wo_vib"] = {
    "name": "link9999 (not vibrational) error",
    "statement": None,
    "default_resubmission": False,
    "from_final_geom": True,
    "paths": list(),
}
errs["link9999_w_vib"] = {
    "name": "link9999 (vibrational ending) error",
    "statement": None,
    "default_resubmission": True,
    "from_final_geom": True,
    "paths": list(),
}
for log_p in errs["link9999"]["paths"]:
    with log_p.open("r") as f:
        log_ls = f.readlines()
    rms_force_vs = []
    for i, log_l in enumerate(log_ls):
        if "RMS     Force" in log_l:
            try:
                rms_force_vs.append(float(log_l.split()[2]))
            except ValueError:
                log_print("ValueError during collecting RMS forces:{}".format(log_l))
    if len(rms_force_vs) < 12:
        log_print(f"Exception in handring link9999 error: {log_p.name}: short iteration")
        continue
    if (statistics.stdev(rms_force_vs[-12::2]) < 0.000001) or (statistics.stdev(rms_force_vs[-12::3]) < 0.000001):
        errs["link9999_w_vib"]["paths"].append(log_p)
    else:
        errs["link9999_wo_vib"]["paths"].append(log_p)
errs.pop("link9999")

for _e in errs.values():
    if len(_e["paths"]) != 0:
        log_print("{} log files terminated by {}.".format(len(_e["paths"]), _e["name"]))
        for _p in _e["paths"]:
            log_print("   {}".format(_p.name))
        if _e["from_final_geom"]:
            if yes_no("Submit again from final coordinate?", default_value=_e["default_resubmission"]):
                for _p in _e["paths"]:
                    if gen_gjf_from_final_stdorient(_p):
                        que_qsub_ps.append(_p.with_suffix(".qsh"))
        else:
            if yes_no("Submit again?", default_value=_e["default_resubmission"]):
                que_qsub_ps.extend([_p.with_suffix(".qsh") for _p in _e["paths"]])
        print()

if len(unknown_err["paths"]) != 0:
    log_print("{} log files with unknown termination or under calculation.".format(len(unknown_err["paths"])))
    for _p in unknown_err["paths"]:
        log_print("   {}".format(_p.name))
    if yes_no("Submit again?", default_value="n"):
        if yes_no("for each?", default_value="y"):
            for _p in unknown_err["paths"]:
                if yes_no("Submit again?: {}".format(_p.name), default_value="y"):
                    if yes_no("from final coordinate?", default_value="y"):
                        if gen_gjf_from_final_stdorient(_p):
                            que_qsub_ps.append(_p.with_suffix(".qsh"))
                    else:
                        que_qsub_ps.append(_p.with_suffix(".qsh"))
        else:
            que_qsub_ps.extend([_p.with_suffix(".qsh") for _p in unknown_err["paths"]])
    print()

if len(que_qsub_ps) != 0:
    with Path("atmqsubs.sh").open("w") as f:
        f.write("#!/bin/sh\n\n")
        if Path("_batch_log").exists():
            f.write("eval date >> ./_batch_log/atmqsubs.log\n\n")
            sh_ls = [
                "eval qsub -e ./_batch_log/ -o ./_batch_log/ {} >> ./_batch_log/atmqsubs.log\n".format(que_qsub_p.name)
                for que_qsub_p in que_qsub_ps
            ]
        else:
            f.write("eval date >> atmqsubs.out\n\n")
            sh_ls = ["eval qsub {} >> atmqsubs.out\n".format(que_qsub_p.name) for que_qsub_p in que_qsub_ps]
        f.writelines(sh_ls)
        log_print("A script file was created. Please excute the following command.")
        log_print("sh atmqsubs.sh")
        with log_path.open("a") as f:
            sh_ls = ["the contents are as follows\n"] + sh_ls
            f.writelines(["{} {}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), _l) for _l in sh_ls])
else:
    log_print("There is no resubmission.")
