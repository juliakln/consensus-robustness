"""
Expected times analysis of reaching and holding a stable consensus

This script runs PRISM to compute expected times to reach and hold consensus,
reads expected times from PRISM output, and plots the results.


Results are saved in PMC_results/times/. and plots in figures/
"""

import os
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path


# change to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# PRISM location
prism_path = os.environ.get("PRISM_PATH")
if prism_path is None:
    raise EnvironmentError("PRISM_PATH environment variable is not set.")

# Base project folder
project_dir = Path(__file__).resolve().parent.parent
model_path_crossinh_zealots = (project_dir / "models" / "crossinh_zealots.pm").resolve()
model_path_crossinh_contrarians = (project_dir / "models" / "crossinh_contrarians.pm").resolve()
props_path = (project_dir / "models" / "props_times.pctl").resolve()

results_dir = project_dir / "PMC_results" / "times"
figures_dir = project_dir / "figures"
(figures_dir).mkdir(parents=True, exist_ok=True)



"""
Run PRISM model with given parameters and save results to file
"""
def run_prism(model_path, props_path, consts, prop_index, resultfile, maxiters = 1000000):
    resultfile.parent.mkdir(parents=True, exist_ok=True)

    prismcommand = [
        str(prism_path),
        str(model_path),
        str(props_path),
        f"-const {consts}",
        f"-prop {prop_index}",
        f"-sparse -maxiters {maxiters}",
        f"-exportresults {resultfile}"
    ]

    try:
        subprocess.run(" ".join(prismcommand), shell=True, check=True,
                       capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print("\n❌ PRISM ERROR")
        print("Command:", e.cmd)
        print("Exit code:", e.returncode)
        print("Stdout:", e.stdout)
        print("Stderr:", e.stderr)
        raise


"""
Generic function to verify expected reaching/holding times
    prefix: e.g. "reach_zealots" or "hold_contrarians"
    const_template: e.g. "Zx={zealots},Zy={zealots},qx={qx},qy={qy}"
"""
def verify_expected_times(model_path, props_path, prefix, rates, const_template):
    ratex, ratey = rates
    values_range = np.arange(1, 46, 1)  # from 2 to 90 disruptive individuals

    for value in values_range:
        resultfile = (results_dir / f"{prefix}_{rate_label(rates)}" / f"prob_{value}.txt").resolve()

        if not resultfile.exists():
            consts = const_template(value, ratex, ratey)
            prop_index = 1 if "reach" in prefix else 2
            run_prism(model_path, props_path, consts, prop_index, resultfile)


"""
Read PRISM result files and return sorted pairs of values and expected times
"""
def read_expected_times(prefix, rates):
    folder = results_dir / f"{prefix}_{rate_label(rates)}"
    if not folder.exists():
        raise FileNotFoundError(f"Results folder '{folder}' does not exist.")
    
    times = {}
    for file in folder.glob("prob_*.txt"):
        value = int(file.stem.split("_")[1])
        with open(file, 'r') as f:
            for last_line in f:
                pass
            try:
                t = float(last_line)
            except ValueError:
                t = np.nan
            times[value] = t

    values = sorted(times.keys())
    outputs = [times[v] for v in values]
    return values, outputs


"""
Plot expected times for zealots and contrarians
"""
def plot_times(rates, values_z, times_z, values_c, times_c, title, filename):
    xvals = [2 * v for v in values_z]  # total disruptive individuals

    fig = plt.figure(figsize=(6,6))
    plt.plot(xvals, np.log(times_z), 'b', linewidth = 1.5, label = 'Zealots')
    plt.plot(xvals, np.log(times_c), 'r', linewidth = 1.5, label = 'Contrarians')

    plt.xlim(0,100)
    plt.ylim(-3, 10)
    plt.xlabel("Amount of disruptive individuals")
    plt.ylabel("Expected time in log scale")
    plt.title(f"{title}: {str(rates)} rates")
    plt.legend()
    fig.savefig(figures_dir / filename)
    plt.close()


"""
Helper function to create constant strings for zealots and contrarians
"""
def const_zealots(value, ratex, ratey):
    return f"Zx={value},Zy={value},qx={ratex},qy={ratey}"

def const_contrarians(value, ratex, ratey):
    return f"cHalf={value},qx={ratex},qy={ratey}"

"""
Helper function to create rate label strings
"""
def rate_label(rates):
    qx, qy = rates
    return f"{str(qx).replace('.', 'p')}x_{str(qy).replace('.', 'p')}y"




def main():

    rate_list = [
    (0.01, 0.01),
    (0.0101, 0.0099),
    (0.0105, 0.0095),
    (0.011, 0.009),
    (0.014, 0.006)
    ]

    for rates in rate_list:
        verify_expected_times(model_path_crossinh_zealots, props_path, "hold_zealots", rates, const_zealots)
        verify_expected_times(model_path_crossinh_contrarians, props_path, "hold_contrarians", rates, const_contrarians)
        values_z, times_z = read_expected_times("hold_zealots", rates)
        values_c, times_c = read_expected_times("hold_contrarians", rates)
        plot_times(rates, values_z, times_z, values_c, times_c, "Expected time to hold consensus", f"holdingtimes_{rate_label(rates)}.png")



if __name__ == "__main__":
    main()
