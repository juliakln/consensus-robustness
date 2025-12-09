"""
Helper functions 
- Write models in RML and properties in BLTL for Plasmalab
    - Cross-inhibition model with zealots, contrarians, or both
    - Property of reaching stable consensus or switching consensus
- Run Plasmalab and compute probability of property, save in txt file
- Read txt data
- Plot probabilities over number of stubborn individuals in the system for different settings
"""

import os
import numpy as np
import subprocess
import time
import random
import matplotlib.pyplot as plt
from collections import defaultdict

from models import *
from properties import *

# change to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# define file for model and property
model = "../../models/consensus_model.rml"
property = "../../models/consensus.bltl"

# define colours and linetypes for plotting
#colours = ['b', 'b', 'b', '#ffaa00', '#ffaa00', 'r', 'r', 'g', 'g', 'm', 'm', 'c', 'c']
colours = ['b','#ffaa00', '#ffaa00', 'b', 'b', 'g', 'g', 'r', 'r', 'g', 'g', 'm', 'm', 'c', 'c']
linetypes = ['-', '-', '--', '-', '--', '-', '--', '-', '--', '-', '--', '-']

colours_switch = ['b', 'b', 'b', 'b', 'b', 'r', 'r', 'r', 'r']
linetypes_switch = ['-', '-', '--', ':', '-.', '-', '--', ':', '-.','-', '--', ':', '-.']

colours_groups = ['k', 'g', 'r', 'b', 'c', 'm', 'y']







"""
Run Plasmalab with commandline and compute probability of stable/switch property using Monte Carlo
    stubborn: z (zealots), c (contrarians)
    N: total population size
    majority: more than m% of population commits to same decision
    distance: difference of at least d between majority and those favouring opposite decision
    transient: consensus is achieved within t minutes from start of dynamics
    holding: group maintains consensus for at least h minutes
    range: range of stubborn individuals to compute robustness
    filename: where to save the result files
    samples: how many samples to use for Monte Carlo

    Return: dictionary with keys = #stubborns and values = probabilities, read from result files
"""
def stableconsensus(stubborn, N, majority, distance, transient, holding, range, filename, samples):
    dir_con = '../inference_results/' + filename + '_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in range:
        result = "./plasmares_" + str(int(s)) + ".txt"
        if not os.path.exists(result):
            write_model(stubborn, N, int(s))
            write_property_stableconsensus(N, stubborn, majority, distance, transient, holding)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
            pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data(os.getcwd())
    os.chdir('../')

    return probs

def stableconsensus_voter(stubborn, N, majority, distance, transient, holding, range, filename, samples):
    dir_con = '../inference_results/' + filename + '_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in range:
        result = "./plasmares_" + str(int(s)) + ".txt"
        if not os.path.exists(result):
            write_votermodel(stubborn, N, int(s))
            write_property_stableconsensus(N, stubborn, majority, distance, transient, holding)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
            pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data(os.getcwd())
    os.chdir('../')

    return probs

def stableconsensus_voter_asym_x(stubborn, N, majority, distance, transient, holding, z_range, filename, samples, ratex, ratey):
    dir_con = '../inference_results/' + filename + '_x_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(int(s)) + ".txt"
        if not os.path.exists(result):
            write_votermodel_asym(stubborn, N, int(s), ratex, ratey)
            write_property_stableconsensus_asym_x(N, stubborn, majority, distance, transient, holding)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
            pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data(os.getcwd())
    os.chdir('../')

    return probs


def stableconsensus_voter_asym_y(stubborn, N, majority, distance, transient, holding, range, filename, samples, ratex, ratey):
    dir_con = '../inference_results/' + filename + '_y_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in range:
        result = "./plasmares_" + str(int(s)) + ".txt"
        if not os.path.exists(result):
            write_votermodel_asym(stubborn, N, int(s), ratex, ratey)
            write_property_stableconsensus_asym_y(N, stubborn, majority, distance, transient, holding)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
            pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data(os.getcwd())
    os.chdir('../')

    return probs



def stableconsensus_ci_asym_x(stubborn, N, majority, distance, transient, holding, z_range, filename, samples, ratex, ratey):
    dir_con = '../inference_results/' + filename + '_x_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(int(s)) + ".txt"
        if not os.path.exists(result):
            write_ci_asym(stubborn, N, int(s), ratex, ratey)
            write_property_stableconsensus_asym_x(N, stubborn, majority, distance, transient, holding)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
            pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data(os.getcwd())
    os.chdir('../')

    return probs



def stableconsensus_combined_asym_x(stubborn, N, majority, distance, transient, holding, z_range, filename, samples, ratex, ratey):
    dir_con = '../inference_results/' + filename + '_x_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(int(s)) + ".txt"
        if not os.path.exists(result):
            write_combinedmodel_asym(stubborn, N, int(s), ratex, ratey)
            write_property_stableconsensus_asym_x(N, stubborn, majority, distance, transient, holding)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
            pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data(os.getcwd())
    os.chdir('../')

    return probs


def stableconsensus_ci_asym_y(stubborn, N, majority, distance, transient, holding, range, filename, samples, ratex, ratey):
    dir_con = '../inference_results/' + filename + '_y_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in range:
        result = "./plasmares_" + str(int(s)) + ".txt"
        if not os.path.exists(result):
            write_ci_asym(stubborn, N, int(s), ratex, ratey)
            write_property_stableconsensus_asym_y(N, stubborn, majority, distance, transient, holding)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
            pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data(os.getcwd())
    os.chdir('../')

    return probs





def switchconsensus(stubborn, N, majority, distance, transient, holding, range, filename, samples):
    dir_con = '../inference_results/' + filename + '_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in range:
        result = "./plasmares_" + str(int(s)) + ".txt"
        if not os.path.exists(result):
            write_model(stubborn, N, int(s))
            write_property_switching(stubborn, distance, transient, holding)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
            pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data(os.getcwd())
    os.chdir('../')

    return probs



def recovery_voter(stubborn, N, majority, distance, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recovery/' + filename + '_('+str(majority)+','+str(distance)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_votermodel_asym(stubborn, N, int(s), 1.05, 0.95)
            write_property_recovery(N, stubborn, majority, distance, G, rec_time)

            success = False
            max_retries = 5
            for attempt in range(1, max_retries + 1):
                seed = random.randint(0, 999999)
                print(f"▶️  Run {s}: Versuch {attempt}/{max_retries}")
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result

                try:
                    proc = subprocess.run(
                        pcommand,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        timeout=100  # 2 Minuten Timeout
                    )
                    output = proc.stdout
                    #print(output)

                    # Nach IndexOutOfBoundsException suchen
                    if "IndexOutOfBoundsException" in output or "Exception" in output:
                        print(f"⚠️  PlasmaLab IndexOutOfBoundsException (Versuch {attempt}) – retrying...")
                        time.sleep(1)
                        continue

                    # Alles ok
                    success = True
                    break

                except subprocess.TimeoutExpired:
                    print(f"⚠️  PlasmaLab TimeoutExpired bei Versuch {attempt} – retrying...")
                    continue

                except Exception as e:
                    print(f"⚠️  Sonstiger Fehler bei Versuch {attempt}: {e}")
                    time.sleep(2)

            if not success:
                print(f"❌ PlasmaLab konnte für s={s} nach {max_retries} Versuchen nicht ausgeführt werden. Überspringe diesen Wert.")
                with open(result, "w") as f:
                    f.write("-1.0\n")
                continue



    probs = read_data(os.getcwd())
    os.chdir('../')

    return probs


def recovery_ci(stubborn, N, majority, distance, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recovery/' + filename + '_('+str(majority)+','+str(distance)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_ci_asym(stubborn, N, int(s), 1.05, 0.95)
            write_property_recovery(N, stubborn, majority, distance, G, rec_time)

            success = False
            max_retries = 5
            for attempt in range(1, max_retries + 1):
                seed = random.randint(0, 999999)
                print(f"▶️  Run {s}: Versuch {attempt}/{max_retries}")
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result

                try:
                    proc = subprocess.run(
                        pcommand,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        timeout=100  # 2 Minuten Timeout
                    )
                    output = proc.stdout
                    #print(output)

                    # Nach IndexOutOfBoundsException suchen
                    if "IndexOutOfBoundsException" in output or "Exception" in output:
                        print(f"⚠️  PlasmaLab IndexOutOfBoundsException (Versuch {attempt}) – retrying...")
                        time.sleep(1)
                        continue

                    # Alles ok
                    success = True
                    break

                except subprocess.TimeoutExpired:
                    print(f"⚠️  PlasmaLab TimeoutExpired bei Versuch {attempt} – retrying...")
                    continue

                except Exception as e:
                    print(f"⚠️  Sonstiger Fehler bei Versuch {attempt}: {e}")
                    time.sleep(2)

            if not success:
                print(f"❌ PlasmaLab konnte für s={s} nach {max_retries} Versuchen nicht ausgeführt werden. Überspringe diesen Wert.")
                with open(result, "w") as f:
                    f.write("-1.0\n")
                continue



    probs = read_data(os.getcwd())
    os.chdir('../')

    return probs



def recoveryandconsensus_voter(stubborn, N, majority, distance, F, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recoveryandconsensus/' + filename + '_('+str(F)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_votermodel_asym(stubborn, N, int(s), 1.05, 0.95)
            write_property_recoveryandconsensus(N, stubborn, majority, distance, F, G, rec_time)

            success = False
            max_retries = 5
            for attempt in range(1, max_retries + 1):
                seed = random.randint(0, 999999)
                print(f"▶️  Run {s}: Versuch {attempt}/{max_retries}")
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result

                try:
                    proc = subprocess.run(
                        pcommand,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        timeout=100  # 2 Minuten Timeout
                    )
                    output = proc.stdout
                    #print(output)

                    # Nach IndexOutOfBoundsException suchen
                    if "IndexOutOfBoundsException" in output or "Exception" in output:
                        print(f"⚠️  PlasmaLab IndexOutOfBoundsException (Versuch {attempt}) – retrying...")
                        time.sleep(1)
                        continue

                    # Alles ok
                    success = True
                    break

                except subprocess.TimeoutExpired:
                    print(f"⚠️  PlasmaLab TimeoutExpired bei Versuch {attempt} – retrying...")
                    continue

                except Exception as e:
                    print(f"⚠️  Sonstiger Fehler bei Versuch {attempt}: {e}")
                    time.sleep(2)

            if not success:
                print(f"❌ PlasmaLab konnte für s={s} nach {max_retries} Versuchen nicht ausgeführt werden. Überspringe diesen Wert.")
                with open(result, "w") as f:
                    f.write("-1.0\n")
                continue



    probs = read_data_recovery(os.getcwd())
    os.chdir('../')

    return probs


def recoveryandconsensus_ci(stubborn, N, majority, distance, F, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recoveryandconsensus/' + filename + '_('+str(F)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_ci_asym(stubborn, N, int(s), 0.105, 0.095)
            write_property_recoveryandconsensus(N, stubborn, majority, distance, F, G, rec_time)

            success = False
            max_retries = 5
            for attempt in range(1, max_retries + 1):
                seed = random.randint(0, 999999)
                print(f"▶️  Run {s}: Versuch {attempt}/{max_retries}")
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result

                try:
                    proc = subprocess.run(
                        pcommand,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        timeout=100  # 2 Minuten Timeout
                    )
                    output = proc.stdout
                    #print(output)

                    # Nach IndexOutOfBoundsException suchen
                    if "IndexOutOfBoundsException" in output or "Exception" in output:
                        print(f"⚠️  PlasmaLab IndexOutOfBoundsException (Versuch {attempt}) – retrying...")
                        time.sleep(1)
                        continue

                    # Alles ok
                    success = True
                    break

                except subprocess.TimeoutExpired:
                    print(f"⚠️  PlasmaLab TimeoutExpired bei Versuch {attempt} – retrying...")
                    continue

                except Exception as e:
                    print(f"⚠️  Sonstiger Fehler bei Versuch {attempt}: {e}")
                    time.sleep(2)

            if not success:
                print(f"❌ PlasmaLab konnte für s={s} nach {max_retries} Versuchen nicht ausgeführt werden. Überspringe diesen Wert.")
                with open(result, "w") as f:
                    f.write("-1.0\n")
                continue



    probs = read_data_recovery(os.getcwd())
    os.chdir('../')

    return probs



def recoveryandconsensus_combined(stubborn, N, majority, distance, F, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recoveryandconsensus/' + filename + '_('+str(F)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_combinedmodel_asym(stubborn, N, int(s), 1.05, 0.95)
            write_property_recoveryandconsensus(N, stubborn, majority, distance, F, G, rec_time)

            success = False
            max_retries = 5
            for attempt in range(1, max_retries + 1):
                seed = random.randint(0, 999999)
                print(f"▶️  Run {s}: Versuch {attempt}/{max_retries}")
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result

                try:
                    proc = subprocess.run(
                        pcommand,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        timeout=100  # 2 Minuten Timeout
                    )
                    output = proc.stdout
                    #print(output)

                    # Nach IndexOutOfBoundsException suchen
                    if "IndexOutOfBoundsException" in output or "Exception" in output:
                        print(f"⚠️  PlasmaLab IndexOutOfBoundsException (Versuch {attempt}) – retrying...")
                        time.sleep(1)
                        continue

                    # Alles ok
                    success = True
                    break

                except subprocess.TimeoutExpired:
                    print(f"⚠️  PlasmaLab TimeoutExpired bei Versuch {attempt} – retrying...")
                    continue

                except Exception as e:
                    print(f"⚠️  Sonstiger Fehler bei Versuch {attempt}: {e}")
                    time.sleep(2)

            if not success:
                print(f"❌ PlasmaLab konnte für s={s} nach {max_retries} Versuchen nicht ausgeführt werden. Überspringe diesen Wert.")
                with open(result, "w") as f:
                    f.write("-1.0\n")
                continue



    probs = read_data_recovery(os.getcwd())
    os.chdir('../')

    return probs




def recoveryandconsensusgood_voter(stubborn, N, majority, distance, F, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recoveryandconsensusgood/' + filename + '_('+str(F)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_votermodel_asym(stubborn, N, int(s), 1.05, 0.95)
            write_property_recoveryandconsensusgood(N, stubborn, majority, distance, F, G, rec_time)

            success = False
            max_retries = 5
            for attempt in range(1, max_retries + 1):
                seed = random.randint(0, 999999)
                print(f"▶️  Run {s}: Versuch {attempt}/{max_retries}")
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result

                try:
                    proc = subprocess.run(
                        pcommand,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        timeout=100  # 2 Minuten Timeout
                    )
                    output = proc.stdout
                    #print(output)

                    # Nach IndexOutOfBoundsException suchen
                    if "IndexOutOfBoundsException" in output or "Exception" in output:
                        print(f"⚠️  PlasmaLab IndexOutOfBoundsException (Versuch {attempt}) – retrying...")
                        time.sleep(1)
                        continue

                    # Alles ok
                    success = True
                    break

                except subprocess.TimeoutExpired:
                    print(f"⚠️  PlasmaLab TimeoutExpired bei Versuch {attempt} – retrying...")
                    continue

                except Exception as e:
                    print(f"⚠️  Sonstiger Fehler bei Versuch {attempt}: {e}")
                    time.sleep(2)

            if not success:
                print(f"❌ PlasmaLab konnte für s={s} nach {max_retries} Versuchen nicht ausgeführt werden. Überspringe diesen Wert.")
                with open(result, "w") as f:
                    f.write("-1.0\n")
                continue



    probs = read_data_recovery(os.getcwd())
    os.chdir('../')

    return probs


def recoveryandconsensusgood_ci(stubborn, N, majority, distance, F, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recoveryandconsensusgood/' + filename + '_('+str(F)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_ci_asym(stubborn, N, int(s), 1.05, 0.95)
            write_property_recoveryandconsensusgood(N, stubborn, majority, distance, F, G, rec_time)

            success = False
            max_retries = 5
            for attempt in range(1, max_retries + 1):
                seed = random.randint(0, 999999)
                print(f"▶️  Run {s}: Versuch {attempt}/{max_retries}")
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result

                try:
                    proc = subprocess.run(
                        pcommand,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        timeout=100  # 2 Minuten Timeout
                    )
                    output = proc.stdout
                    #print(output)

                    # Nach IndexOutOfBoundsException suchen
                    if "IndexOutOfBoundsException" in output or "Exception" in output:
                        print(f"⚠️  PlasmaLab IndexOutOfBoundsException (Versuch {attempt}) – retrying...")
                        time.sleep(1)
                        continue

                    # Alles ok
                    success = True
                    break

                except subprocess.TimeoutExpired:
                    print(f"⚠️  PlasmaLab TimeoutExpired bei Versuch {attempt} – retrying...")
                    continue

                except Exception as e:
                    print(f"⚠️  Sonstiger Fehler bei Versuch {attempt}: {e}")
                    time.sleep(2)

            if not success:
                print(f"❌ PlasmaLab konnte für s={s} nach {max_retries} Versuchen nicht ausgeführt werden. Überspringe diesen Wert.")
                with open(result, "w") as f:
                    f.write("-1.0\n")
                continue



    probs = read_data_recovery(os.getcwd())
    os.chdir('../')

    return probs



def recoveryandconsensusgood_combined(stubborn, N, majority, distance, F, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recoveryandconsensusgood/' + filename + '_('+str(F)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_combinedmodel_asym(stubborn, N, int(s), 1.05, 0.95)
            write_property_recoveryandconsensusgood(N, stubborn, majority, distance, F, G, rec_time)

            success = False
            max_retries = 5
            for attempt in range(1, max_retries + 1):
                seed = random.randint(0, 999999)
                print(f"▶️  Run {s}: Versuch {attempt}/{max_retries}")
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result

                try:
                    proc = subprocess.run(
                        pcommand,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        timeout=100  # 2 Minuten Timeout
                    )
                    output = proc.stdout
                    #print(output)

                    # Nach IndexOutOfBoundsException suchen
                    if "IndexOutOfBoundsException" in output or "Exception" in output:
                        print(f"⚠️  PlasmaLab IndexOutOfBoundsException (Versuch {attempt}) – retrying...")
                        time.sleep(1)
                        continue

                    # Alles ok
                    success = True
                    break

                except subprocess.TimeoutExpired:
                    print(f"⚠️  PlasmaLab TimeoutExpired bei Versuch {attempt} – retrying...")
                    continue

                except Exception as e:
                    print(f"⚠️  Sonstiger Fehler bei Versuch {attempt}: {e}")
                    time.sleep(2)

            if not success:
                print(f"❌ PlasmaLab konnte für s={s} nach {max_retries} Versuchen nicht ausgeführt werden. Überspringe diesen Wert.")
                with open(result, "w") as f:
                    f.write("-1.0\n")
                continue



    probs = read_data_recovery(os.getcwd())
    os.chdir('../')

    return probs



def recoverygood_voter(stubborn, N, majority, distance, F, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recoverygood/' + filename + '_('+str(F)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_votermodel_asym(stubborn, N, int(s), 1.05, 0.95)
            write_property_recoverygood(N, stubborn, majority, distance, F, G, rec_time)

            seed = random.randint(0, 999999)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result

            proc = subprocess.run(
                pcommand,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=100  # 2 Minuten Timeout
            )
            output = proc.stdout

    probs = read_data_recovery(os.getcwd())
    os.chdir('../')

    return probs


def recoverygood_ci(stubborn, N, majority, distance, F, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recoverygood/' + filename + '_('+str(F)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_ci_asym(stubborn, N, int(s), 1.05, 0.95)
            write_property_recoverygood(N, stubborn, majority, distance, F, G, rec_time)

            seed = random.randint(0, 999999)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result

            proc = subprocess.run(
                        pcommand,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        timeout=100  # 2 Minuten Timeout
                )
            output = proc.stdout


    probs = read_data_recovery(os.getcwd())
    os.chdir('../')

    return probs


def recoverygood_combined(stubborn, N, majority, distance, F, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recoverygood/' + filename + '_('+str(F)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_combinedmodel_asym(stubborn, N, int(s), 1.05, 0.95)
            write_property_recoveryandconsensusgood(N, stubborn, majority, distance, F, G, rec_time)

            seed = random.randint(0, 999999)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result
            proc = subprocess.run(
                        pcommand,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        timeout=100  # 2 Minuten Timeout
                    )
            output = proc.stdout

    probs = read_data_recovery(os.getcwd())
    os.chdir('../')

    return probs


def recoverygood_sc(N, majority, distance, F, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recoverygood/' + filename + '_('+str(F)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_scmodel(N, 1.05, 0.95, 1.0)
            write_property_recoverygood(N, 'z', majority, distance, F, G, rec_time)

            seed = random.randint(0, 999999)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result

            proc = subprocess.run(
                pcommand,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=100  # 2 Minuten Timeout
            )
            output = proc.stdout

    probs = read_data_recovery(os.getcwd())
    os.chdir('../')

    return probs


def recoverygood_cc(N, majority, distance, F, G, rec_time, z_range, filename, samples):
    dir_con = '../inference_results_recoverygood/' + filename + '_('+str(F)+','+str(G)+','+str(int(z_range[0]))+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in z_range:
        result = "./plasmares_" + str(rec_time) + ".txt"
        if not os.path.exists(result):
            write_ccmodel(N, 1.05, 0.95, 1.0)
            write_property_recoverygood(N, 'z', majority, distance, F, G, rec_time)

            seed = random.randint(0, 999999)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -A\"Simulation time\"=1000 -A\"Seed\"="+str(seed)+" -f proba -o " + result

            proc = subprocess.run(
                pcommand,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=100  # 2 Minuten Timeout
            )
            output = proc.stdout

    probs = read_data_recovery(os.getcwd())
    os.chdir('../')

    return probs

"""
Run Plasmalab with commandline and compute probability of stable/switch property with both zealots & contrarians using Monte Carlo
    N: total population size
    majority: more than m% of population commits to same decision
    distance: difference of at least d between majority and those favouring opposite decision
    transient: consensus is achieved within t minutes from start of dynamics
    holding: group maintains consensus for at least h minutes
    range_z: range of zealots to compute robustness
    range_c: range of contrarians to compute robustness
    filename: where to save the result files
    samples: how many samples to use for Monte Carlo

    Return: dictionary with keys = #stubborns and values = probabilities, read from result files
"""
def stableconsensus_both(N, majority, distance, transient, holding, range_z, range_c, filename, samples):
    dir_con = '../inference_results/' + filename + '_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for z in range_z:
        for c in range_c:
            result = "./plasmares_" + str(int(z)) + "_" + str(int(c)) + ".txt"
            if not os.path.exists(result):
                write_model_both(N, int(z), int(c))
                write_property_stableconsensus(N, 'b', majority, distance, transient, holding)
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
                pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data_2dim(dir_con)

    return probs

def switchconsensus_both(N, majority, distance, transient, holding, range_z, range_c, filename, samples):
    dir_con = '../inference_results/' + filename + '_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for z in range_z:
        for c in range_c:
            result = "./plasmares_" + str(int(z)) + "_" + str(int(c)) + ".txt"
            if not os.path.exists(result):
                write_model_both(N, int(z), int(c))
                write_property_switching('b', distance, transient, holding)
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
                pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data_2dim(dir_con)

    return probs



"""
Plot robustness analysis
    stubborn: z (zealots), c (contrarians)
    results: dictionaries of probabilities for different consensus settings
    labels: specification of settings for plotting labels
    zealots: list of zealots for which we computed the probabilities
    ylabel: What to write on y axis
    figname: how to name file

    Output: lineplot png of probabilities over stubborn individuals for different property settings
"""
def plot_results(stubborn, results, labels, N, ylabel, figname):
    fig = plt.figure(figsize=(6,6))
    for i in range(1, len(results)):
#        perc = np.linspace(0,70,len(results[i].keys()))
        p1 = [100 * value / N for value in results[i].keys()]
        plt.plot(p1, results[i].values(), linestyle = linetypes[i], color = colours[i], label = labels[i])
    #perc = np.linspace(0,70,len(results[0].keys()))
    p1 = [100 * value / N for value in results[0].keys()]
    plt.plot(p1, results[0].values(), 'k', linewidth = 1.5, label = 'Baseline')
    #calc_labels = (percentages % 10 == 0)
    #zealots_labels, unique_indices = np.unique(percentages[calc_labels], return_index=True)
    #p_indices = np.where(calc_labels)[0][unique_indices]
    #zealots_ticks = np.array(list(results[0].keys()))[p_indices]
    #zealots_labels = percentages[calc_labels]
    #plt.xticks(zealots_ticks, zealots_labels, rotation='horizontal')
    if stubborn == 'z':
#        plt.xlabel('Percentage of Zealots Z = Zx + Zy')
        plt.xlabel('Amount of zealots as % of the total group')
    elif stubborn == 'c':
        plt.xlabel('Amount of contrarians as % of the total group')
    else:
        raise Exception('Type of stubborn individual not supported.')
    plt.ylabel(ylabel)
    plt.legend()
    fig.savefig('../figures/' + figname + '.png')
    plt.close()   



def plot_results_switch(stubborn, results, labels, N, ylabel, figname):
    fig = plt.figure(figsize=(6,6))
    for i in range(1, len(results)):
#        perc = np.linspace(0,70,len(results[i].keys()))
        p1 = [100 * value / N for value in results[i].keys()]
        plt.plot(p1, results[i].values(), linestyle = linetypes_switch[i], color = colours_switch[i], label = labels[i])
    #perc = np.linspace(0,70,len(results[0].keys()))
    p1 = [100 * value / N for value in results[0].keys()]
    plt.plot(p1, results[0].values(), 'k', linewidth = 1.5, label = 'Baseline')
    if stubborn == 'z':
#        plt.xlabel('Percentage of Zealots Z = Zx + Zy')
        plt.xlabel('Amount of zealots as % of the total group')
    elif stubborn == 'c':
        plt.xlabel('Amount of contrarians as % of the total group')
    else:
        raise Exception('Type of stubborn individual not supported.')
    plt.ylabel(ylabel)
    plt.legend()
    fig.savefig('../figures/' + figname + '.png')
    plt.close()   

""" 
Plot 2dim robustness analysis for model with both zealots & contrarians
    results: dictionaries of probabilities for different consensus settings
    range_z: list of zealots for which we computed the probabilities
    range_c: list of zealots for which we computed the probabilities
    figname: how to name file

    Output: contput plot png of probabilities over both zealots & contrarians for 1 setting
"""
def plot_results_2dim(results, range_z, range_c, figname):
    x, y = np.meshgrid(range_z, range_c)
    p = np.array(list(results.values())).reshape(25,-1)

    fig = plt.figure(figsize=(6,6))
    plt.contourf(x, y, p.T, 50, vmin=0, vmax=1)
    plt.xticks([5,10,15,20,25,30,35,40,45], [10,20,30,40,50,60,70,80,90], rotation='horizontal')
    plt.yticks([5,10,15], [10,20,30], rotation='horizontal')
    plt.xlabel('Amount of zealots Z = Zx + Zy')
    plt.ylabel('Amount of contrarians C = Cx + Cy')
    plt.colorbar()
    plt.clim(0,1)
    fig.savefig('../figures/' + figname + '.png')
    plt.close()

"""
Plot robustness analysis
    stubborn: z (zealots), c (contrarians)
    results: dictionaries of probabilities for different consensus settings
    labels: specification of settings for plotting labels
    zealots: list of zealots for which we computed the probabilities
    ylabel: What to write on y axis
    figname: how to name file

    Output: lineplot png of probabilities over stubborn individuals for different property settings
"""
def plot_results_groups(stubborn, results, labels, percentages, ylabel, figname):
    fig = plt.figure(figsize=(6,6))
    for i in range(0, len(results)):
        plt.plot(results[i].keys(), results[i].values(), color = colours_groups[i], label = labels[i])

    calc_labels = (percentages % 10 == 0)
    zealots_labels, unique_indices = np.unique(percentages[calc_labels], return_index=True)
    p_indices = np.where(calc_labels)[0][unique_indices]
    zealots_ticks = np.array(list(results[-1].keys()))[p_indices]
    #zealots_labels = percentages[calc_labels]
    plt.xticks(zealots_ticks, zealots_labels, rotation='horizontal')
    if stubborn == 'z':
#        plt.xlabel('Percentage of Zealots Z = Zx + Zy')
        plt.xlabel('Proportion of Zealots as % of the total group')
    elif stubborn == 'c':
        plt.xlabel('Proportion of Contrarians as % of the total group')
    else:
        raise Exception('Type of stubborn individual not supported.')
    plt.ylabel(ylabel)
    plt.legend()
    fig.savefig('../figures/' + figname + '.png')
    plt.close()   


"""
Read txt files containing probability of satisfying property
    dir: current directory of txt files

    Output: dictionary of probabilities, sorted by #stubborn individuals
"""
def read_data(dir = os.getcwd()):
    probabilities = {}

    for dirpath, dirs, files in os.walk(dir):
        for file in files:
            if file.startswith("plasmares_"):
                zealots = int((file.split("_")[1]).split(".")[0])
                with open(os.path.join(dirpath, file), 'r') as f:
                    data = f.read()
                    probabilities[zealots] = float(data.split('\n')[0])
                
    probs = dict(sorted(probabilities.items()))

    return probs

def read_data_recovery(dir = os.getcwd()):
    probabilities = {}

    for dirpath, dirs, files in os.walk(dir):
        for file in files:
            if file.startswith("plasmares_"):
                zealots = float(file.split("_")[1].replace(".txt", ""))
                with open(os.path.join(dirpath, file), 'r') as f:
                    data = f.read()
                    probabilities[zealots] = float(data.split('\n')[0])
                
    probs = dict(sorted(probabilities.items()))

    return probs


def read_data_2dim(dir = os.getcwd()):
    probabilities = defaultdict(list)

    for dirpath, dirs, files in os.walk(dir):
        for file in files:
            if file.startswith("plasmares_"):
                zealots = int((file.split("_")[1]))
                contrarians = int((file.split("_")[2]).split(".")[0])
                with open(os.path.join(dirpath, file), 'r') as f:
                    data = f.read()
                    probabilities[(zealots, contrarians)] = float(data.split('\n')[0])
                
    probs = dict(sorted(probabilities.items()))

    return probs


"""
Compute range of values for which the probability should be computed
    input:
        N: group size
    output: 
        array of evenly spaced values within a given interval

"""
def compute_range(N = 100):

    # We compute the probabilities for maximum 70% of stubborn individuals
    max_stubborn = int(N * 0.8)

    # We want to have around 120 values of zealots for each group size, if possible
    desired_count = 120 
    step = 2
    # Calculate the step size to try to get close to `desired_count` values
    max_possible_values = (max_stubborn - 0) // step + 1
    while max_possible_values > desired_count:
        step += 2
        max_possible_values = (max_stubborn - 0) // step + 1
    
    # Generate the even numbers within the range
    return np.arange(2, max_stubborn + 1, step)