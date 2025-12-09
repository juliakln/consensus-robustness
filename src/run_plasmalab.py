"""
Generic and model-specific functions to write models in RML and properties in BLTL for Plasmalab

Monte-Carlo analysis of properties using Plasmalab command-line interface

The core function `run_consensus` handles directory creation, model generation,
property file generation, Plasmalab invocation, and result parsing.  
Model-specific wrapper functions (e.g., `stableconsensus_voter`) call it with 
the appropriate model-writing and property-writing functions, supporting both 
symmetric and asymmetric parameter configurations.

Results are saved in ../SMC_results/
- consensus
- consensus_switch
- recovery
- recoverygood
- consensusandrecovery
- consensusandrecoverygood

"""

import os
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from collections import defaultdict

from consensus_models import *
from consensus_properties import *

# change to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)


## Global parameters
# define file for model and property
model = "../models/consensus_model.rml"
property = "../models/consensus.bltl"

# get plasmalab cli path from environment variable
plasmacli = os.environ.get("PLASMALAB_PATH", "plasmacli.sh")



"""
Generic wrapper to tun Plasmalab with commandline and compute robustness analysis for consensus
    Parameters:
        stubborn: z (zealots), c (contrarians), b (both), none
        N: total population size
        majority: more than m% of population commits to same decision
        distance: difference of at least d between majority and those favouring opposite decision
        reaching: consensus is reached within t minutes from start of dynamics
        holding: group maintains consensus for at least h minutes
        stubborn_int: amount of stubborn individuals to compute robustness
        foldername: where to save the result files
        samples: how many samples to use for Monte Carlo
        write_model_fn: function that writes the model RML file
        write_property_fn: function that writes the property BLTL file
        folder_suffix: optional suffix such as '_x_' or '_y_'
        extra_args: optional dict of extra parameters passed to write_model_fn, usually ratex and ratey

    Returns: 
        dict: probabilities per number of stubborn individuals
"""
def run_consensus(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        stubborn_int,
        foldername,
        samples,
        write_model_fn,
        write_property_fn,
        folder_suffix="",
        extra_args=None,
):
    extra_args = extra_args or {}

    rx = str(extra_args['ratex']).replace('.', 'p') 
    ry = str(extra_args['ratey']).replace('.', 'p') 

    foldername += folder_suffix

    # --- directory creation ---
    os.makedirs(f'../SMC_results/{foldername}', exist_ok=True)
    # choose filename according to properties of model and property
    filename = f'../SMC_results/{foldername}/{N}N_{rx}x_{ry}y_{majority}m{distance}d{reaching}t{holding}h_{stubborn_int}{stubborn_type}.txt'

    if not os.path.exists(filename):

        # write model file (pass model-specific arguments)
        write_model_fn(stubborn_type, N, int(stubborn_int), **extra_args)

        # write properties
        write_property_fn(N, stubborn_type, majority, distance, reaching, holding)

        # run Plasmalab
        pcommand = (
            f"sh {plasmacli} "
            f"launch -m {model}:rml -r {property}:bltl "
            f"-a montecarlo -A\"Total samples\"={samples} -f proba --progress "
            f"-o {filename}"
        )
        subprocess.check_call(pcommand, shell=True)

    # --- collect results ---
    with open(filename, 'r') as f:
        prob = f.read()

    # return probability    
    return float(prob.split('\n')[0])



"""
Define specific wrappers for different models
"""
def stableconsensus_voter(stubborn_type, N, stubborn_int, ratex=1.0, ratey=1.0, majority=50, distance=None, reaching=35, holding=40, samples=4239, **kwargs):
    if distance is None:
        distance = int(N/10)
    return run_consensus(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        stubborn_int,
        'consensus/voter',
        samples,
        write_model_fn=write_voter_model,
        write_property_fn=write_property_stableconsensus,
        extra_args={"ratex": ratex, "ratey": ratey, **kwargs}
    )

def stableconsensus_voter_x(stubborn_type, N, stubborn_int, ratex=1.0, ratey=1.0, majority=50, distance=None, reaching=35, holding=40, samples=4239, **kwargs):
    if distance is None:
        distance = int(N/10)
    return run_consensus(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        stubborn_int,
        'consensus_x/voter',
        samples,
        write_model_fn=write_voter_model,
        write_property_fn=write_property_stableconsensus_asym_x,
        extra_args={"ratex": ratex, "ratey": ratey, **kwargs}
    )

def stableconsensus_voter_y(stubborn_type, N, stubborn_int, ratex=1.0, ratey=1.0, majority=50, distance=None, reaching=35, holding=40, samples=4239, **kwargs):
    if distance is None:
        distance = int(N/10)
    return run_consensus(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        stubborn_int,
        'consensus_y/voter',
        samples,
        write_model_fn=write_voter_model,
        write_property_fn=write_property_stableconsensus_asym_y,
        extra_args={"ratex": ratex, "ratey": ratey, **kwargs}
    )

def stableconsensus_crossinh(stubborn_type, N, stubborn_int, ratex=1.0, ratey=1.0, majority=50, distance=None, reaching=35, holding=40, samples=4239, **kwargs):
    if distance is None:
        distance = int(N/10)
    return run_consensus(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        stubborn_int,
        'consensus/crossinh',
        samples,
        write_model_fn=write_crossinh_model,
        write_property_fn=write_property_stableconsensus,
        extra_args={"ratex": ratex, "ratey": ratey, **kwargs}
    )

def stableconsensus_crossinh_x(stubborn_type, N, stubborn_int, ratex=1.0, ratey=1.0, majority=50, distance=None, reaching=35, holding=40, samples=4239, **kwargs):
    if distance is None:
        distance = int(N/10)
    return run_consensus(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        stubborn_int,
        'consensus_x/crossinh',
        samples,
        write_model_fn=write_crossinh_model,
        write_property_fn=write_property_stableconsensus_asym_x,
        extra_args={"ratex": ratex, "ratey": ratey, **kwargs}
    )

def stableconsensus_crossinh_y(stubborn_type, N, stubborn_int, ratex=1.0, ratey=1.0, majority=50, distance=None, reaching=35, holding=40, samples=4239, **kwargs):
    if distance is None:
        distance = int(N/10)
    return run_consensus(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        stubborn_int,
        'consensus_y/crossinh',
        samples,
        write_model_fn=write_crossinh_model,
        write_property_fn=write_property_stableconsensus_asym_y,
        extra_args={"ratex": ratex, "ratey": ratey, **kwargs}
    )

def stableconsensus_combined(stubborn_type, N, stubborn_int, ratex=1.0, ratey=1.0, majority=50, distance=None, reaching=35, holding=40, samples=4239, **kwargs):
    if distance is None:
        distance = int(N/10)
    return run_consensus(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        stubborn_int,
        'consensus/combined',
        samples,
        write_model_fn=write_combined_model,
        write_property_fn=write_property_stableconsensus,
        extra_args={"ratex": ratex, "ratey": ratey, **kwargs}
    )

def stableconsensus_combined_x(stubborn_type, N, stubborn_int, ratex=1.0, ratey=1.0, majority=50, distance=None, reaching=35, holding=40, samples=4239, **kwargs):
    if distance is None:
        distance = int(N/10)
    return run_consensus(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        stubborn_int,
        'consensus_x/combined',
        samples,
        write_model_fn=write_combined_model,
        write_property_fn=write_property_stableconsensus_asym_x,
        extra_args={"ratex": ratex, "ratey": ratey, **kwargs}
    )

def stableconsensus_combined_y(stubborn_type, N, stubborn_int, ratex=1.0, ratey=1.0, majority=50, distance=None, reaching=35, holding=40, samples=4239, **kwargs):
    if distance is None:
        distance = int(N/10)
    return run_consensus(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        stubborn_int,
        'consensus_y/combined',
        samples,
        write_model_fn=write_combined_model,
        write_property_fn=write_property_stableconsensus_asym_y,
        extra_args={"ratex": ratex, "ratey": ratey, **kwargs}
    )


## probability for switching consensus
def switchconsensus_voter(stubborn_type, N, stubborn_int, ratex=1.0, ratey=1.0, majority=50, distance=None, reaching=35, holding=40, samples=4239, **kwargs):
    if distance is None:
        distance = int(N/10)
    return run_consensus(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        stubborn_int,
        'consensus_switch/voter',
        samples,
        write_model_fn=write_voter_model,
        write_property_fn=write_property_switching,
        extra_args={"ratex": ratex, "ratey": ratey, **kwargs}
    )

def switchconsensus_crossinh(stubborn_type, N, stubborn_int, ratex=1.0, ratey=1.0, majority=50, distance=None, reaching=35, holding=40, samples=4239, **kwargs):
    if distance is None:
        distance = int(N/10)
    return run_consensus(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        stubborn_int,
        'consensus_switch/crossinh',
        samples,
        write_model_fn=write_crossinh_model,
        write_property_fn=write_property_switching,
        extra_args={"ratex": ratex, "ratey": ratey, **kwargs}
    )





"""
Generic wrapper to run Plasmalab with commandline and compute robustness analysis for recovery
    Parameters:
        stubborn: z (zealots), c (contrarians), b (both), none
        N: total population size
        majority: more than m% of population commits to same decision
        distance: difference of at least d between majority and those favouring opposite decision
        reaching: consensus is reached within t minutes from start of dynamics
        holding: group maintains consensus for at least h minutes
        stubborn_int: range of stubborn individuals to compute robustness
        filename: where to save the result files
        samples: how many samples to use for Monte Carlo
        write_model_fn: function that writes the model RML file
        write_property_fn: function that writes the property BLTL file
        folder_suffix: optional suffix such as '_x_' or '_y_'
        extra_args: optional dict of extra parameters passed to write_model_fn

    Returns: 
        dict: probabilities per number of stubborn individuals
"""
def run_recovery(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        recovery,
        stubborn_int,
        foldername,
        samples,
        write_model_fn,
        write_property_fn,
        folder_suffix="",
        extra_args=None,
):
    extra_args = extra_args or {}

    rx = str(extra_args['ratex']).replace('.', 'p') 
    ry = str(extra_args['ratey']).replace('.', 'p') 

    foldername += folder_suffix

    # --- directory creation ---
    os.makedirs(f'../SMC_results/{foldername}', exist_ok=True)
    # choose filename according to properties of model and property
    filename = f'../SMC_results/{foldername}/{N}N_{rx}x_{ry}y_{majority}m{distance}d{reaching}t{holding}h_{stubborn_int}{stubborn_type}.txt'

    if not os.path.exists(filename):

        # write model file (pass model-specific arguments)
        write_model_fn(stubborn_type, N, int(stubborn_int), **extra_args)

        # write properties
        write_property_fn(N, stubborn_type, majority, distance, reaching, holding, recovery)

        # run Plasmalab
        pcommand = (
            f"sh {plasmacli} "
            f"launch -m {model}:rml -r {property}:bltl "
            f"-a montecarlo -A\"Total samples\"={samples} -f proba --progress "
            f"-o {filename}"
        )
        subprocess.check_call(pcommand, shell=True)

    # --- collect results ---
    with open(filename, 'r') as f:
        prob = f.read()

    # return probability    
    return float(prob.split('\n')[0])


def recovery_voter(stubborn_type, N, majority, distance, reaching, holding, recovery, stubborn_int, samples, ratex=1.05, ratey=0.95):
    return run_recovery(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        recovery,
        stubborn_int,
        'recovery/voter',
        samples,
        write_model_fn=write_voter_model,
        write_property_fn=write_property_recovery,
        extra_args={"ratex": ratex, "ratey": ratey}
    )

def recovery_crossinh(stubborn_type, N, majority, distance, reaching, holding, recovery, stubborn_int, samples, ratex=1.05, ratey=0.95):
    return run_recovery(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        recovery,
        stubborn_int,
        'recovery/crossinh',
        samples,
        write_model_fn=write_crossinh_model,
        write_property_fn=write_property_recovery,
        extra_args={"ratex": ratex, "ratey": ratey}
    )

def recoverygood_voter(stubborn_type, N, majority, distance, reaching, holding, recovery, stubborn_int, samples, ratex=1.05, ratey=0.95):
    return run_recovery(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        recovery,
        stubborn_int,
        'recoverygood/voter',
        samples,
        write_model_fn=write_voter_model,
        write_property_fn=write_property_recoverygood,
        extra_args={"ratex": ratex, "ratey": ratey}
    )

def recoverygood_crossinh(stubborn_type, N, majority, distance, reaching, holding, recovery, stubborn_int, samples, ratex=1.05, ratey=0.95):
    return run_recovery(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        recovery,
        stubborn_int,
        'recoverygood/crossinh',
        samples,
        write_model_fn=write_crossinh_model,
        write_property_fn=write_property_recoverygood,
        extra_args={"ratex": ratex, "ratey": ratey}
    )

def recoverygood_combined(stubborn_type, N, majority, distance, reaching, holding, recovery, stubborn_int, samples, ratex=1.05, ratey=0.95):
    return run_recovery(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        recovery,
        stubborn_int,
        'recoverygood/combined',
        samples,
        write_model_fn=write_combined_model,
        write_property_fn=write_property_recoverygood,
        extra_args={"ratex": ratex, "ratey": ratey}
    )

def consensusandrecovery_voter(stubborn_type, N, majority, distance, reaching, holding, recovery, stubborn_int, samples, ratex=1.05, ratey=0.95):
    return run_recovery(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        recovery,
        stubborn_int,
        'consensusandrecovery/voter',
        samples,
        write_model_fn=write_voter_model,
        write_property_fn=write_property_recoveryandconsensus,
        extra_args={"ratex": ratex, "ratey": ratey}
    )

def consensusandrecovery_crossinh(stubborn_type, N, majority, distance, reaching, holding, recovery, stubborn_int, samples, ratex=1.05, ratey=0.95):
    return run_recovery(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        recovery,
        stubborn_int,
        'consensusandrecovery/crossinh',
        samples,
        write_model_fn=write_crossinh_model,
        write_property_fn=write_property_recoveryandconsensus,
        extra_args={"ratex": ratex, "ratey": ratey}
    )

def consensusandrecovery_combined(stubborn_type, N, majority, distance, reaching, holding, recovery, stubborn_int, samples, ratex=1.05, ratey=0.95):
    return run_recovery(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        recovery,
        stubborn_int,
        'consensusandrecovery/combined',
        samples,
        write_model_fn=write_combined_model,
        write_property_fn=write_property_recoveryandconsensus,
        extra_args={"ratex": ratex, "ratey": ratey}
    )


def consensusandrecoverygood_voter(stubborn_type, N, majority, distance, reaching, holding, recovery, stubborn_int, samples, ratex=1.05, ratey=0.95):
    return run_recovery(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        recovery,
        stubborn_int,
        'consensusandrecoverygood/voter',
        samples,
        write_model_fn=write_voter_model,
        write_property_fn=write_property_recoveryandconsensusgood,
        extra_args={"ratex": ratex, "ratey": ratey}
    )

def consensusandrecoverygood_crossinh(stubborn_type, N, majority, distance, reaching, holding, recovery, stubborn_int, samples, ratex=1.05, ratey=0.95):
    return run_recovery(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        recovery,
        stubborn_int,
        'consensusandrecoverygood/crossinh',
        samples,
        write_model_fn=write_crossinh_model,
        write_property_fn=write_property_recoveryandconsensusgood,
        extra_args={"ratex": ratex, "ratey": ratey}
    )

def consensusandrecoverygood_combined(stubborn_type, N, majority, distance, reaching, holding, recovery, stubborn_int, samples, ratex=1.05, ratey=0.95):
    return run_recovery(
        stubborn_type,
        N,
        majority,
        distance,
        reaching,
        holding,
        recovery,
        stubborn_int,
        'consensusandrecoverygood/combined',
        samples,
        write_model_fn=write_combined_model,
        write_property_fn=write_property_recoveryandconsensusgood,
        extra_args={"ratex": ratex, "ratey": ratey}
    )

