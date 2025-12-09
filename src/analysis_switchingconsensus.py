"""
Robustness analysis of switching consensus 

This script evaluates how robustly a group switches consensus
under variations of:
- consensus parameters (majority threshold, distance, reaching time, and holding time)
- group sizes

Analyses are performed across voter & cross-inhibition model, and 
systems containing disruptive individuals (zealots, contrarians, or both).

The script uses Plasmalab (Monte Carlo algorithm) to estimate satisfaction probabilities
for the BLTL property defining stable consensus.

Default parameter values: majority=50, distance=(N/10), reaching=35, holding=40 (= switching), samples=4239 (for Monte Carlo)

Resulting plots are saved in figures/

"""

import os
import matplotlib.pyplot as plt
import time

from consensus_models import *
from consensus_properties import *
from run_plasmalab import *

# change to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# define file for model and property
model = "../models/consensus_model.rml"
property = "../models/consensus.bltl"


# use dictionary to choose model function
model_functions = {
    'voter': switchconsensus_voter,
    'crossinh': switchconsensus_crossinh
}

# dictionary of disruptive individuals for plot labels
dict_stubborns = {
    "z": "zealots",
    "c": "contrarians"
}

# parameters for font sizes in plots
plt.rcParams.update({
    'legend.fontsize': 16,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'axes.titlesize': 18
})


"""
Compute range of values for which the probability should be computed
    input:
        N: group size
        fraction: maximum fraction of stubborn individuals
        step: step size between values
    output: 
        array of evenly spaced values within a given interval

"""
def compute_range(N = 100, fraction = 0.7):

    # We compute the probabilities for maximum 70% of stubborn individuals
    max_stubborn = int(N * fraction)

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

"""
Analyse satisfaction probability of switching consensus wrt amount of stubborns in the group
    model: 'voter' or 'crossinh'
    N: total group size
    stubborn_type: z (zealots), c (contrarians)
    ratex: rate for opinion x
    ratey: rate for opinion y
    plot: whether to plot the results

    Output: list of satisfaction probabilities for different amounts of stubborn individuals
"""
def analyse_switchconsensus_stubborns(model, N = 100, stubborn_type = 'z', ratex = 1.0, ratey = 1.0, plot = False, **kwargs):
    
    # compute time taken for analysis
    start_time = time.time()

    # select function according to model
    try:
        analyse_fn = model_functions[model]
    except KeyError:
        raise ValueError(f"Unknown model '{model}'")
    
    # Generate the even numbers within the range
    stubborn_range = compute_range(N = N, fraction = 0.7, step = 2)

    results = []

    for stubborn in stubborn_range:
        results.append(analyse_fn(stubborn_type, N, stubborn_int = stubborn, ratex = ratex, ratey = ratey, **kwargs))

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    if plot:
        fig = plt.figure(figsize=(6,6))
        stubborns = [100 * s / N for s in stubborn_range]
        plt.plot(stubborns, results, 'k', linewidth = 1.5)
        plt.ylim(0,1)
        plt.xlabel(f'Percentage of {dict_stubborns[stubborn_type]} in the group')
        plt.ylabel("Satisfaction probability")
        fig.savefig(f"../figures/consensus_switch_{model}_{N}N_{str(ratex).replace('.', 'p')}x_{str(ratey).replace('.', 'p')}y_{stubborn_type}.png")
        plt.close()
    
    return results

"""
Compute satisfaction probability of switching consensus for different consensus parameters distance and holding time
    model: 'voter' or 'crossinh'
    N: total group size
    stubborn_type: z (zealots), c (contrarians)
    ratex: rate for opinion x
    ratey: rate for opinion y

    Output: lineplot over #stubborn individuals for different settings
"""
def analyse_switchconsensus_stubborns_parameters(model, N = 100, stubborn_type = 'z', ratex = 1.0, ratey = 1.0, plot = False, **kwargs):

    # Generate the even numbers within the range
    stubborn_range = compute_range(N = N, fraction = 0.7, step = 2)

    # run analyses for different parameters
    baseline = analyse_switchconsensus_stubborns(model, N, stubborn_type, ratex, ratey)
    d_low1 = analyse_switchconsensus_stubborns(model, N, stubborn_type, ratex, ratey, distance = int(N*0.05))
    d_low2 = analyse_switchconsensus_stubborns(model, N, stubborn_type, ratex, ratey, distance = int(N*0.07))
    d_high1 = analyse_switchconsensus_stubborns(model, N, stubborn_type, ratex, ratey, distance = int(N*0.12))
    d_high2 = analyse_switchconsensus_stubborns(model, N, stubborn_type, ratex, ratey, distance = int(N*0.15))
    h_low1 = analyse_switchconsensus_stubborns(model, N, stubborn_type, ratex, ratey, holding = 10)
    h_low2 = analyse_switchconsensus_stubborns(model, N, stubborn_type, ratex, ratey, holding = 20)
    h_high1 = analyse_switchconsensus_stubborns(model, N, stubborn_type, ratex, ratey, holding = 40)
    h_high2 = analyse_switchconsensus_stubborns(model, N, stubborn_type, ratex, ratey, holding = 50)

       
    if plot:   
        fig = plt.figure(figsize=(6,6))
        stubborns = [100 * s / N for s in stubborn_range]
        plt.plot(stubborns, baseline, 'k', linewidth = 1.5, label = 'Baseline')
        plt.plot(stubborns, d_low1, linestyle = '-', color = 'b', label = 'd='+str(int(N*0.05)))
        plt.plot(stubborns, d_low2, linestyle = '--', color = 'b', label = 'd='+str(int(N*0.07)))
        plt.plot(stubborns, d_high1, linestyle = ':', color = 'b', label = 'd='+str(int(N*0.12)))
        plt.plot(stubborns, d_high2, linestyle = '-.', color = 'b', label = 'd='+str(int(N*0.15)))
        plt.plot(stubborns, h_low1, linestyle = '-', color = 'r', label = 's=10')
        plt.plot(stubborns, h_low2, linestyle = '--', color = 'r', label = 's=20')
        plt.plot(stubborns, h_high1, linestyle = ':', color = 'r', label = 's=40')
        plt.plot(stubborns, h_high2, linestyle = '-.', color = 'r', label = 's=50')
        plt.ylim(0,1)
        plt.xlabel(f'Percentage of {dict_stubborns[stubborn_type]} in the group')
        plt.ylabel("Satisfaction probability")
        plt.legend()
        fig.savefig(f"../figures/consensus_switch_{model}_parameters_{N}N_{str(ratex).replace('.', 'p')}x_{str(ratey).replace('.', 'p')}y_{stubborn_type}.png")
        plt.close()


    return baseline, d_low1, d_low2, d_high1, d_high2, h_low1, h_low2, h_high1, h_high2



def analyse_switchconsensus_stubborns_groupsizes(model, group_sizes, stubborn_type = 'z', ratex = 1.0, ratey = 1.0, plot = False, **kwargs):
    results = {}
    for N in group_sizes:
        results[N] = analyse_switchconsensus_stubborns(model, N, stubborn_type, ratex, ratey, **kwargs)

    if plot:
        fig = plt.figure(figsize=(6,6))
        colours = ['k', 'g', 'r', 'b', 'c', 'm', 'olive', 'orange', 'darkgreen', 'brown', 'darkviolet']

        for i, (N, probs) in enumerate(results.items()):
            stubborn_range = compute_range(N = N, fraction = 0.7)
            stubborns = [100 * s / N for s in stubborn_range]
            cs = colours[i % len(colours)]
            plt.plot(stubborns, probs, color = cs, label = f'N={str(N)}')

        plt.ylim(0,1)
        plt.xlabel(f'Percentage of {dict_stubborns[stubborn_type]} in the group')
        plt.ylabel("Satisfaction probability")
        fig.savefig(f"../figures/consensus_switch_{model}_groupsizes_{str(ratex).replace('.', 'p')}x_{str(ratey).replace('.', 'p')}y_{stubborn_type}.png")
        plt.close()

    return results




def main():
    # analyse probability to switch consensus wrt amount of stubborns
    analyse_switchconsensus_stubborns('voter', N = 100, stubborn_type = 'c', ratex = 1.0, ratey = 1.0, plot = True)
    # analyse probability to switch consensus wrt amount of stubborns for different consensus parameters
    analyse_switchconsensus_stubborns_parameters('voter', N = 100, stubborn_type = 'c', ratex = 1.0, ratey = 1.0, plot = True)
    # analyse probability to switch consensus wrt amount of stubborns for different group sizes
    group_sizes = [20,50,100,1000,4000]
    analyse_switchconsensus_stubborns_groupsizes('voter', group_sizes, stubborn_type = 'c', ratex = 1.0, ratey = 1.0, plot = True)


if __name__ == "__main__":
    main()

