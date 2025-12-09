"""
Robustness analysis of reaching a stable consensus 

This script evaluates how robustly a group reaches and maintains a stable consensus
under variations of:
- consensus parameters (majority threshold, distance, reaching time, and holding time)
- group sizes
- consensus opinions separately (X vs Y)
- rate combinations favoring opinion X vs opinion Y

Analyses are performed across voter & cross-inhibition model, and 
systems containing disruptive individuals (zealots, contrarians, or both).

The script uses Plasmalab (Monte Carlo algorithm) to estimate satisfaction probabilities
for the BLTL property defining stable consensus.

Default parameter values: majority=50, distance=(N/10), reaching=35, holding=40, samples=4239 (for Monte Carlo)

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


## Global parameters

# define file for model and property
model = "../models/consensus_model.rml"
property = "../models/consensus.bltl"


# use dictionary to choose model function
model_functions = {
    'voter': stableconsensus_voter,
    'crossinh': stableconsensus_crossinh
}
model_functions_x = {
    'voter': stableconsensus_voter_x,
    'crossinh': stableconsensus_crossinh_x
}
model_functions_y = {
    'voter': stableconsensus_voter_y,
    'crossinh': stableconsensus_crossinh_y
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
Analyse stable consensus wrt. the amount of stubborns in the group
    model: 'voter' or 'crossinh'
    N: total group size
    stubborn_type: z (zealots), c (contrarians)
    ratex: rate for opinion x
    ratey: rate for opinion y
    plot: whether to plot the results
    **kwargs: additional parameters for analysis function

    Return: probabilities for results
"""
def analyse_stableconsensus_stubborns(model, N = 100, stubborn_type = 'z', ratex = 1.0, ratey = 1.0, plot = False, **kwargs):
    
    # compute time for analysis
    start_time = time.time()

    # select function according to model
    try:
        analyse_fn = model_functions[model]
    except KeyError:
        raise ValueError(f"Unknown model '{model}'")

    # range of stubborn individuals for which the probability is evaluated
    stubborn_range = compute_range(N = N, fraction = 0.7)

    results = []

    for stubborn in stubborn_range:
        
        results.append(analyse_fn(stubborn_type, N, stubborn_int = stubborn, ratex = ratex, ratey = ratey, **kwargs))
        

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    if plot:
        fig = plt.figure(figsize=(6,6))
        stubborns = [100 * s / N for s in stubborn_range]
        plt.plot(stubborns, results, 'k', linewidth = 1.5, label = 'Stable Consensus')
        plt.ylim(0,1)
        plt.xlabel(f'Percentage of {dict_stubborns[stubborn_type]} in the group')
        plt.ylabel("Satisfaction probability")
        plt.legend()
        fig.savefig(f"../figures/consensus_{model}_{N}N_{str(ratex).replace('.', 'p')}x_{str(ratey).replace('.', 'p')}y_{stubborn_type}.png")
        plt.close()

    return results


"""
Analyse stable consensus wrt. the amount of stubborns in the group
    for different consensus parameters
    model: 'voter' or 'crossinh'
    N: total group size
    stubborn_type: z (zealots), c (contrarians)
    ratex: rate for opinion x
    ratey: rate for opinion y

    Return: probabilities for baseline, m_low, m_high, d_low, d_high, t_low, t_high, h_low, h_high
"""
def analyse_stableconsensus_stubborns_parameters(model, N = 100, stubborn_type = 'z', ratex = 1.0, ratey = 1.0, plot = False, **kwargs):   

    # range of stubborn individuals for which the probability is evaluated
    stubborn_range = compute_range(N = N, fraction = 0.7)

    # run analyses with modified parameters
    baseline = analyse_stableconsensus_stubborns(model, N = N, stubborn_type = stubborn_type, ratex = ratex, ratey = ratey)
    m_low = analyse_stableconsensus_stubborns(model, N = N, stubborn_type = stubborn_type, ratex = ratex, ratey = ratey, majority = 35)
    m_high = analyse_stableconsensus_stubborns(model, N = N, stubborn_type = stubborn_type, ratex = ratex, ratey = ratey, majority = 65)
    d_low = analyse_stableconsensus_stubborns(model, N = N, stubborn_type = stubborn_type, ratex = ratex, ratey = ratey, distance = 1)
    d_high = analyse_stableconsensus_stubborns(model, N = N, stubborn_type = stubborn_type, ratex = ratex, ratey = ratey, distance = 20)
    t_low = analyse_stableconsensus_stubborns(model, N = N, stubborn_type = stubborn_type, ratex = ratex, ratey = ratey, reaching = 20)
    t_high = analyse_stableconsensus_stubborns(model, N = N, stubborn_type = stubborn_type, ratex = ratex, ratey = ratey, reaching = 50)
    h_low = analyse_stableconsensus_stubborns(model, N = N, stubborn_type = stubborn_type, ratex = ratex, ratey = ratey, holding = 25)
    h_high = analyse_stableconsensus_stubborns(model, N = N, stubborn_type = stubborn_type, ratex = ratex, ratey = ratey, holding = 55)


    if plot:
        # plot results
        fig = plt.figure(figsize=(6,6))
        stubborns = [100 * s / N for s in stubborn_range]
        plt.plot(stubborns, baseline, 'k', linewidth = 1.5, label = 'Baseline')
        plt.plot(stubborns, m_low, linestyle = '-', color = '#ffaa00', label = 'm=35')
        plt.plot(stubborns, m_high, linestyle = '--', color = '#ff0000', label = 'm=65')
        plt.plot(stubborns, d_low, linestyle = '-', color = 'b', label = 'd=1')
        plt.plot(stubborns, d_high, linestyle = '--', color = 'b', label = 'd=20')
        plt.plot(stubborns, t_low, linestyle = '-', color = 'g', label = 't=20')
        plt.plot(stubborns, t_high, linestyle = '--', color = 'g', label = 't=50')
        plt.plot(stubborns, h_low, linestyle = '-', color = 'r', label = 'h=25')
        plt.plot(stubborns, h_high, linestyle = '--', color = 'r', label = 'h=55')
        plt.ylim(0,1)
        plt.xlabel(f'Percentage of {dict_stubborns[stubborn_type]} in the group')
        plt.ylabel("Satisfaction probability")
        plt.legend()
        fig.savefig(f"../figures/consensus_{model}_parameters_{N}N_{str(ratex).replace('.', 'p')}x_{str(ratey).replace('.', 'p')}y_{stubborn_type}.png")
        plt.close()


    return baseline, m_low, m_high, d_low, d_high, t_low, t_high, h_low, h_high



"""
Analyse stable consensus wrt. the amount of stubborns in the group
    for different group sizes
    model: 'voter' or 'crossinh'
    group_sizes: list of group sizes
    stubborn_type: z (zealots), c (contrarians)
    ratex: rate for opinion x
    ratey: rate for opinion y
    plot: whether to plot the results

    Return: dictionary of probabilities for all group sizes
"""
def analyse_stableconsensus_stubborns_groupsizes(model, group_sizes, stubborn_type, ratex, ratey, plot = False, **kwargs):
    results = {}
    for N in group_sizes:
        results[N] = analyse_stableconsensus_stubborns(model, N = N, stubborn_type = stubborn_type, ratex = ratex, ratey = ratey, **kwargs)

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
        fig.savefig(f"../figures/consensus_{model}_groupsizes_{str(ratex).replace('.', 'p')}x_{str(ratey).replace('.', 'p')}y_{stubborn_type}.png")
        plt.close()

    return results



"""
Analyse stable consensus of X and Y separately wrt. the amount of stubborns in the group 
    model: 'voter' or 'crossinh'
    N: total group size
    stubborn_type: z (zealots), c (contrarians)
    ratex: rate for opinion x
    ratey: rate for opinion y

    Return: probabilities for results_x, results_y
"""
def analyse_stableconsensus_XY_stubborns(model, N = 100, stubborn_type = 'z', ratex = 1.0, ratey = 1.0, plot = False, **kwargs):

    # compute time for analysis
    start_time = time.time()

    # select function according to model
    try:
        analyse_fn_x = model_functions_x[model]
        analyse_fn_y = model_functions_y[model]
    except KeyError:
        raise ValueError(f"Unknown model '{model}'")

    # range of stubborn individuals for which the probability is evaluated
    stubborn_range = compute_range(N = N, fraction = 0.7)

    results_x, results_y = [], []

    # compute satisfaction probability for all amounts of stubborns
    for stubborn in stubborn_range:

        results_x.append(analyse_fn_x(stubborn_type, N, stubborn_int = stubborn, ratex = ratex, ratey = ratey, **kwargs))
        results_y.append(analyse_fn_y(stubborn_type, N, stubborn_int = stubborn, ratex = ratex, ratey = ratey, **kwargs))

    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    if plot:
        # plot results
        fig = plt.figure(figsize=(6,6))
        stubborns = [100 * s / N for s in stubborn_range]
        plt.plot(stubborns, results_x, 'r', linewidth = 1.5, label = 'X')
        plt.plot(stubborns, results_y, 'b', linewidth = 1.5, label = 'Y')
        plt.ylim(0,1)
        plt.xlabel(f'Percentage of {dict_stubborns[stubborn_type]} in the group')
        plt.ylabel("Satisfaction probability")
        plt.legend()
        fig.savefig(f"../figures/consensus_xy_{model}_{N}N_{str(ratex).replace('.', 'p')}x_{str(ratey).replace('.', 'p')}y_{stubborn_type}.png")
        plt.close()

    return results_x, results_y 



"""
Analyse stable consensus of X and Y separately wrt. the amount of stubborns in the group
    for different rate combinations
    model: 'voter' or 'crossinh'
    N: total group size
    stubborn_type: z (zealots), c (contrarians)
    rate_pairs: list of (ratex, ratey) tuples

    Output: lineplot over #stubborn individuals for X and Y for different rate combinations
    Return: dictionary of probabilities for all rate combinations

NOTE: maybe change colours/linestyles representations
"""
def analyse_stableconsensus_XY_stubborns_rates(model, N, stubborn_type, rate_pairs, plot = False, **kwargs):
 
    stubborn_range = compute_range(N = N, fraction = 0.7)

    results = {}
    # compute satisfaction probability for all rate combinations over stubborn range
    for ratex, ratey in rate_pairs:

        results_x, results_y = analyse_stableconsensus_XY_stubborns(model, N = N, stubborn_type = stubborn_type, ratex = ratex, ratey = ratey, **kwargs)

        results[(ratex, ratey)] = (results_x, results_y)

    if plot:
        # plot results
        fig = plt.figure(figsize=(6,6))
        stubborns = [100 * s / N for s in stubborn_range]
        linestyles = ['-', '--', ':', '-.', (0, (3, 1, 1, 1)), (0, (5, 1))]
        rate_legend_handles = []

        for i, ((ratex, ratey), (probs_x, probs_y)) in enumerate(results.items()):
            ls = linestyles[i % len(linestyles)]

            plt.plot(stubborns, probs_x, 'r', linestyle = ls)
            plt.plot(stubborns, probs_y, 'b', linestyle = ls)
            
            # Dummy line for rates legend
            h = plt.Line2D([0], [0], color='black', linestyle=ls, label=f'qx={str(ratex)}, qy={str(ratey)}')
            rate_legend_handles.append(h)

        legend_A = plt.legend(handles=[plt.Line2D([0], [0], color='r', label='X'),
                                        plt.Line2D([0], [0], color='b', label='Y')],
                                loc='upper right')
        plt.gca().add_artist(legend_A)
        plt.legend(handles=rate_legend_handles, loc = 'upper right', bbox_to_anchor=(1, 0.82))
        plt.ylim(0,1)
        plt.xlabel(f'Percentage of {dict_stubborns[stubborn_type]} in the group')
        plt.ylabel("Satisfaction probability")
        fig.savefig(f"../figures/consensus_xy_{model}_rates_{N}N_{stubborn_type}.png")
        plt.close()

    return results


"""
Analyse stable consensus of X and Y separately wrt. the amount of stubborns in the group
    for different group sizes
    model: 'voter' or 'crossinh'
    group_sizes: list of group sizes
    stubborn_type: z (zealots), c (contrarians)
    ratex: rate for opinion x
    ratey: rate for opinion y

    Output: lineplot over #stubborn individuals for X and Y for different group sizes
    Return: dictionary of probabilities for all group sizes
"""
def analyse_stableconsensus_XY_stubborns_groupsizes(model, group_sizes, stubborn_type, ratex, ratey, plot = False, **kwargs):
    stubborn_range = compute_range(N = N, fraction = 0.7)

    results = {}
    for N in group_sizes:
        results[N] = analyse_stableconsensus_XY_stubborns(model, N = N, stubborn_type = stubborn_type, ratex = ratex, ratey = ratey, **kwargs)

    if plot:
        # plot results
        fig = plt.figure(figsize=(6,6))
        stubborns = [100 * s / N for s in stubborn_range]
        linestyles = ['-', '--', ':', '-.', (0, (3, 1, 1, 1)), (0, (5, 1))]
        size_legend_handles = []

        for i, (N, (probs_x, probs_y)) in enumerate(results.items()):
            ls = linestyles[i % len(linestyles)]

            plt.plot(stubborns, probs_x, 'r', linestyle = ls)
            plt.plot(stubborns, probs_y, 'b', linestyle = ls)
            
            # Dummy line for rates legend
            h = plt.Line2D([0], [0], color='black', linestyle=ls, label=f'N={str(N)}')
            size_legend_handles.append(h)

        legend_A = plt.legend(handles=[plt.Line2D([0], [0], color='r', label='X'),
                                        plt.Line2D([0], [0], color='b', label='Y')],
                                loc='upper right')
        plt.gca().add_artist(legend_A)
        plt.legend(handles=size_legend_handles, loc = 'upper right', bbox_to_anchor=(1, 0.82))
        plt.ylim(0,1)
        plt.xlabel(f'Percentage of {dict_stubborns[stubborn_type]} in the group')
        plt.ylabel("Satisfaction probability")
        fig.savefig(f"../figures/consensus_xy_{model}_groupsizes_{str(ratex).replace('.', 'p')}x_{str(ratey).replace('.', 'p')}y_{stubborn_type}.png")
        plt.close()

    return results




def main():
    # analyse probability to reach a stable consensus wrt amount of stubborns for different consensus parameters
    analyse_stableconsensus_stubborns_parameters(model = 'voter', N = 100, stubborn_type = 'z', ratex = 1.0, ratey = 1.0)
    # analyse probability to reach a stable consensus wrt amount of stubborns for different group sizes
    analyse_stableconsensus_stubborns_groupsizes(model = 'voter', group_sizes = [20, 50, 100, 1000, 4000], stubborn_type = 'z', ratex = 1.0, ratey = 1.0)
    # analyse probability to reach a stable consensus for X and Y separately wrt amount of stubborns
    analyse_stableconsensus_XY_stubborns(model = 'voter', N = 100, stubborn_type = 'z', ratex = 1.05, ratey = 0.95)
    # analyse probability to reach a stable consensus for X and Y separately wrt amount of stubborns for different rate combinations
    rate_pairs = [(1.0, 1.0), (1.01, 0.99), (1.05, 0.95), (1.1, 0.9)]
    analyse_stableconsensus_XY_stubborns_rates(model = 'voter', N = 100, stubborn_type = 'z', rate_pairs = rate_pairs)
    # analyse probability to reach a stable consensus for X and Y separately wrt amount of stubborns for different group sizes
    group_sizes = [20, 50, 100, 1000, 4000]
    analyse_stableconsensus_XY_stubborns_groupsizes(model = 'voter', group_sizes = group_sizes, stubborn_type = 'z', ratex = 1.05, ratey = 0.95)


if __name__ == "__main__":
    main()




