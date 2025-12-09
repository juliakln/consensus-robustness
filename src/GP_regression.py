"""
Gaussian-Process Classification for Model-Checking Data
------------------------------------------------------

Based on: Smoothed Model Checking (Bortolussi et al., 2013), Rasmussen and Williams (2006)

This script estimates satisfaction probabilities produced by PlasmaLab
and fits a Gaussian Process (GP) model to them.

Workflow:
1. Generate training data: for each (p1, p2), run Monte-Carlo simulations and
   obtain empirical probabilities p̂.
2. Convert p̂ to latent values using a probit link (f = Φ⁻¹(p̂)), turning the
   classification problem into a GP regression problem.
3. Estimate observation noise from binomial variance and propagate it to
   the latent space.
4. Optimise GP hyperparameters by minimising the marginal likelihood.
5. Compute the GP posterior on a grid of test points and map predictions
   back to probabilities via Φ(f).

This provides a smooth probability surface and confidence bounds over the
parameter space.

Here: analysis of stable consensus wrt. reaching and holding times
"""


import sys
import warnings
import numpy as np
from scipy.special import erf 
from scipy.linalg import cholesky, solve  # scipy gives upper triangular matrix, numpy lower
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib import cm
import os
from scipy.optimize import minimize

from GP_kernels import *
from run_plasmalab import *

warnings.filterwarnings("ignore")


""" Global Parameters """

# correction for kernel computation to ensure numerical stability
correction = 1e-4

# change to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# define file location for model and property
model = "../models/consensus_model.rml"
property = "../models/consensus.bltl"

# use dictionary to choose model function
model_functions = {
    'voter': stableconsensus_voter,
    'crossinh': stableconsensus_crossinh
}


# parameters for font sizes in plots
plt.rcParams.update({
    'legend.fontsize': 16,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'axes.titlesize': 18
})


def to_latent(p, N=None, eps=1e-6):
    """p: probabilities or (counts) / N if N given. Returns latent f = Phi^{-1}(p)."""
    if N is not None:
        p = (p / float(N)).reshape(-1)
    p = np.clip(p, eps, 1.0 - eps)
    return norm.ppf(p).reshape(-1,1)


def gp_regression_posterior(X_train, y_train_latent, X_test, kernel, params, noise_var=1e-2, jitter=1e-6):
    """
    Classic GP regression posterior on latent targets.
    X_train: (n, d)
    y_train_latent: (n,1)
    X_test: (m, d)
    kernel: callable K(a,b,params)
    params: kernel params dict
    noise_var: observational noise variance (scalar) on latent scale
    Returns: f_mean (m,1), f_var (m,1) (predictive mean and diagonal predictive variance)
    """
    n = X_train.shape[0]
    K = kernel(X_train, X_train, params)  # (n,n)

    if np.isscalar(noise_var):
        K += (noise_var + jitter) * np.eye(n)
    else:
        K += np.diag(noise_var.reshape(-1)) + jitter * np.eye(n)

    # Cholesky
    L = cholesky(K, lower=True)  # L L^T = K
    # compute alpha = K^{-1} y via solves
    alpha = solve(L.T, solve(L, y_train_latent.reshape(-1,1)))
    Ks = kernel(X_test, X_train, params)  # (m,n)
    f_mean = Ks.dot(alpha)  # (m,1)
    # predictive variance diag:
    v = solve(L, Ks.T)      # (n,m)
    kss = kernel(X_test, X_test, params).diagonal().reshape(-1,1)  # (m,1)
    f_var = np.maximum(kss - np.sum(v * v, axis=0).reshape(-1,1), 1e-12)
    # The formula above: var(f*) = kss - Ks K^{-1} Ks^T
    # we computed kss diagonals and subtracted sum of cols of v^2
    return f_mean, f_var


def calc_mse(x, x_s, f, probs):
    """ Compute MSE for posterior function of the training data points
    
    Parameters:
        x: Training data input
        x_s: Test data input 
        f: Training data output
        probs: computed probabilities (posterior)
        
    Returns:  
        mse: Average distance of predictions to true training data output
    """
    predictions = []
    if len(x_s[0,:]) == 1:
        for i in x:
            idx = np.absolute(x_s - i).argmin()
            predictions.append(probs[idx])

        mse = np.square(np.subtract(f.reshape(1,-1), predictions)).mean()
    elif len(x_s[0,:]) == 2:
        probs = probs.reshape(len(x_s),-1)
        for i in x:
            idx = np.absolute(x_s - i).argmin()
            if idx >= len(x_s):
                idx = len(x_s) - 1
            predictions.append(probs[idx])
        mse = np.square(np.subtract(f.reshape(1,-1), predictions)).mean()

    return mse


def negative_log_marginal_likelihood(params_vec, x, y, kernel, noise_latent=0.001, jitter=1e-8):
    """
    Computes the negative log marginal likelihood (to minimize) for GP regression.
    Args:
        params_vec: [log_var, log_ell]
        x: training inputs (N, D)
        y: training targets (N, )
        kernel: kernel function, e.g. kernel_matern
        noise_latent: Gaussian noise in latent space
    """
    log_var, log_ell = params_vec
    var = np.exp(log_var)
    ell = np.exp(log_ell)


    # Compute kernel matrix with noise
    K = kernel(x, x, {'var': var, 'ell': ell}) + np.diag(noise_latent.reshape(-1)) + jitter * np.eye(len(x))
    
    try:
        L = np.linalg.cholesky(K)
    except np.linalg.LinAlgError:
        return np.inf  # Not positive definite → large penalty

    # flatten target
    yv = y.reshape(-1)

    # Solve for alpha = K^-1 y using Cholesky
    alpha = np.linalg.solve(L.T, np.linalg.solve(L, yv))

    # Log marginal likelihood
    log_marg_like = -0.5 * np.dot(yv, alpha)
    log_marg_like -= np.sum(np.log(np.diag(L)))
    log_marg_like -= 0.5 * len(x) * np.log(2 * np.pi)
    
    return -float(log_marg_like)  # negative for minimization



def optimize_hyperparams_gp(x, y, kernel, noise_latent=0.02, var_init=1.0, ell_init=10.0):
    """
    ze GP kernel hyperparameters (var, ell) by maximizing log marginal likelihood.
    """
    print("🔍 Optimizing hyperparameters for GP regression ...")

    res = minimize(
        negative_log_marginal_likelihood,
        x0=np.log([var_init, ell_init]),
        args=(x, y, kernel, noise_latent),
        bounds=[(np.log(1e-6), np.log(10.0)), (np.log(1e-2), np.log(200.0)), (np.log(1e-6), np.log(0.5))],
        method='L-BFGS-B',
        options={'disp': True, 'maxiter': 200}
    )

    var_opt, ell_opt = np.exp(res.x)
    print(f"✅ Optimized hyperparameters: var = {var_opt:.6f}, ell = {ell_opt:.6f}")
    return {'var': var_opt, 'ell': ell_opt, 'noise_latent': noise_latent}



def analyse(modelname, kernel, var, var_bounds, ell, ell_bounds, scale, ptest1, ptest2, trainingpoints, testpoints, zealots, paramValueSet=None, paramValueOutput=None):
    """ 
    For 2 dimensions -> vary 2 parameters p1, p2
    Analyze satisfaction probability to find out if function is robust
    Then compute GPC of satisfaction probabilty
    Args:
        modelname: voter or crossinh
        kernel: kernel function
        var: initial variance for kernel
        var_bounds: bounds for variance optimization
        ell: initial lengthscale for kernel
        ell_bounds: bounds for lengthscale optimization
        scale: number of observations per input point
        ptest1: range of parameter 1 for test set
        ptest2: range of parameter 2 for test set
        trainingpoints: number of training points
        testpoints: number of test points for posterior - should be a square number
        zealots: number of zealots in model
        paramValueSet: training inputs / parameter values (N,2)
        paramValueOutput: training outputs / probabilities (N,1)
        plot: whether to plot the results
    
    Returns:
        p: predictive probabilities

    """
    if paramValueSet is None or paramValueOutput is None:
        paramValueSet, paramValueOutput = create_trainingdata_stableconsensus_reaching_holding(model=modelname, N=100, stubborn_type='z', stubborn_int=zealots, trainingpoints=trainingpoints)


    # Compute satisfaction numbers from probabilities
    counts = (paramValueOutput.reshape(-1) * scale).astype(int)
    counts = np.clip(counts, 0, scale)  # ensure counts are within valid range
    Nsim = scale
    X_gp = paramValueSet #[mask]
    p_hat = (counts / float(Nsim)).reshape(-1,1)

    # transform to latent:
    y_latent = to_latent(p_hat, N=None)  # already p_hat
        
    # observational noise in latent space
    phi = (1/np.sqrt(2*np.pi)) * np.exp(-0.5 * (y_latent.ravel()**2))
    var_p = (p_hat.ravel() * (1-p_hat.ravel())) / float(Nsim)
    est_var_latent = (1.0 / (phi + 1e-12))**2 * var_p
    #noise_latent = np.clip(np.median(est_var_latent), 1e-6, 0.01) --> for single value
    noise_latent = np.clip(est_var_latent, 1e-12, 0.01)

    # Optimize hyperparameters of kernel
    res = minimize(
        negative_log_marginal_likelihood, 
        x0=np.log([var, ell]), 
        args=(X_gp, y_latent, kernel, noise_latent), 
        bounds=[var_bounds, ell_bounds],   
        method='L-BFGS-B', 
        options={'disp': True, 'maxiter': 200})
   
    params = {'var': np.exp(res.x[0]), 
            'ell': np.exp(res.x[1]),        
            'ell_dim': [1,1], 
            'var_b': 1,
            'off': 1}
    


    # define test set for which posterior is derived
    npoints = int(np.ceil(testpoints**(1/2)))  # number of points per dimension
    testset = np.zeros((testpoints,2))  # create grid of equally spaced input points
    lin1 = np.linspace(ptest1[0], ptest1[1], npoints)
    lin2 = np.linspace(ptest2[0], ptest2[1], npoints)
    testset[:,0] = np.repeat(lin1, npoints)
    testset[:,1] = np.tile(lin2, npoints)

    # GP regression on latent targets:
    f_mean_latent, f_var_latent = gp_regression_posterior(X_gp, y_latent, testset, kernel, params, noise_var=noise_latent)

    # map back to probabilities:
    prob_mean = norm.cdf(f_mean_latent)  # (Ns,1)

    # reshape to grid as in your get_posterior
    probabilities = prob_mean.reshape(npoints, npoints).T

    phi_f = (1/np.sqrt(2*np.pi)) * np.exp(-0.5 * f_mean_latent**2)
    var_p = (phi_f**2) * f_var_latent
    lower_p = (np.clip(prob_mean - 1.96 * np.sqrt(var_p), 0, 1)).reshape(npoints, npoints).T
    upper_p = (np.clip(prob_mean + 1.96 * np.sqrt(var_p), 0, 1)).reshape(npoints, npoints).T

    # compute MSE
    mse = calc_mse(paramValueSet, testset, paramValueOutput, probabilities)
    print("MSE: ", mse)


    # plot posterior
    p1 = np.unique(testset[:,0])
    p2 = np.unique(testset[:,1])

    # # Surface Plot with mean posterior and training points
    # l1, l2 = np.meshgrid(p1, p2)
    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(10,10))
    # surf = ax.plot_surface(l1, l2, probabilities, cmap=cm.coolwarm, vmin=0, vmax=1)
    # ax.scatter(X_gp[:,0], X_gp[:,1], p_hat, c='black')
    # plt.xlabel(plabel1 + '$_*$')
    # plt.ylabel(plabel2 + '$_*$', rotation=0)
    # plt.title('Posterior distribution of consensus in ' + modelname + ' model with ' + str(zealots) + ' stubborns\nMean posterior surface and training data points')
    # cbar = fig.colorbar(surf, shrink=0.5)
    # cbar.set_label('Probability')
    # plt.savefig(f'../../figures/GP10002_{modelname}_{str(zealots)}_posteriorSurf.png', dpi=300)

    # # Surface Plot with confidence bounds and training points
    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(10,10))
    # surf = ax.plot_surface(l1, l2, lower_p, cmap=cm.coolwarm, vmin=0, vmax=1)
    # surf = ax.plot_surface(l1, l2, upper_p, cmap=cm.coolwarm, vmin=0, vmax=1)
    # ax.scatter(l1, l2, probabilities, c='black')
    # plt.xlabel(plabel1 + '$_*$')
    # plt.ylabel(plabel2 + '$_*$', rotation=0)
    # plt.title('Posterior distribution of stable consensus in ' + modelname + ' model with ' + str(zealots) + ' stubborns\nPosterior surface with 95% confidence bounds and test data points')
    # cbar = fig.colorbar(surf, shrink=0.5)
    # cbar.set_label('Probability')
    # plt.savefig(f'../../figures/GP10002_{modelname}_{str(zealots)}_posteriorSurfConf.png', dpi=300)


    return X_gp, p_hat, p1, p2, probabilities, lower_p, upper_p
    



    

def create_trainingdata_stableconsensus_reaching_holding(model, N, stubborn_type, stubborn_int, ratex=1.0, ratey=1.0, trainingpoints=100, 
                                         p1_range=[0,20], p2_range=[0,20], trajs=100, plot=False):
    """
    Create training data for 2-dimensional input

    Args:
        model: voter or crossinh
        N: group size
        stubborn_type: 'z' for zealots, 'c' for contrarians
        stubborn_int: number of stubborn individuals
        ratex: rate for opinion x
        ratey: rate for opinion y
        trainingpoints: number of training points
        p1_range: range of parameter 1
        p2_range: range of parameter 2
        trajs: number of trajectories per data point
        plot: whether to plot training data
    Returns:          
        paramValueSet: training inputs / parameter values (N,2)
        paramValueOutput: training outputs / probabilities (N,1)               
    """
    points = int(np.ceil(trainingpoints**(1/2)))  # number of points per dimension
    X = np.zeros((trainingpoints,2))  # create grid of equally spaced input points
    lin1 = np.linspace(p1_range[0], p1_range[1], points)
    lin2 = np.linspace(p2_range[0], p2_range[1], points)
    X[:,0] = np.repeat(lin1, points)
    X[:,1] = np.tile(lin2, points)

    analyse_fn = model_functions[model]

    # save number of satisfactions for each parameter combination
    satisfactions = defaultdict(list)

    for params in X:
        reaching = params[0]
        holding = params[1]

        result = analyse_fn(stubborn_type, N, stubborn_int, ratex, ratey, reaching=int(reaching), holding=int(holding), samples=trajs)
        satisfactions[(np.round(reaching,3), np.round(holding,3))] = result


    # Training data
    f = []
    for key in sorted(satisfactions):
        f.append(satisfactions[key])

    f = np.array(f).reshape(-1,1)

    if plot:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(12,6))
        scat = ax.scatter(X[:,0], X[:,1], f, c=f, cmap=cm.coolwarm, vmin=0, vmax=1)
        ax.set_xlabel('reaching time')
        ax.set_ylabel('holding time', rotation=0)
        cbar = fig.colorbar(scat, shrink=0.5)
        cbar.ax.set_ylabel('probability', rotation=0)
        ax.set_title(f'Satisfaction probability of stable consensus in {model} model with {str(stubborn_int)} {stubborn_type} stubborns (Training set)')
        plt.tight_layout()
        plt.savefig(f'../figures/GP_consensus_reachingholding_{model}_{N}N_{str(ratex).replace('.', 'p')}x_{str(ratey).replace('.', 'p')}y_{stubborn_type}_training.png', dpi=300)
        plt.close()

    return X, f


def main():
    print('-----VOTER MODEL-----\n')

    # --- Configuration parameters ---
    model = "voter"
    stubborn_type = "z"
    zealots = 0
    ratex = 1.0
    ratey = 1.0

    N = 100     # population size for model
    trainingpoints = 324
    testpoints = 400
    observations = 1000
    p1_range = [0, 250] # reaching time
    p2_range = [0, 250] # holding time

    # --- Get training data ---
    X_train, y_train = create_trainingdata_stableconsensus_reaching_holding(
        model,
        N=N,
        stubborn_type=stubborn_type,
        stubborn_int=zealots,
        ratex=ratex,
        ratey=ratey,
        trainingpoints=trainingpoints,
        p1_range=p1_range,
        p2_range=p2_range,
        trajs=observations,
        plot=False,
    )

    # --- Run GP analysis ---
    X_gp, p_hat, p1, p2, probabilities, lower_p, upper_p = analyse(
        modelname=model,
        kernel=kernel_matern32,
        var=1.0,
        var_bounds=(np.log(1e-6), np.log(3.0)),
        ell=10.0,
        ell_bounds=(np.log(1), np.log(400.0)),
        scale=observations,
        ptest1=p1_range,
        plabel1="reaching time",
        ptest2=p2_range,
        plabel2="holding time",
        trainingpoints=trainingpoints,
        Ns=testpoints,
        zealots=zealots,
        paramValueSet=X_train,
        paramValueOutput=y_train,
    )

    # --- Plot posterior probability ---
    fig = plt.figure(figsize=(8, 8))
    plt.contourf(p1, p2, probabilities, levels=np.linspace(0, 1, 11))
    plt.xlabel("reaching$_*$")
    plt.ylabel("holding$_*$")
    plt.title(
        f"Posterior consensus probability — {model} model, {zealots} stubborns"
    )
    cbar = plt.colorbar()
    cbar.set_label("Satisfaction Probability")

    # --- Build clean filename suffix ---
    suffix = (
        f"{model}_{N}N_"
        f"{str(ratex).replace('.', 'p')}x_"
        f"{str(ratey).replace('.', 'p')}y_"
        f"{stubborn_type}{zealots}"
    )

    fig.savefig(f"../figures/GP_consensus_posterior_{suffix}.png", dpi=300)
    plt.close(fig)

    print("Done.")


if __name__ == "__main__":
    sys.exit(main())