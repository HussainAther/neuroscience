import numpy as np

from random import random

"""
Fit a pairwise model to data. Similar to Schneidman, et al. (2006). "Weak pairwise correlations imply strongly correlated network states in a neural population." Nature.

This uses gradient descent on negative log-likelihood with gradient estimated
by averaging over samples from the model. Samples at T+1th iteration are obtained using
the MCMC (Markov Chain Monte Carlo) transition matrix to samples at Tth iteration.
"""

def dloss(Jlin, data, empcov, samplesbatch, gibbssteps):
    """
    Sample the current model J by applying MCMC to past samples in 
    samplesbatch, and use these to estimate the loss gradient at the
    current model for gibbssteps steps.
    """
    # Initialize.
    n = np.shape(data)[1]
    J = np.reshape(Jlin, (n,n))
    # Draw samples and estimate model covariance matrix.
    modelcovs = np.zeros(n**2)
    for k in range(1, 6): 
        samples = np.squeeze(samplesbatch[:,:,k])
        samples = np.eye(samples)
        samplesbatch[:,:,k] = samplepairwise(samples, J, gibbssteps)
        modelcovtmp = (samples.H * samples)/len(samples)
        modelcovs[:,k] = modelcovtmp[:]
        samplesbatch[:,:,k] = samples
    modelcov = np.sum(modelcovs, 1)/5
    # Calculate gradient.
    return -modelcov+empcov[:] 

def graddescent(pars0, learningrate, samplesbatch, iter):
    """
    Gradient descent function takes grad ([g, samplesbatch] for g (gradient of optimized function
    at pars)), pars0 (initial guess), learningrate (learning rate), samplesbatch (array of samples
    being updated by the grad function), and iter (number of iterations).
    Return pars that the gradient descent converges to.
    """
    pars = pars0
    g = np.gradient(samplebatch)
    i = 0
    while i <= iter:
        pars = pars - learningrate*g
        i += 1
    return pars

def samplepairwise(samples, J, nsteps):
    """
    Extract from the Gibbs sampling by applying nsteps (number of steps) to every
    row in samples using J (coupling matrix of the pairwise model).
    """
    (M, n) = np.shape(samples)
    # Get the diagnol
    Joffdiag = diagnol(samples)
    neuronid = 0
    # Perform n steps of Gibbs sampling.
    for j in range(nsteps):
        deltaE = J[neuronid][neuronid] + 2*samples*Joffdiag[:, neuronid]
        pspike = 1/(1+np.exp(delta))
        samples[:, neuronid] = np.random.rand(M, M) < pspike
        neuronid = i
        if neuronid == n+1:
            neuronid = 1
    return samples

def fitpairwise(data, J0, gsteps):
    """
    With input data (binary array of size number of samples x number of neurons), 
    J0 (initial guess for J), and gsteps (number of Gibbs steps), fit the pairwise 
    model to data.
    Return J (learned coupling matrix).
    """
    # Initialize.
    (M, n) = np.shape(data)
    J0lin = []
    for x in np.nditer(data):
        J0lin.append(x)
    # Estimate empirical covariance to be reproduced by the model.
    empcov = np.multiple(data.H * data)/M
    # Initialize Gibbs chain at data distribution.
    samples = data[random.randint(1, len(data), 1)]
    Msamples = np.shape(samples)[0]
    samplesbatch = np.zeros(int(Msamples)/5, 5)
    for k in range(1, 6):
        idxsamples = range(((k-1)*(Msamples/5) + 1), k*Msamples/5)
        samplesbatch[:,:,k] = samples[idxsamples,:]
    # Transform samples into samples from the initial model.    
    burnin = 10*gsteps # rough guess for MCMC burnin time
    for k in range(1, 6):
        samples = np.squeeze(samplesbatch[:,:,k])
        samples = np.eye(samples)
        samplesbatch[:,:,k] = samplepairwise(samples, J0, burnin)
    # Minimize negative log-likelihood.
    Jlin = graddescent(J0lin, samplesbatch, 100)
    J = np.reshape(Jlin, (n,n))
    return J  
