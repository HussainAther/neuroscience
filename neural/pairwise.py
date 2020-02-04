import numpy as np

from random import random

"""
Fit a pairwise model to data. Similar to Schneidman, et al. (2006). "Weak pairwise correlations imply strongly correlated network states in a neural population." Nature.

This uses gradient descent on negative log-likelihood with gradient estimated
by averaging over samples from the model. Samples at T+1th iteration are obtained using
the MCMC transition matrix to samples at Tth iteration.
"""

def fitpairwise(data, J0, options):
    """
    With input data (binary array of size number of samples x number of neurons), 
    J0 (initial guess for J), and options, fit the pairwise model ot data.
    Return J (learned coupling matrix).
    """
    # Initialize.
    M, n = np.shape(data)
    J0lin = []
    for x in np.nditer(datta):
        J0lin.append(x)
    # Estimate empirical covariance to be reproduced by the model.
    empcov = np.multiple(data.H * data)/M
    # Initialize Gibbs chain at data distribution.
    samples = data[random.randint(1, len(data), 1)]
