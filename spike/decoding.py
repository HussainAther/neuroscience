import numpy as np
from scipy.misc import factorial

"""
Estimate static stimulus values on the basis of spike-counting firing rates.
We estimate such a stimulus from the sequence of firing times of the spikes that the stimulus evokes.

If r is the spike-count firing rate of neuron a, we can obtain an estimate of the wind direction
on any given trial from the direction of the vector v_pop.

v_pop = the sum of (r/r_max)*c_a from a=1 to 4.

If we assume each neuron in a popultion fires independently, the firing-rate probability
for the population is the product of the individual probabilities:

"""

def popProb(problist, T, r):
    """
    population probability of firing. problist is the lsit of probabilities
    for each nueron. T is the time period. r is the firing rate.
    """
    result = 1
    for prob in problist:
        factor = exp(-prob*T)
        factor *= (prob*T)**(r*T)
        factor / np.factorial(r*T)
        result *= factor
    return result
"""

The Cramér-Rao bound sets a limit of the variance of any estimate s_est according to:

(sigma_est)^2 >= (1 + b'_est(s))^2 / I_F(s)

in which b'_est(s) is the derivative of b_est(s), the bias on all stimulus variables. I_F(s) is the
Fisher information, a measure of encoding accuracy. Through this, the Fisher infomration limits the
accuracy with which any decoding scheme can extract an estimate of an encoded quantity.
"""

def fisherInformation(xvals, sigmavals):
    """
    xvals is the list of x positions at which the spikes are measured,
    sigmavals is the list of the x position uncertainties.
    """
    npar = len(xvals)
    F = np.zeros([npar,npar])
    for x,sigma in zip(xvals,sigmavals):
        for i in range(npar):
            if i==0:
                dfdpi = x
            else:
                dfdpi = 1
        for j in range(npar):
            if j==0:
                dfdpj = x
            else:
                dfdpj = 1
            F[i,j] += sigma**-2*dfdpi*dfdpj
    return np.mat(F).I
