import itertools
import numpy as np
import scipy as sp

"""
When we observe spike train of a sensory neuron, we learn about many different aspects of the stimulus.
This way, we gain information. Shannon entropy is closely related to information. We can think about information transmission by sensory neurons in terms of
various parts of the spike train as part of stimulus space.

In classical mechanics, particles can take on a continuous range of positions and velocities. There's no natural scale for the precisino of the measurements.
"""

def totalProb(variableProbs):
    """
    The additivity of information mean that for any distribution of N independent variables,
    the probability of each of them occuring is the product of each of the individual variables.
    """
    return itertools.product(variableProbs)

def entropySum(variableProbs, k):
    """
    Summation notion of entropy as logarithm of the number of possible states the systme can occupy.
    k is a constant.
    """
    S = 0
    for prob in variableProbs:
        S += (prob * np.log(prob))
    S *= -k
    return S

def entropyIntIntegrand(p):
    """
    The inside of the integral for the integral notion of entropy.
    """
    return p * np.log(p))

def entropyInt(variableProbs, k):
    """
    Integral notion of entropy.
    """
    S = sp.integrate(entropyIntIntegrand(prob for prob in variableProbs), 0, 1)
    return S *= -k

def entropyDiscrete(deltaV, sigma):
    """
    If we measure voltage V in milliVolts then we can take hte continuous voltage variable and place it in discrete bins of size
    ΔV. This equation works if ΔV is very small. sigma is the variance.
    """
    return (1/2)*np.log2(2*np.pi*exp((sigma**2/(deltaV**2)))

def entropyConditionalIntegrand(P):
    """
    Entropy conditional integrand
    """
    return P*np.log2(P))
    
def entropyConditional(conditionalProbs):
    """
    Conditional entropy uses a conditional distribution of probabilities that measure rrelative
    likelihood of different input signals x given that we have observed a particular output value y.
    """
    return -sp.integrate(entropyConditionalIntegrand(prob for prob in conditionalProbs), 0, 1)
    
def probReadout(y):
    """
    Readout y is porportional to s with some gain g. This is the mutual information of the Gaussian channel. si s a signal and we have an s-detector
    that gives us a readout y with some gain g. y is a list of values of which we calculate the arviacne for this function.
    """
    return (1/np.sqrt(2 * np.pi * np.var(y)**2)) * exp(-y**2)/(2*np.var(y**2))

"""
To use these signals in a biological context we need to generalize a bit and consider signals that vary in time.
We use time-dependent signals. We want to describe a function of time f(t) and confine our attentino to a time interval of size T with
0 < t < T. We can use a fourier series, sum of sine and cosine functions. They later become helpful in constructing power spectral density.
"""
def cn(n):
   c = y*np.exp(-1j*2*n*np.pi*time/period)
   return c.sum()/c.size

def fourier(x, Nh):
   f = np.array([2*cn(i)*np.exp(1j*2*i*np.pi*x/period) for i in range(1,Nh+1)])
   return f.sum()
