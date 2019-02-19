import itertools
import numpy as np
import scipy as sp

"""
When we observe spike train of a sensory neuron, we learn about many different aspects of the stimulus.
This way, we gain information. Shannon entropy is closely related to information. We can think about information transmission by sensory neurons in terms of
various parts of the spike train as part of stimulus space.


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

def integrand(p):
    """
    The inside of the integral for the integral notion of entropy.
    """
    return p * np.log(p))

def entropyInt(variableProbs, k):
    """
    Integral notion of entropy.
    """
    S = sp.integrate(integrand(prob for prob in variableProbs), 0, 1)
    return S *= -k

def totalEntropy(
