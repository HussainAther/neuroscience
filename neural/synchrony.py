import numpy as np

"""
Linear measures of neuronal signal synchrony estimate the synchrony between two or sometimes
more continuous time series of brain activitiy which yield low values for independent time
series and high values for correlated time series. 

We can look at both cross-correlation and coherence measures of this synchrony.
"""

def crosscorrelation(tau, x, n):
    """
    For some time lag tau derived from normalized signals x and y of length N and with zero mean, we
    can calculate the unit variance as linear cross correlation and use that as a measure of synchronization.
    """
    summ = 0
    for i in range(N-tau):
        summ += x[i+tau] * y[n] # account for time lag tau
