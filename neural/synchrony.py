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
    The absolute value is symmetric in the x and y directions and reaches a maxixmum of 1 for complete lag 
    synchronization.
    """
    if tau < 0:
        tau = -tau # flip the sign so we're looking at time lag in the same direction every time we run this function
    summ = 0
    for i in range(N-tau):
        summ += x[i+tau] * y[n] # account for time lag tau
    return (1/(N-tau)) * summ

def crossspec(w):
    """
    Return cross-spectrum of some signal w by performing Discrete Fourier transform as the x direction and its inverse for the y.
    """
    return np.fft.fft(w)*np.fft.ifft(w) # multiply the discrete fourier transform of the signal by its complex conjugate.

def coherence(w):
    """
    We can quantify linear correlations in the frequency domain with the cross spectrum for some signal w.
    """
    cross = crossspec(w)
    """
    If we normalize the amplitdue of the power spectrum of both the x and y systems, the we can calculate the 
    cross-spectrum as the coherence function:
    """
    num = abs(crossspec(w))**2 # numerator
    den = abs(np.fft.fft(w)*np.fft.fft(w)) * abs(np.fft.ifft(w) * np.fft.ifft(w)) 
    return num / den
