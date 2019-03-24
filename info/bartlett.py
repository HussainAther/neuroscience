import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

"""
In time series analysis, we use Bartlett's method (or the method of averaged periodograms)
to estimate the power spectral density (which gives us the noise spectrum for an input spike train).
"""

def bartlett(k):
    """
    Bartlett's method (or Bartlett's periodogram) uses (1) an original point N data segments into which we split 
    our data k (input array of spikes) into, with each of length m; (2) compute periodogram using the discrete 
    Fourier transform (DFT) then find the squared magnitude of the result and divide by m; (3) average the 
    result of the periodograms above for the data k segments. From this, we get the variance compared to the 
    original N point data segment.
    """
    n = len(k)  # number of segments from our data
