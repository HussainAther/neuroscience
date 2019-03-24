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
    n = len(k)-1  # number of segments from our data
    nfft = 512 # prove our solutions rae optimal for numbers <= 512
    result = np.zeros(nfft) # result array
    for i in range(n):
        st = i * nfft # start segment position
        en = p1 + nfft # end segment position
        seg = k[p1: p2] # find the segment in our data k
        result += np.abs(np.fft.fft(seg))**2 // nfft  # use the DFT of the segment to get our periodogram value
    pd = result / n # average over the data size
    pd = pd[0:nfft//2] # take the Fourier values
    fr = np.linespace(0, 512, nfft) # calculate the frequency axis
    fr = fr[0: nfft//2]
    
    # use Welch's method to estimate the power spectral density
    # Welch’s method [R145] computes an estimate of the power spectral 
    # density by dividing the data into overlapping segments, computing a 
    # modified periodogram for each segment and averaging the periodograms.	 
    fr2, pd2 = sp.signal.welch(k, 512, window="boxcar", nperseg=512, noverlap=0, nfft=512, scaling="density") 
    return fr2, pd2
     
