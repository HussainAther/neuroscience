import numpy as np

"""
We makea a simple Python simulation of Adaptive Line Enhancement using a single
sinusoid at normalized frequency plus additive white Gaussian noise. We assume
we have a narrowband signal buried in broadband additive noise. For statistically
stationary inptus, we use a Wiener filter.

The term adaptive filter implies changing the characteristic of a filter in some automated
fashion to obtain the best possible signal quality in spite of changing signal/system conditions.
Adaptive filters are usually associated with the broader topic of statistical signal processing.

The operation of signal filtering by definition implies extracting something desired from a signal
containing both desired and undesired components. With linear FIR and IIR filters the filter output is
obtained as a linear function of the observation (signal applied) to the input.

An optimum linear filter in the minimum mean square sense can be designed to extract a signal from noise by
minimizing the error signal formed by subtracting the filtered signal from the desired signal. For noisy
signals with time varying statistics, this minimization process is often done using an adaptive filter.
"""


def lms_ale(SNR,N,M,mu,sqwav=False,Nfft=1024):
    """
    Least Mean Squares Adaptive
    Line Enhancement Algorithm using an IIR filter.
    n,x,x_hat,e,ao,F,Ao = lms_ale(SNR,N,M,mu)
    *******LMS ALE Simulation************
    SNR = Sinusoid SNR in dB
    N = Number of simulation samples
    M = FIR Filter length (order M-1)
    mu = LMS step-size
    mode = 0 <=> sinusoid, 1 <=> squarewave

    n = Index vector
    x = Noisy input
    x_hat = Filtered output
    e = Error sequence
    ao = Final value of weight vector
    F = Frequency response axis vector
    Ao = Frequency response of filter in dB
    **************************************

    n = arange(0,N+1) # length N+1
    if not(sqwav):
        x = 1*cos(2*pi*1/20*n) # A = 1, Fo/Fs = 1/20
        x += sqrt(1/2/(10**(SNR/10)))*randn(N+1)
    else:
        x = 1*sign(cos(2*pi*1/20*n)); # square wave. A = 1, Fo/Fs = 1/20

    
