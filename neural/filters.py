import numpy as np
from np.random import randn
from scipy import signal
import matplotlib.pyplot as plt

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
    SNR = Sinusoid signal-to-noise ratio in dB
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
    """
    
    n = arange(0,N+1) # length N+1
    if not(sqwav):
        x = 1 * np.cos(2 * np.pi * 1/20 * n) # A = 1, Fo/Fs = 1/20
        x += np.sqrt(1 / 2/(10 **(SNR/10)))* randn(N+1)
    else:
        x = 1 * np.sign(np.cos(2*np.pi*1/20*n)); # square wave. A = 1, Fo/Fs = 1/20
        x += np.sqrt(1/1/(10**(SNR/10)))*randn(N+1)

    mu /= M + 1
    y = signal.lfilter([0, 1],1,x)
    x_hat = np.zeros_like(x)
    e = np.zeros_like(x)
    ao = np.zeros(M+1)
    zi = signal.lfiltic(ao, 1, y=0)
    ym = np.zeros_like(ao)

    for k, yk in enumerate(y): # filter
        x_hat[l], zi = signal.lfilter(ap, 1, [yk], zi=zi)
        e[k] = x[k] - x_hat[k]
        ao = ao + 2*mu*e[k]*ym
        ym = np.hstack((array(yk), ym[0:-1]))

    F = np.arange(0,0,5,1/Nfft)
    w, Ao = signal.frqez(ao, 1, 2*np.pi*F)
    Ao = 20 * np.log10(abs(Ao))
    return n, x, x_hat, e, ao, F, Ao

Ntaps = 128
n, x, x_hat, e, ao, F, Ao = lms_ale(10, 1000, Ntaps, 0.01, sqwav=False)

# Plot how the ALE filters noisy input and clean output
plt.plot(n,e**2)
plt.ylabel(r"$e^2[n]$")
plt.xlabel(r"Index $n$")
plt.title(r"Squared Error")
plt.grid()
plt.savefig("ALE_mse.pdf")

# Plot frequency response of the approximately optimum filter
plt.plot(F, Ao)
plt.ylim([-40, 2])
plt.plot([.05, .05], [-40, 0], "r--")
plt.xlabel(r"Normalized Frequency $f/f_s$")
plt.ylabel("r$|W_o(e^{j2\pi f/f_s})|$ (dB)")
plt.title(r"ALE Freq. Response for SNR = 10 dB, $\mu = .01/64$")
plt.grid()
plt.savefig("ALE_fresp.pdf")

"""
With Adaptive Interference Canceling, we implement canceling
in the scikit-dpscomm module sigsys.
"""

def lms_ic(r, M, mu, delta=1):
    """
    Least Mean Square (LMS) interference canceller adaptive filter.
    Complete LMS adaptive filter simulation function for the case
    of interference cancellation.
    
    M : FIR Filter length (order M-1)
    delta : Delay used to generate the reference signal
    mu : LMS step-size
    delta : decorrelation delay between input and FIR filter input
    
    n : ndarray Index vector
    r : ndarray noisy (with interference) input signal
    r_hat : ndarray filtered output (NB_hat[n])
    e : ndarray error sequence (WB_hat[n])
    ao : ndarray final value of weight vector
    F : ndarray frequency response axis vector
    Ao : ndarray frequency response of filter
    """
    N = len(r)-1;
    y = signal.lfilter(np.hstack((np.zeros(delta), np.array([1]))),1,r)
    r_hat = np.zeros(N+1)
    e = np.zeros(N+1)
    ao = np.zeros(M+1)
    z = np.zeros(M)
    ym = np.zeros(M+1)
    for k in range(N+1):
        r_hat[k],z = signal.lfilter(ao, np.array([1]), np.array([y[k]]), zi=z)
        e[k] = r[k] - r_hat[k]
        ao = ao + 2*mu*e[k]*ym
