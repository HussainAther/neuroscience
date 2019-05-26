import numpy as np
import scipy.signal as signal
import nitime.algorithms as nt_alg
import nitime.utils as nt_ut
import matplotlib.pyplot as plt

"""
Another application of the Slepian functions is to estimate the complex demodulate of a 
narrowband signal. This signal is normally of interest in neuroimaging when finding the 
lowpass power envelope and the instantaneous phase. The traditional technique uses the 
Hilbert transform to find the analytic signal. However, this approach suffers problems of 
bias and reliability, much like the periodogram suffers in PSD estimation. Once again, a 
multi-taper (multitaper) approach can provide an estimate with lower variance.
"""

N = 10000
nfft = np.power( 2, int(np.ceil(np.log2(N))) )
NW = 40
W = float(NW)/N

# Create lowpass band-limited signal
s = np.cumsum(np.random.randn(N))
(b, a) = signal.butter(3, W, btype='lowpass')
slp = signal.lfilter(b, a, s)

# Modulate
s_mod = s * np.cos(2*np.pi*np.arange(N) * float(200) / N)
slp_mod = slp * np.cos(2*np.pi*np.arange(N) * float(200) / N)
fm = int( np.round(float(200) * nfft / N) )

# Create Slepians with the desired badnpass resolution
(dpss, eigs) = nt_alg.dpss_windows(N, NW, 2*NW)
keep = eigs > 0.9
dpss = dpss[keep]; eigs = eigs[keep]

"""
Compare multitaper baseband power estimation with regular Hilbert transform method under
actual narrowband conditions.
"""

# MT method
xk = nt_alg.tapered_spectra(slp_mod, dpss, NFFT=nfft)
mtm_bband = np.sum( 2 * (xk[:,fm] * np.sqrt(eigs))[:,None] * dpss, axis=0 )

# Hilbert transform method
hb_bband = signal.hilbert(slp_mod, N=nfft)[:N]

plt.figure()
plt.subplot(211)
plt.plot(slp_mod, "g")
plt.plot(np.abs(mtm_bband), color="b", linewidth=3)
plt.title("Multitaper Baseband Power")

plt.subplot(212)
plt.plot(slp_mod, "g")
plt.plot(np.abs(hb_bband), color="b', linewidth=3)
plt.title("Hilbert Baseband Power")
plt.gcf().tight_layout()

"""
The Slepian sequences of the multitaper spectral estimation method can also be used to perform a hypothesis 
test regarding the presence of a pure sinusoid at any analyzed frequency. The F-test is used to assess 
whether the power at a given frequency can be attributed to a single line component. In this case, the 
power would be given by the summed spectral convolutions of the Slepian frequency functions with the line 
power spectrum, which is a dirac delta. The complex Fourier coefficient of the putative sinusoid is estimated 
through a linear regression of the Slepian DC components, and the strength of the regression coefficient is 
tested against the residual spectral power for the F-test.
"""

# Set up test with 3 harmonic components with Gaussian noise. 
N = 10000
fft_pow = int( np.ceil(np.log2(N) + 2) )
NW = 4
lines = np.sort(np.random.randint(100, 2**(fft_pow-6), size=(3,)))
while np.any( np.diff(lines) < 2*NW ):
    lines = np.sort(np.random.randint(2**(fft_pow-6), size=(3,)))
lines = lines.astype('d')

