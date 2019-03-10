
"""
We makea a simple Python simulation of Adaptive Line Enhancement using a single
sinusoid at normalized frequency plus additive white Gaussian noise. We assume
we have a narrowband signal buried in broadband additive noise. For statistically
stationary inptus, we use a Wiener filter.
"""


def lms_ale(SNR,N,M,mu,sqwav=False,Nfft=1024):
 """
 lms_ale lms ALE adaptation algorithm using an IIR filter.
 n,x,x_hat,e,ao,F,Ao = lms_ale(SNR,N,M,mu)
