import numpy as np

"""
Steady-state (steady state) of a correlation detector for periodic sine gratings
moving at a constant velocity. For sine waves of wavelength lambda (lamd) and contrast
deltaI traveling at velocity v, the inputs are spaced by deltaphi, the temporal filter
has an amplitude and phase response A and phi, respectively. 
"""

def R(A, phi, deltaI, deltaphi, v, lambd, lowpass=False):
    """
    Tme-averaged detector response R.
    When the sinusoidal waves in steady state have an additional amplitude factor A(omega),
    and phase shift phi(omega), we can multiply the respective signals and subtract the result of
    the left and right multiplier to get the time-averaged detector response for amplitude A,
    phase response phi, contrast deltaI, input spacing deltaphi, velocity v, wavelength lambd, and 
    boolean of whether we're assuming it's a low-pass filter. 
    """
    tau = 1 # time constant
    if lowpass == True: # If the temporal filter is a low-pass of first order
        return deltaI**2 * ((tau*2*np.pi*v/lamd)/(1+tau**2*(2*np.pi*v/lamd)**2)) * np.sin(2*np.pi*deltaphi/lambd)
    return deltaI**2 * A*(2*np.pi*v/lamd) * np.sin(-phi*(2*np.pi*v/lamd)) * np.sin(2*np.pi*deltaphi/lamd)
