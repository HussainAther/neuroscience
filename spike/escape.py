import numpy as np
import scipy as sp

'''
In the escape noise model, the integrate-and-fire neuron with a threshold (theta) that fires even though the formal
threshold hasn't been reached or may stay quiescent even though the formal threshold has transiently been passed.
We use the 'escape rate' or 'firing intensity' that depends on the momentary state of the neuron.
This Spike Response Model (SRM) is a type of generalized integrate-and-fire model that we can calculate.
We're going to talk about an escape rate that follows from it.
'''

def kappa(s):
    '''
    Response of the membrane to an incoming pulse at time s.
    '''
    return np.sin(np.pi*s)

def eta(t):
    '''
    Response of the membrane to an outgoing spike at time t.
    '''
    return np.cos(np.pi*t)

def SRM(t, times):
    '''
    Membrane potential of an SRM at some time point t with past firing times times.
    '''
    summ = 0
    Idet = 10 # known deterministic driving current
    for time in times:
        summ += eta(t-time)
    integrand = lambda s: kappa(s) # intergate the kappa function over all possible vvalues
    return sp.integrate.quad(integrand, 0, np.inf) * Idet

def escaperate(t, times, delay):
    '''
    Measure the spike times with a probability density of the escape rate
    that depends on the momentary distance between the noiseless membrane potential
    and the threshold.
    '''
    return SRM(t, times) - delay

def escapeexp(beta, tau0, t, times, delay):
    '''
    For some parameter beta and time constant tau0, we choose the exponential escape function.
    It fits as beta approaches infinity, the threshold turns into a sharp one so that we return
    to the noiseless model. 
    '''
    return (1/tau0)*np.exp(beta(escape(t, times, delay)))
