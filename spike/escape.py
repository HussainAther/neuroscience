import numpy as np

"""
In the escape noise model, the integrate-and-fire neuron with a threshold (theta) that fires even though the formal
threshold hasn't been reached or may stay quiescent even though the formal threshold has transiently been passed.
We use the "escape rate" or "firing intensity" that depends on the momentary state of the neuron.
This Spike Response Model (SRM) is a type of generalized integrate-and-fire model that we can calculate.
We're going to talk about an escape rate that follows from it.
"""

def kappa(t):
    """
    Response of the membrane to an incoming pulse.
    """
    return np.sin(np.pi*t)

def eta(t):
    """
    Response of the membrane to an outgoing spike.
    """
    return np.cos(np.pi*t)

def SRM(t):
    """
    Membrane potential of an SRM at some point t.
    """
    summ = 0
    Idet = 10 # known deterministic driving current
    kappa =

