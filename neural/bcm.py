import numpy as np

"""
During neural development, sensory stimulation induces long-term changes in the receptive 
field of the neurons that encode the stimuli. We use hte Bienenstock-Cooper-Munro (BCM bcm) 
model to analyze this process computationallly. We model the BCM plasticity model for a neuron
with N synapses. 
"""

def bcm(x_i, y, theta, tau_w, tau_theta):
    """
    We create an (N+1)-dimensional system with two equations. We assume N is even, use input strength
    at input i that is x_i, the neuron's output rate y, and we use a linear relation between input
    and output y = summation of i for (w_i*x_i). w(t) is the synaptic weight vector. 
    """
    dwdt = x_i*y(y-theta)/tau_w
    dthetadt = -(theta+y**2) /tau_theta 
    """
    1/tau_theta gives us the threshold update rate, and 1/tau_w is the learning rate.
    """
    return dwdt, dthetadt

"""
In contrast to Hebbian learning, which depends only on 2nd-order statitistics and fails to converge
in whitened environments, BCM depends on higher-order statistics and can find structure in nearly
any environment. With Hebbian, more input of uncorrelated noise into the closed eye for deprivation
slows the deprivation effect. With BCM, more of this input increases the deprivation effect. Hebbian learning
requires an asymmetry in either the neural architecture or the temporal structure to get direction
selectivity, not simultaneously. With BCM, spatio-temporal correlations on the order of the receptive
field size are required to achieve both strong orientation and direction selectivity.
"""
