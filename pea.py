import numpy as np

from scipy.stats import bernoulli

"""
From "A coarse-graining framework for spiking neuronal networks: from strongly-
coupled conductance-based integrate-and-fire neurons to augmented systems of ODEs"
by Zhang, et al.

We use the partition ensemble average (PEA pea) framework in reducing a system.
At each step we introduce a transition rate from nucleation to the release
of an MFE (multiple-firing event) to capture the dynamic nature of their emergence. 

We use a conductance-based integrate-and-fire (IF if) model as they're more realistic
versions of the Hodgkin-Huxley type of neuronal model. We introduce an algorithm
to describe this IF network that has homogeneity and total synchrony.
"""

def mastereq():
    """
    Master equaton to evolve a probability.
    """

def me(i):
    """
    Single-neuron firing rate.
    """
    return np.sin(i)

def peaalg(ni, ne, mit, mey, siy, sey, sii, sie, sei, see):
    """
    For network sizes of inhibitory and excitatory (ni and ne) neurons, feedfoward input rates (miy, mey) and input
    strengths (siy, sey), and network synaptic coupling strengths (sii, sie, sei, see) we perform an algorithm for our
    PEA method. We can calculate PEA for MFE neurons with a probability from the Bernoulli random variable pi_k over the
    time interval [kdt, kdt+dt].
    """
    rho_vt = (miy * siy)/ni + (mey * sey)/ne # probability of a single neuron spiking dependent upon velocity and time
    for i in range(10): # time step range
        pi_k = bernoulli.rvs(rho_vt, size=1) # if MFE has occured: 0 or 1
        if pi_k == 0:
            rho_vt = mastereq(rho_vt, k) 
            k += 1
        pmfe = .1 * ne * me(k*.1) * (1 - q(k*i))
