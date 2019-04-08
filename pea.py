import numpy np

"""
We use the partition ensemble average (PEA pea) framework in reducing a system.
At each step we introduce a transition rate from nucleation to the release
of an MFE (multiple-firing event) to capture the dynamic nature of their emergence. 

We use a conductance-based integrate-and-fire (IF if) model as they're more realistic
versions of the Hodgkin-Huxley type of neuronal model. We introduce an algorithm
to describe this IF network that has homogeneity and total synchrony.
"""

def peaalg(ni, ne, mit, mey, siy, sey, sii, sie, sei, see):
    """
    For network sizes of inhibitory and excitatory (ni and ne) neurons, feedfoward input rates (miy, mey) and input
    strengths (siy, sey), and network synaptic coupling strengths (sii, sie, sei, see) we perform an algorithm for our
    PEA method.
    """ 
