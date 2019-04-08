import numpy as np

"""
Each action potential at the presynaptic termianl causes the binomially distributed
release of one or more vesicles that pour the total amount of transmitter molecules
into the synaptic cleft. These transmitter molecules react with with ionotropic or with 
metabotropic receptors which open ion channels with a postsynaptic current (PSC).
"""

def IPSC(O, N, g, U, E):
    """
    Postsynaptic current as described above. O is the number of open channels, N
    is the total number of channels, g is the channel conductance at maximum conductnace, 
    U is the synaptic potential, and E is the reversal potential of the G protein-gated
    potassium channels. 
    """ 
    return 
