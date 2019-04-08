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
    return O*g*(U-E)/N

def theta(t):
    """
    Postsynaptic potential (PSP) we use ise Heaviside (heaviside) jump function.
    """
    if t <=0:
        return 0
    return 1

def sigma(t):
    """
    Dirac delta function at some time t.
    """
    if t == 0:
        return 1
    return 0

def UPSC(t, t1, I0, tau, C_m):
    """
    We can easily compute the postsynaptic potential evoked by an arbitrary spike train.
    We solve the following equation using Green's function (U_PSP = G(t, t')):
    tau*dUdt + U = - I_PSC / C_m .
    with times t, t1, capacitcance C_m, time constant tau, and initial current I0.  
    """
    G = theta(t-t1) * I0 * np.exp(-t-t1/tau) / C_m
    I = I0*sigma(t-t1)
    return np.convolve(G, I)
     
