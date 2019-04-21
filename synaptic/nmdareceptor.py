import numpy as np

"""
The NMDA (nmda) receptor conductance effects the postsynaptic potnetial not normally
seen in other conductances. The NMDA receptor current uses a factor that
depends on postsynpatic potential V. We write it as g_NMDA * G_NMDA(V) * P(V-E_NMDA) in which
the factor G_NMDA(V) describes an extra voltage dependence due ot the fact that when the
postsynaptic neuron is near its resting poetntial, NMDA receptors are blocked by Mg2+ ions.

To activate the condutctance, the postsynaptic neuron must be depolarized to knock out the blocking ions.
According to Jahr and Stevens (1990) paper, the dependence can be described by the following function:
"""

def G_NMDA(Mg2, V):
    """
    Compute G_NMDA factor due to the way NMDA receptors conduct Ca2+ ions as well as monovalent cations.
    Mg2 is the concentration of Mg2+ (Magnesium ion) and V is the membrane potential. 
    """
    return (1 + ([Mg2]/3.57) * np.exp(-V/16.13))**-1
