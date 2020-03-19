import numpy as np

"""
Membrane potential increases monotonically in the suprathreshold regime.
"""

# Given a positive derivative, the current term in the rate equation is
def Vt(I, R, phi, tau_m):
    """
    Current term V(t) depends on mean input current and reacts to changes
    in presynaptic spike rate on the time scale of hte synaptic time constant
    tau_syn.
    """
