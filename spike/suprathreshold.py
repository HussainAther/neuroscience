import numpy as np

"""
Membrane potential increases monotonically in the suprathreshold regime.
"""

# Given a positive derivative, the current term in the rate equation is
def Vt(I, R, phi, tau_m):
    """
    Current term V(t) depends on mean input current and reacts to changes
    in presynaptic spike rate on the time scale of hte synaptic time constant
    tau_syn. We use an ensemble avergae of the current values I, spike threshold phi,
    membrane time constant tau_m, and resistance R.
    """
    return (np.average(I)*R-phi)/tau_m

def rate(P, I, R, phi, tau_m):
    """
    With voltage term P(t) and the current term V(t) (defined above), return
    the PSTH (peristimulus time histogram) rate.
    """
    return rate*Vt(I, R, phi, tau_m)
