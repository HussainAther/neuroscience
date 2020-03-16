import numpy as np

"""
Central pattern generators (cpg) neurons for cyclic task control
in biological systems create an oscillatory activity.
"""

# Membrane dynamics
def dVdt(V, t):
    """
    Based on the synaptic current (Isyn), resting membrane potential current
    (Ipm), and injected current (Iinj) alongside voltage (V), time (t), membrane
    resistance (R) and capacitance (C).
    """
    return (Isyn(V, t) + Ipm(V, t) + Iinj - V(t)/R)/C

# Output function (firing rate)
def gamma(V, Vo, t):
    """
    Using gamma we can model how the neuron itself fires using a form
    of the Poisson equation with an initial resting potential (Vo),
    time-dependent voltage (V), and inverted electrons-volts (m). 
    """
    return 1/(1+np.exp(-2)*m(V(t) - Vo))
