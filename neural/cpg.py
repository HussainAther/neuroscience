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
    """
