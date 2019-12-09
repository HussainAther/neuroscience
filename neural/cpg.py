import numpy as np

"""
Central pattern generators (cpg) neurons for cyclic task control
in biological systems create an oscillatory activity.
"""

# Membrane dynamics
dVdt = (Isyn(V, t) + Ipm(V, t) + Iinj - V(t)/R)/C
