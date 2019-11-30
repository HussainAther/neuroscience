import numpy as np

"""
Short-term (short term) memory equations.
"""

def additive(x, i):
    """
    In the additive model, we add the terms, possibly nonlinear, that
    determine the rate of change of neuronal activites, or potentials, x,
    with a single potential term i. Return the rate of change for that
    particular potential.
    """
    pd = -5*x[i] # Passive decay term using a nonlinear signal
              # with a weight value (5 here).
    pf = 0 # Positive feedback term summed over the individual 
           # potentials
    for i in range(len(x)):
        pf +=  
