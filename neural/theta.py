import numpy as np

"""
The theta model, or Ermentrout–Kopell canonical model, is a biological 
neuron model originally developed to model neurons in the animal Aplysia, 
and later used in various fields of computational neuroscience. The model 
is particularly well suited to describe neuron bursting, which are rapid 
oscillations in the membrane potential of a neuron interrupted by periods of relatively little oscillation.
"""

def I(t):
    """
    Adjust current input over time t.
    """
    return 0 

def dthetadt(theta, t):
    """
    dtheta/dt = 1 - cos(theta) + (1 + cos(theta)) I(t)
    """
    return 1 - np.cos(theta) + (1 + np.cos(theta)) * I(t)
