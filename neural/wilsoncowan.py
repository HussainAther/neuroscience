import numpy as np
import matplotlib.pyplot as plt

"""
The Wilson-Cowan model describes dynamics of interactions between simple exchitatory
and inhibitory model neurons. It's used in modeling neuronal popuatlions with phase
plane methods and numerical solutions to describe the responses of neuronal populations
to stimuli.
"""

def sig(x)
    """
    Sigmoid function of some value x.
    """
    return 1 / (1 + np.exp(-x))

def fire():
    """
    Each variable represents a different neuron and function in our firing model.
    """
