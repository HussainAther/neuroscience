import numpy as np

"""
Long short-term memory (long short term lstm) model
of a recurrent neural network.
"""

ct = [0, 0, 0] # candidate layer
ht = [0, 0, 0] # hidden layer

def lstmcell(prevct, prevht, input):
    """
    For prevct previous candidate layer, prevht previous hidden layer,
    and input input, update!
    """
