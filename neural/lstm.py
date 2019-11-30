import numpy as np

"""
Long short-term memory (long short term lstm) model
of a recurrent neural network.
"""

def forgetlayer(x):
    """
    Sigmoid function to determine if we forget.
    """
    return 1/(1+np.exp(-x))

def candidatelayer(x):
    """
    Lol.
    """
    return np.tanh(x)   

def inputlayer(x):
    """
    Pass the input gate?
    """ 
    return np.sin(x)

ct = [0, 0, 0] # candidate layer
ht = [0, 0, 0] # hidden layer

def lstmcell(prevct, prevht, input):
    """
    For prevct previous candidate layer, prevht previous hidden layer,
    and input input, update!
    """
    combine = prevht + input
    ft = forgetlayer(combine)
    candidate = candidatelayer(combine)
    it = inputlayer(combine)
     
