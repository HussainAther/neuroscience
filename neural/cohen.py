import numpy as np

"""
Cohen-Grossberg Energy Function (cohen grossberg) used in Hopfield
networks following the Hebb rule consistent with McCulloch-Pitts
model.
"""

def sigmoid(x):
    """
    Sigmoid activation function
    """
    return 1 / (1 + np.exp(-x)) 

def sigmoidderiv(x):
    """
    Derivative of sigmoid function  
    """
    return -(1 + np.exp(-x))**-2 - np.exp(-x) 

def T(u, N):
    """
    Connection strengths between neurons by Hebbian
    outer product rule for N neurons and random binary
    pattern of activity by the condition matrix u that
    details the memory of the model.
    """
    summ = 0
    y, x = u.shape
    for i in range(x):
        summ += (i*g(i))**i
    return (1/N) * summ

def cgef(u, N):
    """
    Cohen and Grossberg converged on the energy function
    from different points of view for N neurons over u matrix of 
    memory conditions.
    """
    summ = 0
         
