import numpy as np

"""
We can measure how much one random variable tells us
about another with mutual information.
"""

def P(x, y):
    """
    Joint probability distribution P for x and y.
    """
    return x*y      

def mutI(x, y):
    """
    Mutual information between two arrays of discrete variables x and y
    with a joint probability distribution Pxy.
    """
    summ = 0
    for i in x:
        for j in y:
            summ += P(i, j) 
