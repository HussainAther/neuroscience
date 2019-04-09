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
            summ += P(i, j) * np.log(P(i, j)/i*j) 
    return summ
   
def H(x,y):
    """
    Uncertainty that (1) should be maximal when Px(x) is uniform and should
    increase with the number of possible values x can take, (2) should be 
    same if we reorder probabilities assigned to different values of x,
    (3) uncertainty about two independent random variables should be the 
    sum of hte uncertainties about each of them. 
    """
    Px = np.zeros(len(x)) # marginals
    Py = np.zeros(len(y))
    Hx = 0 # uncertainty of random variable x
    Hy = 0 # of y
    for i in range(len(y)):
        Px += P(i, y[i])
        Hx -= Px * np.log(P(i, y[i]))
    for i in range(len(x)):
        Py += P(i, x[i])
        Hy -= Py * np.log(P(i, x[i]))
    return Hx, Hy  
