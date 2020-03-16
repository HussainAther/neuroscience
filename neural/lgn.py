import numpy as np

"""
Lateral geniculate nucleus (LGN) neurons for ocular dominance,
the way cells receive more input fomr one eye.

In ocular dominance, the two lateral geinculate nuclei contribute 
to the two weights wL and wR (one weight from each) of the visual cortex. 
"""

eps = .5 # the correlation between activity in the two populations
        # is epsilon (eps) times as strong as the mean squared activity
        # in each population

C = np.array([[1, eps], 
              [eps, 1]])

"""
Find eigenvectors of the matrix.
"""

vec = C*np.array([1, 1]).T 
assert vec == (1+eps)*np.array([1, 1]).T

vec2 = C*np.array([1, -1]).T
assert vec2 == (1-eps)*np.array([1, -1]).T
# eigenvalues 1+eps and 1-eps 

"""
First, the eigendirection [1, 1].T represents the sum of wL and wR, i.e. 
this direction represents the total synaptic weight onto our V1 neuron. 
The eigendirection [1, −1].T represents the difference between wL and wR, 
with positive values meaning that the contribution from the left eye is 
dominant and negative value meaning that the right eye is dominant. 

A reasonable definition of ocular dominance (OD) is the difference between 
the two weights normalized by the sum, i.e. OD = wL − wR/(wL + wR). 
The fact that the two eigenvalues are perpendicular means that these 
two components – the sum of the weights and their difference – 
develop independently. 

For example, the growth in the difference between wL and wR depends only 
on the current difference – the difference grows just the same if wL = 1 
and wR = 2 or if wL = 101 and wR = 102. Likewise, the total weight grows 
the same whether wL = 10 and wR = 10 or wL =1 and wR = 19.
""" 
