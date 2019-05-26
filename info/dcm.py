import numpy as np
import math

"""
Dynamic causal modelling (DCM) concerns the relation existing between cognitive functions and their 
neurobiological “signature” (the spatio-temporal properties of brain activity) requires an understanding 
of how information is transmitted through brain networks. The ambition here is to ask questions such as: 
"what is the nature of the information that region A passes on to region B"? This stems from the notion 
of functional integration, which views function as an emergent property of brain networks. Dynamic causal 
modelling or DCM was developed specifically to address this question.
"""

def deltaz(z, u, theta):
    """
    Change in neural activity z with respect to time. Experimental inputs u and parameters theta that
    control the neural function.
    """

def y(z, theta, eps):
    """
    Timeseries y generated from observation function with parameters theta, neural activity z, and additive
    noise eps (epsilon). 
    """

"""
Neural model in DCM for fMRI is a Taylor approximation that captures the gross causal influences between
brain regions and their change due to experimental inputs.
"""

# Taylor expansion
x = 2
e_to_2 = 0
for i in range(5):
    e_to_2 += x**i/math.factorial(i)
