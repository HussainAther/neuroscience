import numpy as np
import math

"""
Dynamic causal modelling (DCM) concerns the relation existing between cognitive functions and their 
neurobiological “signature” (the spatio-temporal properties of brain activity) requires an understanding 
of how information is transmitted through brain networks. The ambition here is to ask questions such as: 
"what is the nature of the information that region A passes on to region B"? This stems from the notion 
of functional integration, which views function as an emergent property of brain networks. Dynamic causal 
modelling or DCM was developed specifically to address this question.

See "Dynamic causal modelling revisited" by Friston et al.
"""

# Neuronal parameters
r = 10 # number of regions
p = 15 # number of populations 
i = 10 # number of inputs

thetakappa = np.pi # postsynaptic firing prior
kappai = [256, 128, 16] # postsynaptic rate constant for the i-th neuronal population in N regions

thetaa = np.pi # connectivity prior
a = [[[1 for x in range(p) for y in range(p) for z in range(r)]]] # intrinsic connectivity to population i from population k in each region j

thetab = np.pi # intrinsic connectivity change prior
b = [[[[1 for x in range(p) for y in range(r) for z in range(i) for w in range(r)]]]] # change in intrinsic connectivity by m-th input in region j 

thetaA = np.pi # extrinsic connectivity prior
A = [[[[1 for x in range(p) for y in range(r) for z in range(p) for w in range(r)]]]] # extrinsic connectivity to population i in region j from population k in region l 

thetaB = np.pi # extrinsic connectivity change prior
B = [[[[[1 for x in range(p) for y in range(r) for z in range(p) for w in range(r) for v in range(i)]]]]] # change in extrinsic connectivity by m-th input in region j

thetaC = np.pi # direct driving effect prior
C = [[[1 for x in range(p) for y in range(r) for z in range(i)]]] # direct driving effect of the m-th input on population i in region j 
 
# The parameterization for each of the parameters
params = {"kappai": [], 
          "a": [], 
          "b": [],
          "A": [],
          "B": [],
          "C": []}

for constant in kappai:
    params["kappai"].append(np.exp(thetakappa)*constant
for x, y, z in a:
    params["a"].append(np.exp(thetaa))
for x, y, z, w in b:
    params["b"].append(np.exp(thetab))
for x, y, z, w in A:
    params["A"].append(np.exp(thetaA))
for x, y, z, w, v in B:
    params["B"].append(np.exp(thetaB))
for x, y, z in C:
    params["C"].append(np.exp(thetaC))

# Biophysical parameters
thetanu = np.pi # rate of vasodilatory signal decay prior 
nu = .64 * thetanu # rate of vasodilatory signal decay per second

V0 = .08 # blood volume fraction

k1 = 6.9*psi # intravascular coefficient 


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
