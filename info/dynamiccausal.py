import scipy.sparse as sps
import numpy as np

"""
Dynamic causal modelling (DCM) concerns the relation existing between cognitive functions and their 
neurobiological “signature” (the spatio-temporal properties of brain activity) requires an understanding 
of how information is transmitted through brain networks. The ambition here is to ask questions such as: 
"what is the nature of the information that region A passes on to region B"? This stems from the notion 
of functional integration, which views function as an emergent property of brain networks. Dynamic causal 
modelling or DCM was developed specifically to address this question.
"""

def erp(x, u, P, M):
    """
    Event-related potential for state vector x with the following:
    x[1] - voltage (spiny stellate cells)
    x[2] - voltage (pyramidal cells) +ve
    x[3] - voltage (pyramidal cells) -ve
    x[4] - current (spiny stellate cells)    depolarizing
    x[5] - current (pyramidal cells)         depolarizing
    x[6] - current (pyramidal cells)         hyperpolarizing
    x[7] - voltage (inhibitory interneurons)
    x[8] - current (inhibitory interneurons) depolarizing
    x[9] - voltage (pyramidal cells)
    with P as true connectivity parameters, u as potential energy, and M as the dictionary of
    types of models. 
    Return f (dx(t)/dt = f(x(t))), J (df(t)/dx(t)), and D (delay operator dx(t)/dt)
    """
    n = len(P["A"][0]) # number of sources
    E = np.array([1., 1/2., 1/8.,])*32 # extrinsic rates (forward, backward, lateral)
    G = np.array([1, 4/5., 1/4., 1/4*128.]) # intrinsic rates (g1 g2 g3 g4)    
    D = np.array([2, 16])  #  propogation delays (intrinsic, extrinsic)
    H = np.array([4, 32])  #  receptor densities (excitatory, inhibitory)
    T = np.array([8, 16])  #  synaptic constants (excitatory, inhibitory)
    R = np.array([2, 1])/3  # parameters of static nonlinearity
    try:  # test for free parameters on intrinsic connections
        G = G*np.exp(P["H"])
    except: _
        G = np.ones([n,1])*G
    A = [] # exponential transofrm to ensure positivity constraints
    A[0] = np.exp(P["A"][0])*E[0]
    A[1] = np.exp(P["A"][1])*E[1]    
    A[2] = np.exp(P["A"][2])*E[2]
    C = np.exp(P["C"])
    # Intrinsic connectivity and parameters
    Te = T[0]/1000*np.exp(P["T"][0]) # excitatory time constants
    Ti = T[1]/1000*np.exp(P["T"][1]) # inhibitory time constants
    He = H[0]*np.exp(P["G"][0]) # excitatory receptor density
    Hi = H[1]*np.exp(P["G"][1]) # inhibitory receptor density
 
