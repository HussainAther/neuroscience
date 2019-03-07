import numpy as np
from scipy.integrate import odeint

"""
McCulloch-Pitts replaces the involved Hodgkin-Huxley system by a threshold device
with only two states (0 and 1) in which 0 denotes the inactivated, silent condition
and 1 denotes the activiated, firing state. We use the equation:

X(t+1) = Theta(I_i - theta_i)

in which t is the discretized time, theta_i is the activation threshold for unit i,
and I_i = X_i are the data that has been identified. In this case,

Postsynpatic current I_iPSC = - summation of j=1 to n of (w_ij * I_j - I-i)

in which w_ij is the synaptic weight of the connection from unit j to unit i,
dependent on the synpatic gain g_ij that evolves during learning and reflecting the synaptic plasticity
and the intracelular resistances. The capacitance is deliberately neglected and I_iext is the
externally controlled input to neuron i.
"""

def simplemodel(w, I):
    """
    Simple conductance model for some neuron i.
    w is a list of weights of connection from unit j to i for each neuron.
    I is a list of currents from each neuron j.
    """
    I_iPSC = 0 # post-synaptic current
    I_ext = 5 # external current to neuron i
    for j in range(len(w)): # simple summation
        I_iPSC += w[j] * I[j]
    I_iPSC -= I_iext # exteranl current
    I_iPSC *= -1 # for the opposite direction
    return I_iPSC

def mp(simplemodel, theta_i):
    """
    McCulloch-Pitts model. theta_i is the activation threshold for unit i.
    """
    X = [] # state of neuron i
    for j in simplemodel:
        if j - theta_i >= 0:
            X.append(1)
        else:
            X.append(0)
    return X

def fn(s):
    """
    FitzHugh and Nagumo approximated the Hodgkin-Huxley equations using a general linear relation (h = a-bn)
    used in combination with coordinate transformation and rescaling to arrive at the Bonhoeffer-Van-der-Pol
    or FitzHugh-Nagumo equations. Takes in s, an array of states of the voltage for each neuron. It must have
    the four states: voltage of first synapse, work of the first synapse, voltage of the second, and work of the second.
    
    Solve for two synapses using ordinary differential equations.
    """
    (v, w, v2, w2) = (state[0], state[1], state[3], state[4])
    
    # these are our constants
    x = 0.08
    y = 0.7
    z = 0.8

    theta = 0 # voltage of synapse
    Vs = 2 # applied voltage
    Iapp = 1.2 # applied current
    gsyn = 30 # synaptic conductance in pS
    S = 1
    lam = -10 # factor to normalize
 
    # Synaptic currents
    Isyn = gsyn*((v-Vs)*S)/(1+np.power(np.e,(lam*(v2-theta))))
    Isyn2 = gsyn*((v2-Vs)*S)/(1+np.power(np.e,(lam*(v-theta))))
    
    # cell 1
    vd = v - power(v, 3)/3 - w + Iapp + Isyn
    wd = x*(v + y - z*w)
      
    # cell 2
    v2d = v2 - power(v2, 3)/3 - w2 + Iapp + Isyn2
    w2d = x*(v2 + y - z*w2)

    # return state derivatives that odeint uses
    return [vd, wd, v2d, w2d]

state0 = ([-1.2, 1.2, -1.2, 1.2])
t = arange(0.0, 2800.0, 0.01)

odeint(fn, state0, t,rtol=1.49012e-13,atol=1.49012e-13)