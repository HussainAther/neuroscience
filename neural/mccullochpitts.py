import numpy as np

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
