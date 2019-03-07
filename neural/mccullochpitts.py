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
dependent on the synaptci gain g_ij that evolves during learning and reflecting the synaptic plasticity
and the intracelular resistances. The capacitance is deliberately neglected and I_i is the
externally controlled input to neuron i.
"""

def mp():
    """
    McCulloch-Pitts model.
    """
    I_iPSC =  # postsynpatic current
