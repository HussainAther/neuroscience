import numpy as np
import pylab as plt

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
    
    # constants used for fitting to general linear relation
    x = 0.08
    y = 0.7
    z = 0.8

    theta = 0 # voltage phase shift
    Vs = 2 # applied voltage
    Iapp = 1.2 # applied current
    gsyn = 30 # synaptic conductance in pS
 
    # Synaptic currents
    Isyn = gsyn*((v-Vs))/(1+np.power(np.e,(v2-theta)))
    Isyn2 = gsyn*((v2-Vs))/(1+np.power(np.e,(v-theta)))
    
    # synapse 1
    vd = v - np.power(v, 3)/3 - w + Iapp + Isyn
    wd = x*(v + y - z*w)
      
    # synapse 2
    v2d = v2 - np.power(v2, 3)/3 - w2 + Iapp + Isyn2
    w2d = x*(v2 + y - z*w2)

    # return state derivatives that odeint uses
    return [vd, wd, v2d, w2d]

s = ([-1.2, 1.2, -1.2, 1.2])
t = np.arange(0.0, 2800.0, 0.01)

odeint(fn, s, t, rtol=1.49012e-13, atol=1.49012e-13)

"""
Morris-Lecar model described spiking dynamics of potassium- and calcium-controlled muscle fibers.
"""

# Constants
C_m = 1.0 # membrane capacitance, in uF/cm^2
g_Ca = 1.1 # maximum conducances, in mS/cm^2
g_K = 2.0
g_L = 0.5
E_Ca = 100.0 # Nernst reversal potentials, in mV
E_K = -70.0
E_L = -50.0

def m_infty(V):
    """
    Membrane voltage derived from Fourier transform of the derivative of the signal. 
    Returns the open-state probability function of the open state of the channel. They're
    partitioned according to a Boltzmann distribution. 
    """
    return (1.0 + sp.tanh((V + 1.0) / 15.0)) / 2.0

def w_infty(V):
    """
    Same but for the closed state of the channel.
    """
    return (1.0 + sp.tanh(V / 30.0)) / 2.0

def tau_w(V):
    """
    Memnbrane time constant.
    """
    return 5.0 / sp.cosh(V / 60.0)  # in ms

def I_ext(t):
    """
    Input voltage.
    """
    return 10*sp.floor(t/100)

t = sp.arange(0.0, 400.0, 0.1)
I = I_ext(t)

def ml(V, w,t):
    """
    Morris-Lecar using odeint again. V and w are initial conditions.
    """
    V, w = X # initial conditions for
    dVdt = (I_ext(t) - I_Ca(V) - I_K(V, w) - I_L(V)) / C_m
    dwdt = (w_infty(V) - w) / tau_w(V)
    return dVdt, dwdt

X = odeint(ml(-44, .05) t)
