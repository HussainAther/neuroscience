import numpy as np
import math
import matplotlib.pyplot as plt

from scipy.integrate import odeint
from gekko import GEKKO

# Startup GEKKO for differential equations
m = GEKKO()

"""
The sign of a synaptic action can be predicted by knowledge of the relationship between resting potential (V) and
reversal potential (E). WE can determine amplitude between synpatic conductance and the extra synaptic
conductances. Creating equivalent electrical cirucits we can move from a single-channel conductacne to that of macroscopic
conductances and currents. We can determine the total conductance change produced by simultaneous activation of many
synpases. Gamma is the single-channel conductance, P is the probabilty of opening a single channel, and N is the total
number of ligand-gated channels in the post-synaptic conductacne
"""
gamma = 5
P = 10
N = 15

def total_conductance(gamma, P, N)
    """
    Return the total conductance by simultaneous activation of many synapses.
    """
    g_syn = gamma*P*N # Add conductances in parallel for simultaneous activation
    return g_syn

# Similarly, we can find the current of this setup.
g_syn = 20
V = 20
E = 30
def I_syn(g_syn, V, E):
    return g_syn*(V-E) # From Ohm's law
"""
Closure of the switch lets the membrane change from a potential that we observe
so that we can create a more complete analyitcal decsription of postsynaptic factors
underlying the gneeration of a PSP. This allows us to use the capacitative branch of the circuit,
and by conservation of current, the sum of the branches must equal zero.

We can create a differnetial equation to explain this behavior.
"""
# Set variable values
I_L = 5
I_syn = 5
C = 10

# integration time points
m.time = np.linspace(0,10)

# integration time points
m.time = np.linspace(0,10)

#define Zero parameter to normalize the equation
eq = m.Param(value=0)

# Equation for
m.Equation(C(dVdt) + I_L + I_syn == eq)

#Objective
m.Obj(V = E) # If the ligand-gated channels open by release of a transmitter
# from a presynaptci neuron, this equalizes

# Set global options
m.options.IMODE = 3 # Steady state optimization

# Results
m.solve() # use public server

# Solve
print("C: " + str(C.value))
