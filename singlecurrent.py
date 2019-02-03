import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt

"""
The sign of a synaptic action can be predicted by knowledge of the relationship between resting potential (V) and
reversal potential (E). WE can determine amplitude between synpatic conductance and the extra synaptic
conductances. Creating equivalent electrical cirucits we can move from a single-channel conductacne to that of macroscopic
conductances and currents. We can determine the total conductance change produced by simultaneous activvation of many
synpases. Gamma is the single-channel conductance, P is the probabilty of opening a single channel, and N is the total
number of ligand-gated channels in the post-synaptic conductacne
"""
gamma = 5
P = 10
N = 15
# function that returns total conductance
def total_conductance(gamma, P, N)
    g_syn = gamma*P*N # Add conductances in parallel for simultaneous activation
    return g_syn

# Similarly, we can find the urrent of this setup.
g_syn = 20
V = 20
E = 30
def I(g_syn, V, E):
    return g_syn*(V-E)


