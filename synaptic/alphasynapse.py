import time
import numpy as np
from numpy import concatenate as cc
import matplotlib.pyplot as plt

"""
Compute the firing of a neuron via alpha function synapse and a (random) input spike train.

The alpha function is often used for desrcribing synaptic conductance with the expression

P_s = (P_max*t / tau_s) * exp((1-t)/tau_s)

in which P_s is the opening probability of a postsynaptic channel.

for an isolated snynapse at time t = 0.
"""
