import time
import numpy as np
from numpy import concatenate as cc
import matplotlib.pyplot as plt

"""
Compute the firing of a neuron via alpha function synapse and a (random) input spike train.

The alpha function is often used for desrcribing synaptic conductance with the expression

P_s = (P_max*t / tau_s) * exp((1-t)/tau_s)

in which P_s is the opening probability of a postsynaptic channel. For an isolated snynapse at time t = 0,
we can generate random spike inputs and compute the membrane voltage using an I & F implementation of dV/dt = - V/RC + I/C.
"""
np.random.seed(123)

h = 1. # step size, Euler method, = dt ms
t_max= 200 # ms, simulation time period
tstop = int(t_max/h) # number of time steps
ref = 0 # refractory period counter

thr = 0.9 # threshold for random spikes
spike_train = np.random.rand(tstop) > thr

def alpha_():
    """
    Alpha function for synaptic conductance
    """
    t_a = 100 # Max duration of syn conductance
    t_peak = 1 # ms
    g_peak = 0.05 # nS (peak synaptic conductance)
    const = g_peak / (t_peak*np.exp(-1));
    t_vec = np.arange(0, t_a + h, h)
    alpha_func = const * t_vec * (np.exp(-t_vec/t_peak))

    plt.plot(t_vec[:80], alpha_func[:80])
    plt.xlabel("t (in ms)")
    plt.title("Alpha Function (Synaptic Conductance for Spike at t=0)")
    plt.draw()


C = 0.5 # capacitance nF
R = 40 # leak resistance M ohms

g_ad = 0 # conductance
G_inc = 1/h # spike rate adaptation 
tau_ad = 2
