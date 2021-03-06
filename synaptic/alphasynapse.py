import numpy as np
import matplotlib.pyplot as plt
import time

from numpy import concatenate as cc

"""
Compute the firing of a neuron via alpha function synapse and a (random) input spike train. The alpha function is often
used for describing synaptic conductance with the expression

P_s = (P_max*t / tau_s) * exp((1-t)/tau_s)

in which P_s is the opening probability of a postsynaptic channel. For an isolated snynapse at time t = 0,
we can generate random spike inputs and compute the membrane voltage using an I & F implementation of dV/dt = - V/RC + I/C.
Since we"re using this I & F implementation, we"ll be characterizing the alpha function with a slightly
different method of finding the conductance.
"""

np.random.seed(123)

h = 1. # step size, Euler method, = dt ms
t_max= 200 # ms, simulation time period
tstop = int(t_max/h) # number of time steps
ref = 0 # refractory period counter

thr = 0.9 # threshold for random spikes
spike_train = np.random.rand(tstop) > thr

def alphasynapse():
    """
    Alpha function for synaptic conductance
    """
    t_a = 100 # Max duration of syn conductance
    t_peak = 1 # ms
    g_peak = 0.05 # nS (peak synaptic conductance)
    const = g_peak / (t_peak*np.exp(-1));
    t_vec = np.arange(0, t_a + h, h) # time used for plotting
    alpha_func = const * t_vec * (np.exp(-t_vec/t_peak)) # calculate alpha function

    plt.plot(t_vec[:80], alpha_func[:80])
    plt.xlabel("t (in ms)")
    plt.title("Alpha Function (Synaptic Conductance for Spike at t=0)")
    plt.draw()


C = 0.5 # capacitance nF
R = 40 # leak resistance M ohms

g_ad = 0 # conductance
G_inc = 1/h # spike rate adaptation
tau_ad = 2

E_leak = -60 # mV, equilibrium potential
E_syn = 0 # Excitatory synapse (why is this excitatory?)
g_syn = 0 # Current syn conductance
V_th = -40 # spike threshold mV
V_spike = 50 # spike value mV
ref_max = 4/h # Starting value of ref period counter

t_list = np.array([], dtype=int)
V = E_leak
V_trace = [V]
t_trace = [0]

fig, axs = plt.subplots(2, 1)
axs[0].plot(np.arange(0,t_max,h), spike_train)
axs[0].set_title("Input spike train")

"""
Perform the simulation on the computed parameters to see how membrane voltage from
the synpatic current vary over time.
"""

for t in range(tstop):
    # Compute input spike train
    if spike_train[t]: # check for input spike
        t_list = cc([t_list, [1]])

    # Calculate synaptic current due to current and past input spikes
    g_syn = np.sum(alpha_func[t_list])
    I_syn = g_syn*(E_syn - V)

    # Update spike times
    if np.any(t_list):
        t_list = t_list + 1
        if t_list[0] == t_a: # Reached max duration of syn conductance
            t_list = t_list[1:]

    # Compute membrane voltage using Euler method: V(t+h) = V(t) + h*dV/dt
    if not ref:
        V = V + h*(-((V-E_leak)*(1+R*g_ad)/(R*C)) + (I_syn/C))
        g_ad = g_ad + h*(-g_ad/tau_ad) # spike rate adaptation
    else:
        ref -= 1
        V = V_th - 10 # reset voltage after spike
        g_ad = 0

    # Generate spikes
    if (V > V_th) and not ref:
        V = V_spike
        ref = ref_max
        g_ad = g_ad + G_inc


    V_trace += [V]
    t_trace += [t*h]

axs[1].plot(t_trace,V_trace)
plt.draw()
axs[1].set_title("Output spike train")
plt.show()
