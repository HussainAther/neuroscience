import numpy as np
import matplotlib.pyplot as plt

"""
Most neurons receive their inputs through the synapses. A synapse consists of a presynaptic terminal, synaptic cleft and post-synaptic terminal.
On an arrival of an action potential to a pre-synaptic terminal a neurotransmitter is released, which diffuses through the cleft and activates
directly or through a series of chemical reactions ion channels embedded in the membrane of the post-synaptic terminal. These ion channels selectively
let some ions flow through the membrane, thus increasing the membrane conductance for this ion.

In contrast to the current-based inputs, the response to the conductance change depends on the membrane potential. This leads to several
consequences for the dynamics and the excitability of the neuron, which we will investigate in detail in this tutorial.
"""
# Neuron model
#SI base units
s = 1
kg = 1
m = 1
A = 1

#derived units
S = s**3*A**2/(kg*m**2)
V = kg*m**2*s**-3*A**-1
F = s**4 * A**2 * m**-2 * kg ** -1
Hz = 1/s

#with prefixes
nS = 1e-9 * S
uS = 1e-6 * S
mV = 1e-3 * V
pF = 1e-12 * F
ms = 1e-3 * s

dt = 0.01 * ms # integration time step (ms)
E_e = 0*mV # reversal potential
E_i = -75*mV # reversal potential
E_l = -70*mV # reversal potential
g_l = 1./60*uS # leak conductance
C = 250*pF # membrane capacitance
v_reset = -60*mV # resetting potential
threshold = -50*mV # spike threshold
tau_ref = 2*ms # refractory period

"""
Model neuron with leaky integrate-and-fire model. It uses a one-dimensional differential equation that defines
the evolution of the membrane potential V.
"""

def _lif_update(v, refractory, ge, gi):
    if refractory:
        refractory -= 1
        v = v_reset
        return v, refractory
    v = v + (g_l * (E_l - v) +
             ge * (E_e - v) +
             gi * (E_i - v)) * dt / C
    if v > threshold:
        v = -45 * mV
        refractory = int(tau_ref/dt)
    return v, refractory

def lif_run(g_e, g_i):
    n = np.minimum(len(g_e), len(g_i))
    v = np.zeros(n) * mV
    v[0] = E_l
    refractory = 0
    for i in range(1, n):
        v[i], refractory = _lif_update(v[i-1], refractory, g_e[i], g_i[i])
    return v

"""
Response to a single spike from an excitatory neuron.
Alpha-function shaped post-synpatic conductance.
"""

def alpha_psp(a, tau, tmax=None):
    """Return alpha-function-shaped PSC
    Arguments:
    - a - amplitude
    - tau - time scale
    - tmax - duration (optional)
    """
    if tmax is None:
        tmax = 10 * tau
    time = np.arange(int(tmax / dt)) * dt
    return a * time / tau * np.exp(1 - time / tau)

g_e = np.r_[np.zeros(10*ms/dt), alpha_psp(15*nS, 5*ms, 100*ms)]
g_i = np.zeros(len(g_e))


# vmem - membranne potential
t = np.arange(len(g_e)) * dt
vmem = lif_run(g_e, g_i)
plt.plot(t / ms, vmem / mV)
plt.xlabel('time (ms)')
plt.ylabel('membrane potential (mV)')
plt.axhline(threshold/mV, ls='--', color='k')
plt.text(100, -49.5, 'threshold');

v = -70 * mV
gi_const = 10 * nS
ge_const = (-g_l * (E_l - v) - gi_const * (E_i - v)) / (E_e - v)
vmem_high_conductance = lif_run(g_e + ge_const, g_i + 10*nS)
plt.plot(t / ms, vmem / mV)
plt.plot(t / ms, vmem_high_conductance / mV)
plt.axhline(threshold/mV, ls='--', color='k')
plt.xlabel('time (ms)')
plt.ylabel('membrane potential (ms)');

# Response under synaptic bombardment
tmax = 10 * s   # simulation time (ms)
tau_e = 0.2 * ms  # width of excitatory PSC (ms)
tau_i = 2 * ms    # width of inhibitory PSC (ms)
B_e = 7.1 * nS    # peak excitatory conductance (nS)
B_i = 3.7 * nS    # peak inhibitory conductance (nS)
fr_e = 9655 * Hz    # total firing rate of excitatory population (Hz)
fr_i = 4473 * Hz     # total firing rate of excitatory population (Hz)

epsp = alpha_psp(B_e, tau_e)
ipsp = alpha_psp(B_i, tau_i)


t = np.arange(len(epsp)) * dt
plt.plot(t / ms, epsp / nS, 'r-')
t = np.arange(len(ipsp)) * dt
plt.plot(t / ms, ipsp / nS, 'b-')
plt.xlabel('time (ms)')
plt.ylabel('conductance (nS)');

poisson_spikes = np.random.rand(100*ms/dt) < 100 * Hz * dt
plt.figure(figsize=(10, 0.5))
plt.subplot(111, frameon=False)
plt.plot(poisson_spikes)
plt.xticks([])
plt.yticks([]);

conductance_trace = np.convolve(poisson_spikes, epsp, 'valid')
plt.figure(figsize=(10, 0.5))
plt.subplot(111, frameon=False)
plt.plot(conductance_trace / nS)
plt.xticks([])
plt.yticks([]);

def shot_noise(fr, kernel, tmax, dt):
    poisson_spikes = np.random.rand(tmax/dt) < fr * dt
    shot_noise_trace = np.convolve(poisson_spikes, kernel, 'full')
    return shot_noise_trace[:len(poisson_spikes)]

g_i = shot_noise(fr_i, ipsp, tmax, dt)
g_e = shot_noise(fr_e, epsp, tmax, dt)
t = np.arange(len(g_i)) * dt
plt.plot(t / ms, g_i / nS, 'b')
plt.plot(t / ms, g_e / nS, 'r')
plt.xlim([20, 80])
plt.xlabel('time (ms)')
plt.ylabel('conductance (nS)');

# Simulate LIF neuron
vmem = lif_run(g_e, g_i)

def estimate_firing_rate(vmem):
    return np.mean(vmem>threshold) / dt
def estimate_mean_free_pot(vmem):
    return np.mean(vmem[vmem<threshold])
