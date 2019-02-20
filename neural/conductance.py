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
