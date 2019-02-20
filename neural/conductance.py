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
