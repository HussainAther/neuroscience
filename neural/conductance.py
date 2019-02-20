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
