import matplotlib.pyplot as plt
import neuron
import numpy as np

from neuron import h

"""
Stochastic synapse Tsodyks-Markram Model
"""

# Load external files & initialize.
neuron.h.load_file("stdrun.hoc")
neuron.h.stdinit()

soma = neuron.h.Section()
soma.L = 40
soma.diam = 40
soma.insert("pas")

# Configure the passive biophysics.
for sec in h.allsec():
    sec.Ra = 100
    sec.cm = 1

synapse_list = []
rng_list = []
num_synapses = 10
for i in range(num_synapses):
    synapse = h.StochasticTsodyksMarkram_AMPA_NMDA(soma(0.5))
    rng = h.Random()                                                          
    rng.Random123(1) # Configure the random number generator (rng) type, and the "seed".                     
    rng.uniform(0,1) # Configure the rng to emit uniformly distributed random numbers between 0 and 1
                      # as required by the synapse MOD file.
    synapse.setRNG(rng)
    synapse_list.append(synapse)
    rng_list.append(rng)


