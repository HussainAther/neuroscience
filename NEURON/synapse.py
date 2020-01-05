import neuron
import numpy
import matplotlib.pyplot as plt

from neuron import h

"""
Excitatory and inhibitory synapsae
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

synapse = h.Exp2Syn(soma(0.5))
synapse.tau1 = 0.5 # [ms]
synapse.tau2 = 10.0 # [ms]
synapse.e = -80.0 

stimulator = h.VecStim()
spike_times = [100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 950.0]
spikes_vector = h.Vector(spike_times)
stimulator.play(spikes_vector)

connection = h.NetCon(stimulator, synapse)
connection.weight[0] = 0.001        # [uS]
