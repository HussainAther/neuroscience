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
spiketimes = [100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 950.0]
spikesvector = h.Vector(spiketimes)
stimulator.play(spikesvector)

connection = h.NetCon(stimulator, synapse)
connection.weight[0] = 0.001 # [uS]

g_syn = h.Vector()
g_syn.record(synapse._ref_g)
i_syn = h.Vector()
i_syn.record(synapse._ref_i)
v_soma = h.Vector()
v_soma.record(soma(0.5)._ref_v)
time = h.Vector()
time.record(neuron.h._ref_t)

h.tstop = 1000.0 # ms
neuron.h.run()

def plottimecourse(time_array, dependent_var, newfigure=True, show=True, label=None, ylabel="Membrane voltage (mV)", constants=[]):
    """
    Convenience function to plot time courses of dependent variables.
    """
    if newfigure:
        plt.figure()
    plt.plot(time_array, dependent_var, label=label)
    for constant in constants:
        plt.plot(time_array, constant*numpy.ones(len(time_array)))
    plt.xlabel('Time (ms)')
    plt.ylabel(ylabel)
    if show:
        plt.show()
