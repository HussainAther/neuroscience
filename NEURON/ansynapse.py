import neuron
import numpy as np
import matplotlib.pyplot as plt

from neuron import h

"""
AMPA-NMDA synapse model
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

synapse = h.SimpleAMPA_NMDA(soma(0.5))

stimulator = h.VecStim()
spike_times = [100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 950.0]
spikes_vector = h.Vector(spike_times)
stimulator.play(spikes_vector)

connection = h.NetCon(stimulator, synapse)
connection.weight[0] = 1.0 # In units of [nS] due to the gmax scaling factor in our .mod file.

g_syn = h.Vector()
g_syn.record(synapse._ref_g)
i_syn = h.Vector()
i_syn.record(synapse._ref_i)
v_soma = h.Vector()
v_soma.record(soma(0.5)._ref_v)
time = h.Vector()
time.record(neuron.h._ref_t)

h.tstop = 1000.0 # ms
synapse.gmax_NMDA = 0.0
synapse.gmax_AMPA = 0.001 # uS
neuron.h.run()

def plot_timecourse(time_array, dependent_var, newfigure=True, show=True, label=None, ylabel="Membrane voltage (mV)", constants=[]):
    """
    Convenience function to plot time courses of dependent variables.
    """
    if newfigure:
        plt.figure()
    plt.plot(time_array, dependent_var, label=label)
    for constant in constants:
        plt.plot(time_array, constant*np.ones(len(time_array)))
    plt.xlabel("Time (ms)")
    plt.ylabel(ylabel)
    if show:
        plt.show()
       
# Plot. 
plot_timecourse(time, v_soma)
plot_timecourse(time, g_syn, ylabel="Conductance (uS)", label="NEURON")

def dual_exp(t, tau_r, tau_d, t_start):
    """Compute the dual exponential time course using the closed form expression."""
    t = np.array(t)
    time_to_peak = (tau_r*tau_d)/(tau_d-tau_r)*np.log(tau_d/tau_r)
    factor = -np.exp(-time_to_peak/tau_r)+np.exp(-time_to_peak/tau_d)
    f_dual_exp = lambda t: (np.exp(-t/tau_d) - np.exp(-t/tau_r))/factor
    dual_exp = np.zeros_like(t)
    dual_exp[t>=t_start] = f_dual_exp(t[t>=t_start]-t_start)
    return dual_exp
    
plt.plot(time, 0.001*connection.weight[0]*dual_exp(time, synapse.tau_r_AMPA, synapse.tau_d_AMPA, 
                                                   t_start=100.0+connection.delay), "r--", lw=2, label="math. expr.")
plt.legend()

synapse.gmax_NMDA = 0.001 # uS
synapse.mg = 0.0 # mM
synapse.gmax_AMPA = 0 # uS
neuron.h.run()

plot_timecourse(time, g_syn, ylabel="Conductance (uS)", label="NEURON")
plt.plot(time, 0.001*connection.weight[0]*dual_exp(time, synapse.tau_r_NMDA, synapse.tau_d_NMDA, 
                                                   t_start=100.0+connection.delay), "r--", lw=2, label="math. expr.")
plt.legend()

synapse.gmax_AMPA = 0.001 # uS
synapse.gmax_NMDA = 0.7 * 0.001 # uS - 0.7 is a biologically typical ratio of NMDA to AMPA conductance
synapse.mg = 1.0 # mM
g_NMDA = h.Vector()
g_NMDA.record(synapse._ref_g_NMDA)
g_AMPA = h.Vector()
g_AMPA.record(synapse._ref_g_AMPA)
neuron.h.run()

plot_timecourse(time, g_syn, ylabel="Conductance (uS)", label="NEURON - g")
plot_timecourse(time, g_NMDA, ylabel="Conductance (uS)", label="NEURON - g_NMDA", newfigure=False)
plot_timecourse(time, g_AMPA, ylabel="Conductance (uS)", label="NEURON - g_AMPA", newfigure=False)
plt.axis([80.0, 150.0, 0.0, 0.0011])
plt.legend()

g_NMDA_1mM = np.zeros_like(g_NMDA)
g_NMDA_1mM[:] = g_NMDA

plot_timecourse(time, g_NMDA_1mM, ylabel="Conductance (uS)", label="[Mg2+] = 1mM")
mgs = [0.5, 0.25, 0.1, 0.0]
for mg in mgs:
    synapse.mg = mg
    neuron.h.run()
    plot_timecourse(time, g_NMDA, ylabel="Conductance (uS)", label="[Mg2+] = %fmM" % mg, newfigure=False)
plt.axis([80.0, 150.0, 0.0, 0.0011])
plt.legend()

# Now, let's assess the voltage block curve of NMDA for [Mg2+]=1.0 in an _in silico_ reproduction the seminal experiment in Jahr and Stevens, 1990.

# 1) Block the AMPA component of the conductance.
synapse.mg = 1.0 # [mM]
synapse.gmax_AMPA = 0.0 # Apply an "in silico AMPA blocker"
                        # Some things are easy in simulation ...

# 2) Voltage clamp the soma at a given holding voltage.
voltage_clamp = h.VClamp(0.5, sec=soma) # Create a voltage clamp electrode object and place it in the soma.
voltage_clamp.amp[0] = -80.0 # Assign a clamping voltage.
voltage_clamp.dur[0] = h.tstop # Clamp for the whole simulation duration.

# 3) Run the stimulation simulation and extract the peak synaptic conductance.
def extract_peaks(time, trace, event_times, window=10):
    """
    Computes the peak between event_times and returns the times of occurence and the maximums
    Useful for finding PSP or conductance peaks due to synaptic events.
    kwarg "window" defines the time in ms after the event to consider when searching for the peak.
    """
    peaks_list = []
    peaks_times_list = []
    for event_time in event_times:
        i_start = time.searchsorted(event_time)
        i_end = time.searchsorted(event_time+window)
        # Find the index where the max occurs.
        i_max = np.argmax(trace[i_start:i_end])
        # Append the time and value at the max to the respective lists.
        peaks_times_list.append(time[i_start:i_end][i_max])
        peaks_list.append(trace[i_start:i_end][i_max])
    return np.array(peaks_times_list), np.array(peaks_list)
