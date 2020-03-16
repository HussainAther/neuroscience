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

spike_times = [100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 950.0]

conn_list = []
stimulator = h.VecStim()
spikes_vector = h.Vector(spike_times)
stimulator.play(spikes_vector)

for synapse in synapse_list:
    connection = h.NetCon(stimulator, synapse)
    connection.weight[0] = 1.0        # In units of [nS] due to the gmax scaling factor in our .mod file
    conn_list.append(connection)

g_syn_list = []
Use_syn_list = []
R_syn_list = []
for synapse in synapse_list:
    g_syn = h.Vector()
    g_syn.record(synapse._ref_g)
    g_syn_list.append(g_syn)
    Use_syn = h.Vector()
    Use_syn.record(synapse._ref_Use)
    Use_syn_list.append(Use_syn)
    R_syn = h.Vector()
    R_syn.record(synapse._ref_R)
    R_syn_list.append(R_syn)
v_soma = h.Vector()
v_soma.record(soma(0.5)._ref_v)
time = h.Vector()
time.record(neuron.h._ref_t)

for synapse in synapse_list:
    synapse.gmax_AMPA = 0.001 # uS
    synapse.gmax_NMDA = 0.7 * 0.001 # uS - 0.7 is a biologically typical ratio of NMDA to AMPA conductance
    synapse.mg = 1.0 # mM
    synapse.U1 = 0.5 # Baseline release probability
    synapse.tau_rec = 700 # ms - recovery from depression
    synapse.tau_facil = 10 # ms - relaxation from facilitation

h.tstop = 1000.0 # ms
neuron.h.run()

def plot_timecourse(time_array, dependent_var, newfigure=True, show=True, label=None, ylabel="Membrane voltage (mV)", constants=[]):
    """Convenience function to plot time courses of dependent variables"""
    if newfigure:
        plt.figure()
    plt.plot(time_array, dependent_var, label=label)
    for constant in constants:
        plt.plot(time_array, constant*np.ones(len(time_array)))
    plt.xlabel("Time (ms)")
    plt.ylabel(ylabel)
    if show:
        plt.show()
        
plot_timecourse(time, v_soma)
plt.axis([0, 1000, -70, -60])

def plot_synapse_traces():
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.collections import PolyCollection
    
    fig = plt.figure()
    ax = fig.gca(projection="3d")

    # use the hsv colormap (https://matplotlib.org/users/colormaps.html)
    colormap = plt.get_cmap("hsv")

    verts = []
    for i, g_syn in enumerate(g_syn_list):
        verts.append(list(zip(time, 1000.*np.array(g_syn))))

    # produce 10 different colors with 60% transparency (alpha) using
    # list comprehension (https://docs.python.org/2/tutorial/datastructures.html#list-comprehensions)
    facecolors = [colormap(x, alpha=0.6) for x in np.linspace(0,1,10)]
    poly = PolyCollection(verts, facecolors=facecolors, edgecolors=facecolors)
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=range(num_synapses), zdir="y")

    ax.set_xlabel("time [ms]")
    ax.set_xlim3d(0, 1000)
    ax.set_zlabel("conductance [nS]")
    ax.set_zlim3d(0, 1.)
    ax.set_ylabel("synapse #")
    ax.set_ylim3d(0, num_synapses)
    ax.view_init(30, -80)

plot_synapse_traces()

for i in range(len(rng_list)):
    rng_list[i].Random123(i) # seed each rng with its index in the rng_list
    rng_list[i].uniform(0,1)

h.tstop = 1000.0 # ms
neuron.h.run()
plot_synapse_traces()

plot_timecourse(time, np.mean(g_syn_list, axis=0), ylabel="conductance [us]")
plt.axis([0, 1000, 0, 0.001])
plot_timecourse(time, v_soma)
plt.axis([0, 1000, -70, -60])

mean_gsyn_list = []
mean_Rsyn_list = []
mean_Usesyn_list = []
v_list = []
num_trials = 100
for i in range(num_trials):
    neuron.h.run()
    v_list.append(np.array(v_soma))
    mean_gsyn_list.append(np.mean(g_syn_list, axis=0))
    mean_Rsyn_list.append(np.mean(R_syn_list, axis=0))
    mean_Usesyn_list.append(np.mean(Use_syn_list, axis=0))

plt.figure()
for v in v_list:
    plt.plot(time, v, "-", color="0.7")
plt.plot(time, np.mean(v_list, axis=0), "b-", lw=2)
plt.axis([50, 1000, -70, -60])
plt.ylabel("voltage [mV]")
plt.xlabel("time [ms]")

def extract_peaks(time, trace, event_times, window=10):
    """
    Computes the peak between event_times and returns the times of occurence and the maximums
    Useful for finding PSP or conductance peaks due to synaptic events.
    kwarg "window" defines the time in ms after the event to consider when searching for the peak
    """
    
    peaks_list = []
    peaks_times_list = []
    for event_time in event_times:
        i_start = time.searchsorted(event_time)
        i_end = time.searchsorted(event_time+window)
        # find the index where the max occurs
        i_max = np.argmax(trace[i_start:i_end])
        # append the time and value at the max to the respective lists
        peaks_times_list.append(time[i_start:i_end][i_max])
        peaks_list.append(trace[i_start:i_end][i_max])
        
    return np.array(peaks_times_list), np.array(peaks_list)

peak_times, peaks = extract_peaks(np.array(time), np.mean(v_list, axis=0), spike_times)
plt.figure()
for v in v_list:
    plt.plot(time, v, "-", color="0.7")
plt.plot(time, np.mean(v_list, axis=0), "b-", lw=2)
plt.plot(peak_times, peaks, "r.", ms=5)
plt.axis([50, 1000, -70, -60])
plt.ylabel("voltage [mV]")
plt.xlabel("time [ms]")
psps = []
for v in v_list:
    peak_times, peaks = extract_peaks(np.array(time), v, spike_times)
    psps.append(peaks)
# turn it into an array so we can take column slices, to compute the histograms
psps = np.vstack(psps) - np.min(v)
plt.figure()
bins = np.linspace(0., 10., 50)
plt.hist(psps[:,0], bins=bins, facecolor="blue", alpha=0.5, label="PSP dist of $1^\mathrm{st}$ spike")  # 1st spike is the 0th column
plt.hist(psps[:,7], bins=bins, facecolor="red", alpha=0.5, label="PSP dist of $8^\mathrm{th}$ spike")  # 1st spike is the 0th column
