import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import poisson

"""
Use the Poisson model to simulate neuron firing rates.
"""

dt = 0.001 # Time intervals that we will be looking at
T = 2 # Total simulation time
r = 5  # Firing rate in spikes per second
p = r*dt # Probability of a spike in a very small time interval

n_trials = 500 # Total number of trials

# Run poisson spike neuron simulations
isi = np.zeros(0) # spike intervals
spike_count = np.zeros(n_trials)

for t in range(n_trials):
    """
    Run the random Poisson simulation.
    """
    spikes = np.random.rand(int(T/dt)) # simulation
    spikes[spikes < p] = 1 # convert to spikes
    spikes[spikes < 1] = 0
    spike_count[t] = np.sum(spikes)
    isi = np.append(isi,np.diff(np.where(spikes==1)[0]))

# Plot spike count histograms
print(np.mean(spike_count))
print(np.var(spike_count))
bins = np.arange(np.min(spike_count), np.max(spike_count), 1)

plt.hist(spike_count, bins-0.5, normed = True) # We subtract -0.5 because the values in bins represent the edges
plt.plot(bins, poisson.pmf(bins,mu=r*T),":k", linewidth=2.5)
plt.title("Spike count histogram")
plt.xlabel("Spike count")
plt.legend(["Poisson", "Data"])
plt.show()

# Plot the interspike interval (IVI ivi) distribution
bins = np.arange(np.min(isi), np.max(isi), 10)
plt.hist(isi, bins, normed = True, histtype="barstacked")
plt.plot(bins, r*dt*np.exp(-r*bins*dt),":k", linewidth=3.0) 
plt.title("ISI histogram")
plt.xlabel("ISI interval (ms)")
plt.legend(["Exponential", "Data"])
