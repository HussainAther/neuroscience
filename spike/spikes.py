import numpy as np
import math
import random

from neuronpy.graphics import spikeplot

"""
Simple, approximate descriptions of spike statistics. They can be useful for analytic statements about the
qimplications of whole classes of such models for idfferent interesting quantities such as the reliability of code
or the possibility of decoding the spike train.

We can begin with a Poisson model. One spike firing with some probability per unit time can depend on time
but no explicitly on the occurrence times of the other spikes.
"""

# Poisson theta values: these are decided such that the model is Poisson-like
spikes = []
num_cells = 10
num_spikes_per_cell = 20
frequency = 20

for i in range(num_cells):
    isi = np.random.poisson(frequency, num_spikes_per_cell)
    spikes.append(numpy.cumsum(isi))

# spikes is now a list of lists where each cell has a list of spike
# times. Now, let's plot these spikes.
sp = spikeplot.SpikePlot()
sp.plot_spikes(spikes)

"""
Real spike trains can't be exactly Poisson processes. Spikes can't come too close together in time because the
spike-generating mechanism of all cells is refractory.

With a stochastic process of N spikes over T time period, we can create a distribution of event.
"""

N = 1000.0
T = 2.0

count = int(1E6)
x = np.arange(count)
y = -np.log(1.0 - np.random.random_sample(len(x))) / T
np.average(y)

plt.hist(y, 10)
plt.show()

# Intervals over large sum
sum([random.expovariate(T) for i in range(count)])/count

# Distribute time events
intervals = [expovariate(lmbda) for i in range(1000)]
timestamps = [0.0]
timestamp = 0.0
for t in intervals:
    timestamp += t
    timestamps.append(timestamp)
timestamps[:10]

deltas = [y - x for x, y in zip(timestamps, timestamps[1:])]
plt.figure(figsize=(16, 4))
plt.plot(deltas, 'r+')
plt.show()
