import numpy as np
from neuronpy.graphics import spikeplot
import math
import random

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
"""

def nextTime(rateParameter):
    return -math.log(1.0 - random.random()) / rateParameter
