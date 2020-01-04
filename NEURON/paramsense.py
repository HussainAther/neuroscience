import matplotlib.pyplot as plt
import neuron
import numpy as np

from simexample import plottv

"""
Parameter sensitivity
"""

for gnabar in [0.1, 0.15]:
    soma.gkbar_hh = 0.01
    soma.gnabar_hh = gnabar

    neuron.h.tstop=30

    neuron.h.run()

    plt.plot(time, max(voltage)*numpy.ones(len(time)), 'r')
    plottv(time, voltage, show=False)

plt.show()
