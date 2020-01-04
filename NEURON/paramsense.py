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

soma.gkbar_hh = 0.01

# definitely mention critical value where state changes
# show something with percentages
# show value we were using before on plot

max_voltages = []
import numpy
gnabar_range = numpy.arange(.05, 0.2, 0.001)
for gnabar in gnabar_range:
    soma.gnabar_hh = gnabar

    neuron.h.run()

    max_voltages.append(max(voltage))

plt.plot(gnabar_range, max_voltages, 'oC0')
plt.xlabel("gnabar (S/cm2)")
plt.ylabel("Maximum AP voltage")
for xs in [0.1, 0.15]:
    plt.axvline(x=xs, color="r")
plt.show()

# Extend model with dendrite.
dend = neuron.h.Section(name='dend')
dend.connect(soma)
dend.L = 400 # micron
dend.diam = 2.0 # micron
dend.nseg = 9 # number of segments in the dendritic section

dend.insert("hh")
dend.el_hh = -65 # Reversal potential leak current, mV
dend.gl_hh = 5e-4 # Leak conductance, S/cm^2

dend.gkbar_hh = 0.0
dend.gnabar_hh = 0.0
