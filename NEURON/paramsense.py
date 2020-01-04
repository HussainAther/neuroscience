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

for with_dend in [False, True]:
    neuron.h.tstop = 100

    soma.gkbar_hh = 0.01
    soma.gnabar_hh = 0.1

    if with_dend:
        dend.connect(soma)
    else:
        neuron.h.disconnect(sec=dend) # disconnect dend for now
        
    neuron.h.run()

    plottv(time, voltage, show=False, label="with dend" if with_dend else "without dend")

plt.legend()
plt.show()

dend_ra = 100
dend_gl = 5e-4

for with_dend in [False, True]:
    # For every addition of mechanism create figure show newer model
    # Lines for reversal potentials Na, K and leak
    neuron.h.tstop = 100

    soma.gkbar_hh = 0.01
    soma.gnabar_hh = 0.1

    dend.el_hh = -65 # Reversal potential leak current, mV
    dend.gl_hh = dend_gl # Leak conductance, S/cm^2
    dend.Ra = dend_ra
    
    if with_dend:
        dend.connect(soma)
    else:
        neuron.h.disconnect(sec=dend) # disconnect dend for now
        
    neuron.h.run()

    # Convert the NEURON vectors to numpy arrays.
    time_py = time.to_python()
    voltage_py = voltage.to_python()

    plottv(time_py, voltage_py, show=False, label="with dend" if with_dend else "without dend")

plt.legend()
plt.show()
