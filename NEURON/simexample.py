import matplotlib.pyplot as plt
import neuron
import numpy as np

"""
Simulation example
"""

neuron.h.load_file("stdrun.hoc")

neuron.h.stdinit()

soma = neuron.h.Section(name="soma")

print("Soma object:", soma)
print("Soma object name: ", soma.name())
print("Number of segments in the soma:", soma.nseg)

soma.L = 40
soma.diam = 40
print("Soma length: %f micron" % soma.L)
print("Soma diameter: %f micron" % soma.diam)

soma_area_eq = 2 * neuron.h.PI * soma.L * soma.diam / 2
print "Soma area according to cylinder surface area equation: %f micron^2" % soma_area_eq

# The 0.5 refers to the segment in the middle of the soma
# Because there is only one segment, in this case it refers to the entire soma
soma_area = neuron.h.area(0.5, sec=soma)
print("Soma area according to NEURON: %f micron^2" % soma_area)
print("Both values match: %s" % (soma_area_eq == soma_area))

soma_sphere_area_eq = 4 * neuron.h.PI * pow(soma.diam / 2, 2)
print("Soma area according to sphere surface area equation: %f micron^2" % soma_sphere_area_eq)
print("Specific capacitance: %f uf/cm2" % soma.cm)

soma_tcap = (soma.cm * (soma_area / pow(1e4, 2)))
print("Total soma capacitance: %f uf" % (soma.cm * (soma_area / pow(1e4, 2))))

# Run simulation.
print("Membrane voltage soma: %f mV" % soma(.5).v) # mV
print("Current time: %f ms" % neuron.h.t) # ms

neuron.h.tstop = 100
print("Simulation stop time: %f ms" % neuron.h.tstop)
print("Integration time step: %f ms" % neuron.h.dt)

time = neuron.h.Vector()
voltage = neuron.h.Vector()

time.record(neuron.h._ref_t)
voltage.record(soma(.5)._ref_v)

# Run the simulation.
neuron.h.run()

def plottv(time_array, voltage_array, show=True, label=None, constants=[]):
    plt.plot(time_array, voltage_array, label=label)
    for constant in constants:
        plt.plot(time_array, constant*np.ones(len(time_array)))
    plt.xlabel("Time (ms)")
    plt.ylabel("Membrane voltage (mV)")
    if show:
        plt.show()
plottv(time, voltage)

# Inject current.
iclamp = neuron.h.IClamp(.5, sec=soma)
iclamp.amp = 0.1 # nA
iclamp.delay = 10 # ms
iclamp.dur = 50 # ms
neuron.h.run()
plot_tv(time, voltage)

V1 = -65 # Voltage before stimulus, mV
V2 = soma.v # Voltage after stimulus, mV
deltaV = V2 - V1 # Voltage difference, mV
Im = iclamp.amp # nA
deltaT = iclamp.dur # ms
soma_tcap # total soma membrane capacitance, uF

deltaV_eq = Im * deltaT / soma_tcap # in nA * ms / uF == microvolt
deltaV_eq /= 1e3 # Correction factor to get mV

print("Observed dV: %f mV" % deltaV)
print("Calculated dV: %f mV" % deltaV_eq)
print("Simulated dV matches equation dV: %s" % (deltaV - deltaV_eq < 1e-6))

# Add leak conductance.
soma.insert("hh")

# Set conductance to zero for Na and K.
soma.gkbar_hh = 0.0
soma.gnabar_hh = 0.0

soma.gl_hh = 5e-4 # Leak conductance, S/cm^2

el = soma.el_hh = -65 # Reversal potential leak current, mV
print("Reveral of leak current: %f mV" % el)

soma.gkbar_hh = 0.0
soma.gnabar_hh = 0.0

neuron.h.tstop = 100

neuron.h.run()

plottv(time, voltage)

soma.gkbar_hh = 0.01 # in S/cm^2
soma.gnabar_hh = 0.1

neuron.h.run()

plottv(time, voltage)

ena = soma.ena # explain e_na
ek = soma.ek
el = soma.el_hh

print("Na reversal: %f" % ena)
print("K reversal: %f" % ek)
print("Leak reversal: %f" % el)

plottv(time, voltage, constants=[ena, ek, el])

linearcoef = np.polyfit(gkbar_range, maxvoltages, 1)

print("Linear equation maxvoltage = %f*gkbar + %f" % tuple([x for x in linearcoef]))
