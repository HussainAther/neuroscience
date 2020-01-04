import matplotlib.pyplot as plt
import neuron
import numpy

"""
Simulation example
"""

neuron.h.load_file("stdrun.hoc")

neuron.h.stdinit()

soma = neuron.h.Section(name='soma')

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

neuron.h.run()

def plot_tv(time_array, voltage_array, show=True, label=None, constants=[]):
    plt.plot(time_array, voltage_array, label=label)
    for constant in constants:
        plt.plot(time_array, constant*numpy.ones(len(time_array)))
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane voltage (mV)')
    if show:
        plt.show()
    
plot_tv(time, voltage)
