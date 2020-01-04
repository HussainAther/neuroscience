import neuron

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
