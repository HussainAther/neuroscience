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
