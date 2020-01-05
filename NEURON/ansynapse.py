import neuron
import numpy
import matplotlib.pyplot as plt

from neuron import h

"""
AMPA-NMDA synapse model
"""

# Load external files & initialize
neuron.h.load_file("stdrun.hoc")
neuron.h.stdinit()
