import neuron
import numpy
import matplotlib.pyplot as plt

from neuron import h

"""
Stochastic synapse Tsodyks-Markram Model
"""

# Load external files & initialize.
neuron.h.load_file("stdrun.hoc")
neuron.h.stdinit()
