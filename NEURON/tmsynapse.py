import neuron
import numpy
import matplotlib.pyplot as plt

from neuron import h

# Load external files & initialize
neuron.h.load_file("stdrun.hoc")
neuron.h.stdinit()

"""
Tsodyks-Markram (tsodyks markram) model.
"""
