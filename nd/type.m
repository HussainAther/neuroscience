import brian2 as b2
import matplotlib.pyplot as plti
import numpy as np

from neurodynex.tools import input_factory, plot_tools, spike_tools
from neurodynex.neuron_type import neurons

%{
Type I and type II neuron models
}%

# Create an input current.
input_current = input_factory.get_step_current(50, 150, 1.*b2.ms, 0.5*b2.pA)
