import brian2 as b2
import matplotlib.pyplot as plt
import numpy as np

from neurodynex.hodgkin_huxley import HH
from neurodynex.tools import input_factory

%{
Hodgkin-Huxley (hodgkin huxley) model of squid axon.
}%

HH.getting_started()

current = input_factory.get_step_current(5, 100, b2.ms, I_min *b2.uA)
state_monitor = HH.simulate_HH_neuron(current, 120 * b2.ms)
HH.plot_data(state_monitor, title="HH Neuron, minimal current")
