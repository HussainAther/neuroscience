import brian2 as b2
import matplotlib.pyplot as plti
import numpy as np

from neurodynex.tools import input_factory, plot_tools, spike_tools
from neurodynex.neuron_type import neurons

%{
Type I and type II neuron models
%}

% Create an input current.
input_current = input_factory.get_step_current(50, 150, 1.*b2.ms, 0.5*b2.pA)

% Get one instance of class NeuronX and save that object in the variable 'a_neuron_of_type_X'.
a_neuron_of_type_X = neurons.NeuronX()  % we do not know if it's type I or II

% Simulate it and get the state variables.
state_monitor = a_neuron_of_type_X.run(input_current, 200*b2.ms)

% Plot state vs. time.
neurons.plot_data(state_monitor, title='Neuron of Type X')

% Get an instance of class NeuronY.
a_neuron_of_type_Y = neurons.NeuronY()  % we do not know if it's type I or II
state_monitor = a_neuron_of_type_Y.run(input_current, 200*b2.ms)
neurons.plot_data(state_monitor, title='Neuron of Type Y')
