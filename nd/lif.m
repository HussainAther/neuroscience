import brian2 as b2
import matplotlib.pyplot as plt
import numpy as np

from neurodynex.leaky_integrate_and_fire import LIF
from neurodynex.tools import input_factory, plot_tools

%{
Leaky-integrate-and-fire (leaky integrate fire)
%}

LIF.getting_started()
LIF.print_default_parameters()

V_REST = -70*b2.mV
V_RESET = -65*b2.mV
FIRING_THRESHOLD = -50*b2.mV
MEMBRANE_RESISTANCE = 10. * b2.Mohm
MEMBRANE_TIME_SCALE = 8. * b2.ms
ABSOLUTE_REFRACTORY_PERIOD = 2.0 * b2.ms

% Create a step current with amplitude= i_min.
step_current = input_factory.get_step_current(
    t_start=5, t_end=100, unit_time=b2.ms,
    amplitude= i_min)  % set i_min to your value

% Run the LIF model.
% Note: As we do not specify any model parameters, the simulation runs with the default values
(state_monitor,spike_monitor) = LIF.simulate_LIF_neuron(input_current=step_current, simulation_time = 100 * b2.ms)

% Plot I and vm.
plot_tools.plot_voltage_and_current_traces(
state_monitor, step_current, title="min input", firing_threshold=LIF.FIRING_THRESHOLD)
print("nr of spikes: {}".format(spike_monitor.count[0])) % Should be 0.

% Get a random parameter. Provide a random seed to have a reproducible experiment.
random_parameters = LIF.get_random_param_set(random_seed=432)

% Define your test current.
test_current = input_factory.get_step_current(
    t_start=..., t_end=..., unit_time=b2.ms, amplitude= ... * b2.namp)

% Probe the neuron. Pass the test current AND the random params to the function.
state_monitor, spike_monitor = LIF.simulate_random_neuron(test_current, random_parameters)

% Plot.
plot_tools.plot_voltage_and_current_traces(state_monitor, test_current, title="experiment")

% Print the parameters to the console and compare with your estimates.
% LIF.print_obfuscated_parameters(random_parameters)

% Note the higher resolution when discretizing the sine wave: we specify unit_time=0.1 * b2.ms
sinusoidal_current = input_factory.get_sinusoidal_current(200, 1000, unit_time=0.1 * b2.ms,
                                            amplitude= 2.5 * b2.namp, frequency=250*b2.Hz,
                                            direct_current=0. * b2.namp)

% Run the LIF model. By setting the firing threshold to to a high value, we make sure to stay in the linear (non spiking) regime.
(state_monitor, spike_monitor) = LIF.simulate_LIF_neuron(input_current=sinusoidal_current, simulation_time = 120 * b2.ms, firing_threshold=0*b2.mV)

% Plot the membrane voltage.
plot_tools.plot_voltage_and_current_traces(state_monitor, sinusoidal_current, title="Sinusoidal input current")
print("nr of spikes: {}".format(spike_monitor.count[0]))
