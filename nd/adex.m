import brian2 as b2

from neurodynex.adex_model import AdEx
from neurodynex.tools import plot_tools, input_factory

%{
Adaptive Exponential Integrate-and-Fire mode
(integrate and fire)
}%

current = input_factory.get_step_current(10, 250, 1. * b2.ms, 65.0 * b2.pA)
state_monitor, spike_monitor = AdEx.simulate_AdEx_neuron(I_stim=current, simulation_time=400 * b2.ms)
plot_tools.plot_voltage_and_current_traces(state_monitor, current)
print("nr of spikes: {}".format(spike_monitor.count[0]))
# AdEx.plot_adex_state(state_monitor)

# Get a random parameter. provide a random seed to have a reproducible experiment.
random_parameters = LIF.get_random_param_set(random_seed=432)

# Define your test current.
test_current = input_factory.get_step_current(
    t_start=..., t_end=..., unit_time=b2.ms, amplitude= ... * b2.namp)

# Probe the neuron. pass the test current AND the random params to the function.
state_monitor, spike_monitor = LIF.simulate_random_neuron(test_current, random_parameters)

# Plot.
plot_tools.plot_voltage_and_current_traces(state_monitor, test_current, title="experiment")

# Print the parameters to the console and compare with your estimates.
# LIF.print_obfuscated_parameters(random_parameters)
