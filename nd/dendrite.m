import brian2 as b2
import matplotlib.pyplot as plt

from neurodynex.cable_equation import passive_cable
from neurodynex.tools import input_factory

%{
Dendrites and passive cable equation
%}

passive_cable.getting_started()
CABLE_LENGTH = 500. * b2.um  % length of dendrite
CABLE_DIAMETER = 2. * b2.um  % diameter of dendrite
R_LONGITUDINAL = 0.5 * b2.kohm * b2.mm  % Intracellular medium resistance
R_TRANSVERSAL = 1.25 * b2.Mohm * b2.mm ** 2  % cell membrane resistance (->leak current)
E_LEAK = -70. * b2.mV  % reversal potential of the leak current (-> resting potential)
CAPACITANCE = 0.8 * b2.uF / b2.cm ** 2  % membrane capacitance

% spatial indexing
voltage_monitor, cable_model = passive_cable.simulate_passive_cable(...)
probe_location = 0.123 * b2.mm
v = voltage_monitor[cable_model.morphology[probe_location]].v

t_spikes = [10, 15, 20]
l_spikes = [100. * b2.um, 200. * b2.um, 300. * b2.um]
current = input_factory.get_spikes_current(t_spikes, 100*b2.us, 0.8*b2.namp, append_zero=True)
voltage_monitor_ABC, cable_model = passive_cable.simulate_passive_cable(..., current_injection_location=l_spikes, input_current=current, ...)
