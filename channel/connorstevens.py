import numpy as np

"""
The Hodgkin-Huxley model was developed on the basis of data from the giant squid axon.
The Connor-Stevens model provides an alternate descirption of action-potential general using
fast Na+, delayed rectifier K+, and leakaage conductsances (as the HH model does), but the Na+ and K+
conductances have properties somewhat different. The C-S model als ohas an extra K+ conductance
(know nas the A-current) that is transient.
"""

dt = 0.000005
tmax = 2

iclamp_flag = 1 # if this is 1, run under current clamp conditions
vclamp_flag = 0 # otherwise this should be 1, for voltage clamp conditions

istart = 0.5 # time applied current starts
ilength = 1 # length of applied current pulse
Ie=2.405e-7 # magnitude of applied current pulse
Ie = 0.78e-7 # threshold for constant spiking with no A-current

vstart = 0.25 # time to step voltage
vlength = 0.5 # length of voltage step
V0 = -0.080 # initial voltage before step
Ve = -0.040 # value of votage stepped to


def I_m():
    """
    Membrane current for Connor-Stevens Model.
    """
    g_L =.003
    g_Na =
