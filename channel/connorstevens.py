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

def I_m():
    """
    Membrane current for Connor-Stevens Model.
    """
    g_L =.003
    g_Na =
