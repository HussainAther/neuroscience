import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from gekko import GEKKO

"""
When modeling circuits using cable theory, we assume the segments are cylinders with constant radii.
The electrotonic potential is due to a change in the membrane potential.
At any instant of time, the resting membrane potential (E) can be changed
by several means. We can inject current, use extracellular currents, and
discover changes in membrane conductance caused by driving forces.

As such, the current in the neuron behaves in accordance with mathematical theory.
We also posit that the electrotonic current is ohmic. By this, the passive electrotonic current flow us assuemd to be
in acccordance with the simple linear equation E = IR (current times resisntace), as dictated by Ohm's Law.
"""

# Electrotonic potential measured as a differnece of membrane potential (V_m) and
# electromagnetic potential (E)
def V(V_m, E):
    return V_m - E
