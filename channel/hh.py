import matplotlib.pyplot as plt
import numpy as np

from scipy.integrate import odeint

"""
Hodgkin-Huxley (Hodgkin Huxley) using scipy
"""

# Constants
C_m = 1.0  # uF/cm^2
g_Na = 120.0  # mS/cm^2
g_K = 36.0
g_L = 0.3
V_Na = 50.0  # mV
V_K = -77.0
V_L = -54.402
