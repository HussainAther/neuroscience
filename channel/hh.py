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

def alpha_m(V): return 0.1 * (V + 40.0) / (1.0 - np.exp(-0.1 * (V + 40.0)))
def beta_m(V):  return 4.0 * np.exp(-0.0556 * (V + 65.0))
def alpha_h(V): return 0.07 * np.exp(-0.05 * (V + 65.0))
def beta_h(V):  return 1.0 / (1.0 + np.exp(-0.1 * (V + 35.0)))
def alpha_n(V): return 0.01 * (V + 55.0) / (1.0 - np.exp(-0.1 * (V + 55.0)))
def beta_n(V):  return 0.125 * np.exp(-0.0125 * (V + 65))

def I_Na(V,m,h): return g_Na * m**3 * h * (V - V_Na)
def I_K(V, n): return g_K * n**4 * (V - V_K)
def I_L(V): return g_L * (V - V_L)
