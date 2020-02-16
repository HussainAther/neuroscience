import matplotlib.pyplot as plt
import numpy as np

from scipy.integrate import odeint

"""
Fitzhugh-Nagumo Half-Adder (fitzhugh nagumo) 
"""

# Input current
def input_1(t): return 1 * (t > 500) - 1 * (t>1000) + 1 * (t > 1500)
def input_2(t): return 1 * (t > 1000)

# Constants
theta = gamma = epsilon = 0.1
tmax, m, c = 2000, -100, 60

t = np.arange(0.0, 2000.0, 0.1)
