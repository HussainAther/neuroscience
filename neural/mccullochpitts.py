import numpy as np

"""
McCulloch-Pitts replaces the involved Hodgkin-Huxley system by a threshold device
with only two states (0 and 1) in which 0 denotes the inactivated, silent condition
and 1 denotes the activiated, firing state. We use the equation:

X(t+1) = Theta(I_i - theta_i)

in which t is the discretized time, theta_i is the activation threshold for unit i,
and I_i = X_i are the data that has been identified.
"""
