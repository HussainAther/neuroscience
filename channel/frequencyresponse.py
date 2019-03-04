import numpy as np
import matplotlib.pyplot as plt

"""
A passive neural membrane is a shift-invariant linear system. According to Eq. 8, the membrane
potential fluctuations (above and below E) may be predicted by convolving the injected current
with a low-pass filter. The parameters of that low-pass filter, g and C, completely determine the
behavior of the membrane.

We use the membrane potential equation dV/dt = - V/RC + I/C. Using the Euler method of derivation,
V(t+h) = V(t) + h*dV/dt .

We can create the charging and discharging curves as follows.
"""

I = 10 # input current (nA)

C = .1 # capacitance (nF)
R = 100 # leak resistance (ohms)
tau = R*C  # = time constant

tstop = 150 # ms
V_inf = I*R # peak V (mV)
tau = 0 # experimental (ms)
h = .2 # ms (step size)
V = 0 # mV
V-trace = [V] # mV

# begin injecting current
for t in np.arange(h, tstop, h):
    V = V +h*(- (V/(R*C)) + (I/C)) # Euler up
    
    # verfiy the membrane time constant if experimetnal conditions are different
    if (not tau and (V > 0.6321*V_inf)):
        tau = t

    # stop injecting current
    if t >= .6*stop:

