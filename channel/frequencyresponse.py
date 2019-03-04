import numpy as np

"""
A passive neural membrane is a shift-invariant linear system. According to Eq. 8, the membrane
potential fluctuations (above and below E) may be predicted by convolving the injected current
with a low-pass filter. The parameters of that low-pass filter, g and C, completely determine the
behavior of the membrane. We use the membrane potential equation dV/dt = - V/RC + I/C.

We can create the charging and discharging curves as follows.
"""

I = 10 # input current (nA)

C = .1 # capacitance (nF)
R = 100 # leak resistance (ohms)
tau = R*C  # = time constant

tstop = 150 # ms
V_inf = I*R
tau = 0
