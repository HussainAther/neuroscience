import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

"""
We can use steady-state current-voltage relationships of axonal sodium and potassium channels
to detect electrochemical properties of neural systems. We use N, the number of channels that can open,
i, the amplitude of current flowing through a single channel, and P, the average channel open probability,
which is a function of voltage and time.
"""

# function that returns dy/dt
def model(y,t):
    k = 0.3
    dydt = -k * y
    return dydt

# initial condition
y0 = 5

# time points
t = np.linspace(0,20)

# solve ODE
y = odeint(model,y0,t)

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()
