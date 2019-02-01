import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

"""
We can use steady-state current-voltage relationships of axonal sodium and potassium channels
to detect electrochemical properties of neural systems. We use N, the number of channels that can open,
i, the amplitude of current flowing through a single channel, and P, the average channel open probability,
which is a function of voltage and time.
"""

# function that returns current
def current(N,P,i):
    I_calcium = N*P*i # Calculate the calcium current as a function of N, P, and i
    return I_calcium

"""
We get the steady state P from macroscopic currents by measuring the other variables and holding them
constant to observe the relationship between them. Returning membrane potential to a hyperpolarized value between
depolarizaiton can allow VGICs to recover from inactivation and return to a closed state.

"""

# initial condition
y0 = 5

# time points
t = np.linspace(0,20)

# solve ODE
y = odeint(current,y0,t)

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()


"""
We can use the current product of single channel conductance (gamma) and the driving force (V-E) between the voltage
and the electromagnetic force to get the direction of current flow through a single channel. The reversing
potential for calcium channels is because they permit the flow fo cations outwards through the central
ion pore.
"""
def ionpore(gamma, V, E):
    current = gamma*(V-E)
    return current


"""
We can look at tail currents and max currents to determine the probability of pores open.
"""
def P(I_tail, I_max):
    return I_tail/Imax

# solve ODE
y = odeint(ionpore,y0,t)

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()
