import numpy as np
import math
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
V = 5 # Voltage
gamma = .5 # conductance
E = 10 # EM potential.

def ionpore(gamma, V, E):
    current = gamma*(V-E)
    return current

"""
We can look at tail currents and max currents to determine the probability of pores open.
"""
I_tail = 10
I_max = 15

def P(I_tail, I_max):
    return I_tail/I_max

"""
The relationship between P_open and the test pulse is the activation curve for Calcium channels. Similar
protocls can be used to detemrmine the voltage-dependence of activation of VGICs using the Boltzmann function.
"""
V_50 = 20
V_m = 14
k = 1e-23

def Pv2(V_50, V_m,k):
    den = (1+exp(V_50 - V_m)/k)
    retrn 1/den

# solve ODE
y = odeint(ionpore,y0,t)

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()

"""
The interdependence of g, V, and I when VGICs are activated posed a significant challenge to
early biophysiicsts in their quest to deterine the ionic basis of suprathreshold changes in the
membrane potential. Biophysicists wanted to observe adn quantify the voltage and time-dependent changes in
conductance of the plasma membrane to sodium and potassium ions (g_Na and g_K) they believed gave rise
to the action potential.
"""



def g_Na(I_Na, V_m, E_Na):
    conductance = I_Na / (V_m - E_Na)
    return conductance

def g_K(I_K, V_m, E_K):
    conductance = I_K / (V_m - E_K)
    return conductance

