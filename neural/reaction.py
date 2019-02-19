import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

"""
Enzymes catalyze biochemical reactions. Many reactions have complex kinematic mechnisms
and specialized qequations are needed to desrcibe their rates in detail. To model a
series of reactions, we use simplified equations.
"""

# Basic differential equation

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

# Second example

# function that returns dy/dt
def model(y,t,k):
    dydt = -k * y
    return dydt

# initial condition
y0 = 5

# time points
t = np.linspace(0,20)

# solve ODEs
k = 0.1
y1 = odeint(model,y0,t,args=(k,))
k = 0.2
y2 = odeint(model,y0,t,args=(k,))
k = 0.5
y3 = odeint(model,y0,t,args=(k,))

# plot results
plt.plot(t,y1,'r-',linewidth=2,label='k=0.1')
plt.plot(t,y2,'b--',linewidth=2,label='k=0.2')
plt.plot(t,y3,'g:',linewidth=2,label='k=0.5')
plt.xlabel('time')
plt.ylabel('y(t)')
plt.legend()
plt.show()

# Michaelis-Menten reaction

# function that returns dy/dt
def mmreaction(v, s, k, t):
    dpdt = (v*s)/(s+k)
    return dpdt

# initial condition
p0 = 0

# time points
t = np.linspace(0,5)

# solve ODE
y = odeint(mmreaction,p0,t)

# plot results
plt.plot(t,p)
plt.xlabel('time')
plt.ylabel('p(t)')
plt.show()

# Increasing Hill function
def hill(ke, x, s, km, l, v, t):
    dpdt = (l/(l+k)*(v*s)/(s+k))
    return dpdt

t = np.linspace(0,20)
y = odeint(hill)



