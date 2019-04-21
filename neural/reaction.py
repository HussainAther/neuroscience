import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import odeint

"""
Enzymes catalyze biochemical reactions. Many reactions have complex kinematic mechnisms
and specialized qequations are needed to desrcibe their rates in detail. To model a
series of reactions, we use simplified equations.
"""

def model(y,t):
    """
    Differential equation of the reaction rate.
    """
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

"""
With the Michaelis-Menten reaction, we model enzyme kinematics.
"""

def mmreaction(v, s, k, t):
    """
    For a given maximum rate of the system (v) at saturating substrate concentration,
    substrate concentration s, and Michaelis constant k, we determine the rate of reaction. 
    """
    dpdt = (v*s)/(s+k)
    return dpdt

p0 = 0 # initial condition
t = np.linspace(0,5) # time points

# solve ODE
y = odeint(mmreaction,p0,t)

# plot results
plt.plot(t,p)
plt.xlabel("time")
plt.ylabel("p(t)")
plt.show()

"""
We use the Hill-Langmuir equation as a special case of a rectangular hyperboola
to determine the fraction of the receptor protein concentration bound by the ligand.
"""

def hill(ke, x, s, km, l, v, t):
    """
    With the same variabels as above, we introduce l as the free, unbound ligand
    concentration, ke as the dissociation constant derived from the law of mass action.
    """
    dpdt = (l/(l+k)*(v*s)/(s+k))
    return dpdt

t = np.linspace(0,20)
y = odeint(hill)



