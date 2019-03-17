import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from gekko import GEKKO

# Startup GEKKO for differential equations
m = GEKKO()

"""
When modeling circuits using cable theory, we assume the segments are cylinders with constant radii.
The electrotonic potential is due to a change in the membrane potential.
At any instant of time, the resting membrane potential (E) can be changed
by several means. We can inject current, use extracellular currents, and
discover changes in membrane conductance caused by driving forces.

As such, the current in the neuron behaves in accordance with mathematical theory.
We also posit that the electrotonic current is ohmic. By this, the passive electrotonic current flow us assuemd to be
in acccordance with the simple linear equation E = IR (current times resisntace), as dictated by Ohm's Law.
"""

# Electrotonic potential measured as a differnece of membrane potential (V_m) and
# electromagnetic potential (E)
def V(V_m, E):
    return V_m - E

"""
In the steady-state, we ignore membrane capacitance and usually the resting membrane potential.
In the simplest case, the electrotonic potential (V) is relative to a uniform resting potential (E).
Electrotonic current divides between internal and membrane resistance. Axial resistance is inversely
proportional to diameter while membrane resistance is inversely proportional to membrane surface area.
The external medium along the process is assumed to have zero resistivity.
We assume driving forces on membrane conductances are constant, and cables have different boundary conditions

With these assumptions and rules, we rerpesent internal resistance (r_i) connected to the r_i of neighboring segments and
through the membrane resistance (r_m) to the ground. We use differential equations to
describe the spread of electrotonic potential under steady-state conditions.
Electrotonic spreads depends on the characteristci length constant.
"""

# Set variable values
r_i = 5
r_m = 10
V = [m.Var(value = np.sin(2*xpos[i])) for i in range(npx)]

# Discretization of space
xpos = np.linspace(0,10,1)
dx = xpos[1]−xpos[0]

# Equation
m.Equation(V == (r_m/r_i) * (V[0].dt() == c**2))

# Solve it!
m.solve()


"""
When constructing a compartmental model of the passive electrical properties of a dendritic branch,
we need to understand how a endritic segments interacts with its organelles. From an abstract model,
we can look at the membrane capacitance (c_m), membrane resistance (r_m), resting membrane potentail (E),
and internal resistance (r_i). If we have a steady-state current input at point x = 0, the electrotonic potential (V)
along the cable is proportional to the second derivative of the potential (d^2V) with respect to the ratio of membrane resistance (r_m) to
the internal resistance (r_i) and distance to this membrane. We obtain an exponential euler relationship.
"""

# Set variable values
x = 5
lmbda = 10 # Square root of (r_m/r_i)
V = [m.Var(value = np.sin(2*xpos[i])) for i in range(npx)]

# Discretization of space
xpos = np.linspace(0,10,1)
dx = xpos[1]−xpos[0]

# Equation
m.Equation(V == V_0 * exp((-x)/lmbda))

# Solve it!
m.solve()

"""
When x = lambda, the ration of V to V_0 is exp(-1) or about .37. Lambda is the characteristc
length constant of the cable at this point. The decay of membrane potential along an infinte dendritic
cable is described by the length constant.

Current flows within a neuron due to voltage gradients. At any point along a cable of radius a and intracellular
resistivity r_L, the longitudinal current I_L flowing in the direction of increasing x is as follows.
"""
a = 5
r_L = 3

# Equation
m.Equation(I_L == -(np.pi*a**2/r_L) * (V.dx()))

# Solve it!
m.solve()

"""
We determine the membrane potential V(x,t) by solving a partial differential equation of the cable equaiton
that describes how the currents entering, leaving, adn flowing within a neuron affect the rate of change of the
membrane potential.

We can create the second-order partial differential linearized cable equation by grouping the constants together
and expressing tau_m (time constant) and lambda (length constant, or sqrt((a*r_m)/(2*r_L)) ) as our constants.

tau_m * dv/dt = lambda^2 * d^2v/dx^2 - v + r_m * i_e

in which i_e is the inhibitory current.
"""

tau_m = 5
lmbda = 10 # lenght constant: indicates how far a stationary current will influence the voltage along the cable
r_m = 4
i_e = 8
t = 6
v0 = 10

# Second-order differential equation
v = v0 * exp((-x*t)/lmbda)
dvdt = -v0 * t * exp((-x*t)/lmbda) / lmbda
d2vdx2 = x**2 * t**2 * v0 * exp((-x*t)/lmbda) / (lmbda**2)


"""
Solutions to the linear cable equaiton are functions of both positon and time.
If the current injected is constant, the membrane potential settles to a steady-
state solution independent of time. This is an ordinary differential equation
we can solve.
"""

def dv2dx2(v, r, i, lmbda, solve=True):
    """
    0 is everywehre except in small region of size delta x around the injection site.
    General solution to this equaiton is v(x) = B1*exp(-x/lmbda) + B2*exp*x/lmbda) with
    coefficients undetermined B1 and B2.
    
    lmbda (lambda) is the length over which the localized electrode currents decay.
    v (velocity) is the speed at which ions travel across the membrane.
    r (resistance) is the membrane resistance and
    i (current) is the localized current we are studying.
    """
    if solve: # solve the equation
        return (np.exp(-x/lmbda), np.exp(x/lmbda)) # the terms in front of the coefficients B1 and B2
    return  (v-r*i)/lmbda**2 # just return the d2vdx2 term


"""
Rall model is a highly simplified model that captures important elements that affect the responses
of real neurons. Most neurons receive their synaptic inputs over complex dendritic trees.
In the Rall model, a compartment is a compact soma region that connects to a single equivalent
cylindrical cable replacing the entire dendritic region of the neuron. We can create equations that
relate the radii of the various branches with other parts of the tree. We use the geometry of the
various lengths and radii to construct relationships among the resistances of them.
"""

Rlambda = 4 # resistance of the cable of length constant lambda
L = 1.5*lmbda # for some cable length
x = 5 # at some distance along the cable

R1 = (Rlambda *(np.cosh(L/lmbda) - np.cosh(L-x)/lmbda)) / np.sinh(L/lmbda)
R2 = (Rlambda * np.cosh((L-x)/lmbda)) / np.sinh(L/lmbda)

# Because Rsoma acts in series with the two cable branches
Rsoma = R1 + R2
