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
lambda = 10 # Square root of (r_m/r_i)
V = [m.Var(value = np.sin(2*xpos[i])) for i in range(npx)]

# Discretization of space
xpos = np.linspace(0,10,1)
dx = xpos[1]−xpos[0]

# Equation
m.Equation(V == V_0 * exp((-x)/lambda))

# Solve it!
m.solve()

"""
When x = lambda, the ration of V to V_0 is exp(-1) or about .37. Lambda is the characteristc
length constant of the cable at this point. The decay of membrane potential along an infinte dendritic
cable is described by the length constant.
"""