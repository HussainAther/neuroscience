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