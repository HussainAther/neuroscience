import numpy as np

"""
In states of nonequilibrium the Nernst (nernst) potential (V_n) doesn"t need to equal the actual
potential jump across a membrane. The flux of ions through a membrane will be dissipative and
proportional to a net driving force.

Ohmic conductance hypothesis dictates:
j_qi = z_ie*j_i = (deltaV - V_n) * g_i

in which j_i is the number of ions of type i per area per time crossing the membrane,
j_qi is the electric charge flux, z_ie is the charge on one ion, and g_i is the constant
of proportionality of conductance per area.
"""

g_K = 1 # conductance per area normalize with respect to potassium
g_Na = 25 # sodium, based on  Hodgkin and B. Katz (1948)
e = -1 # charge of electron
V_n = -.09 # mV Nernst sodium

def V(i):
    """
    Some voltage of the membrane that depends upon time.
    """
    return np.sin(i)

def pump(i):
    """
    We define some current source in parallel with the other elemnets that must do
    work to pump sodium ions uphill (against their electrochemical potential gradient). 
    """
    return 4*i

def ionPump(ion):
    """
    Return the electric chrage flux of an ion over time as it pumps through a membrane.
    """
    result = []
    for i in range(0, 180):
        result.append((g_Na/e) * (V(i) - V_n) + pump(i))
    return (g_Na / e) * V(i)
