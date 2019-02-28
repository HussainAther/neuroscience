import numpy as np

"""
In the leaky integrate-and-fire model of the neuron, we add a "leak" to the membrane potential
that reflects the diffusion of ions that occurs throguh the membrane when some equilibrium is not
reached in the cell. we model:

I(t) - V_m(t)/R_m = C_m * dV_m(t) / dt

in which I is the current, V_m is the membrane voltage, R_m is the membrane resistance, and C_m
is the membrane capacitance. If the cell's input current exceeds some threshold (I_th), it fires.
Otherwise it just leaks out any change in potential. The firing frequency is:

f(I) = (0 if I <= I_th) or ([t_ref - R_m*C_m*log(1 - V_th/(I*R_m)) if I > I_th)

We can integrate this over time to find the number of times the cell fired.
"""

R_m = 5

def V_m(t): # potential function
    result = []
    for i in range(t):
        result.append(i*3 + 2)
    return result



