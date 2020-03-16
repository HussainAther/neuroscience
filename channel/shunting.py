"""
Shunting equation for membrane dynamics

The Shunting equation has the form of the membrane equation on 
which cellular neurophysiology is based. This membrane equation 
is the voltage equation that appears in the equations of Hodgkin 
and Huxley. In other words, the gedanken experiment shows how the 
noise-saturation dilemma is solved by using the membrane, or shunting, 
equation of neurophysiology to describe cells interacting in on-center 
off-surround anatomies. Because on-center off-surround anatomy and 
shunting dynamics work together to solve the noise-saturation dilemma, 
it is reasonable to predict that they coevolved during evolution.
"""

def shunting(A, B, I):
    """
    For I total input, B spontaneous decay rate of the excitation, 
    and B number of excitable states, calculate the equilibrium value 
    of the number of excited states.
    """
    return (B*I)/(A+I)

def voltagerate(V, g, Vt, C):
    """
    For Vt voltage rate dependent on time, tuples V and g that are
    three numbers in length that describe excitatory, inhibitory, and
    passive saturation voltages and conductances, calculate the 
    voltage rate of the cell with a given capacitance C.
    """ 
    return ((V[0]-Vt)*g[0] + (V[1]-Vt)*g[1] + (V[2]-Vt)*g[2])/C 
