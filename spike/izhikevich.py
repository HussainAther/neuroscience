import numpy as np

"""
Izhikevich model for firing
"""

class IzhNeuron:
    """
    Initialize a class for variables of the neuron
    for differential equation use.
    """
    def __init__(self, label, a, b, c, d, v0, u0=None):
        self.label = label
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.v = v0
        self.u = u0 if u0 is not None else b*v0

class IzhSim:
    """
    Simulate the neuron.
    """
    def __init__(self, n, T, dt=0.25):
         self.neuron = n
         self.dt = dt
         self.t = t = arange(0, T+dt, dt)
         self.stim = zeros(len(t))
         self.x = 5
         self.y = 140
         self.du = lambda a, b, v, u: a*(b*v - u)

    def integrate(self, n=None):
        if n is None: n = self.neuron
        trace = zeros((2,len(self.t)))
        for i, j in enumerate(self.stim):
            n.v += self.dt * (0.04*n.v**2 + self.x*n.v + self.y - n.u + self.stim[i])
            n.u += self.dt * self.du(n.a,n.b,n.v,n.u)
            if n.v > 30:
                trace[0,i] = 30
                n.v = n.c
                n.u += n.d
            else:
                trace[0,i] = n.v
                trace[1,i] = n.u
        return trace
