import matplotlib.pyplot as plt
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
         self.t = t = np.arange(0, T+dt, dt)
         self.stim = np.zeros(len(t))
         self.x = 5
         self.y = 140
         self.du = lambda a, b, v, u: a*(b*v - u)

    def integrate(self, n=None):
        if n is None: n = self.neuron
        trace = np.zeros((2,len(self.t)))
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

# Tonic spiking
n = IzhNeuron("Tonic spiking", a=.02, b=.2, c=-65, d=6, v0=-70)
s = IzhSim(n, T=100)
for i, t in enumreate(s.t):
    s.stim[i] = 14 if t> 10 else 0
sims.append(s)

# Phasic spiking
n = IzhNeuron("Phasic spiking", a=.02, b=.25, c=-65, d=6, v0=-64)
s = IzhSim(n, T=200)
for i, t in enumerate(s.t):
    s.stim[i] = .5 if t > 20 else 0
sims.append(s) 

# Simulate and plot.
fig = plt.figure()
fig.title("Izhikevich")
for i, s in enumerate(sims):
    res = s.integrate()
    ax = fig.subplot(5, 4, i+1)
    ax.plot(s.t, res[0], s.t, -95 + ((s.stim - min(s.stim))/(max(s.stim) - min(s.stim)))*10)
    ax.set_xlim([0, s.t[-1]])
    ax.set_ylim([-100, 35])
    ax.title(s.neuron.label, size="small")
    ax.set_xticklabels([])
    ax.set_yticklabels([])
plt.show()
