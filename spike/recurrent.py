import numpy as np

"""
Recurrent spiking networks
"""
# Initialize vairables.
N = 10 # matrix dimension
T = 50 # total time to simulate (msec)
dt = 0.125 # simulation time step (msec)
time = arange(0, T+dt, dt) # time array
trest = 0 # initial refractory time
Vm = zeros(len(time)) # potential (V) trace over time
Rm = 1 # resistance (kOhm)
Cm = 10 # capacitance (uF)
taum = Rm*Cm # time constant (msec)
tauref = 4 # refractory period (msec)
taupsc = 3 # presynaptic current time constant (msec)
Vth = 1 # spike threshold (V)
Vspike = 0.5 # spike delta (V)

# Synapse weight matrix equally weighted ring connectivity
synapses = np.eye(N)
synapses = np.roll(synapses, -1, 1)

# Synapse current model
def Isyn(t):
    """
    t is an array of times since each neuron's 
    last spike event.
    """
    t[np.nonzero(t < 0)] = 0
    return t*np.exp(-t/taupsc)

lastspike = np.zeros(N) - tauref

# Simulate the network.
raster = np.zeros([N,len(time)])*np.nan
for i, t in enumerate(time[1:],1):
    active = np.nonzero(t > lastspike + tauref)
    Vm[active,i] = Vm[active,i-1] + (-Vm[active,i-1] + I[active,i-1]) / taum * dt
    spiked = np.nonzero(Vm[:,i] > Vth)
    lastspike[spiked] = t
    raster[spiked,i] = spiked[0]+1
    I[:,i] = Iext + synapses.dot(Isyn(t - lastspike))
