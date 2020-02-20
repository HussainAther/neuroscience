import matplotlib.pyplot as plt
import numpy as np

'''
Leaky integrate and fire (leaky-integrate-and-fire LIF) 
'''

# Initialize vairables.
T = 50 # total time to simulate (msec)
dt = 0.125 # simulation time step (msec)
time = arange(0, T+dt, dt) # time array
trest = 0 # initial refractory time
Vm = zeros(len(time)) # potential (V) trace over time
Rm = 1 # resistance (kOhm)
Cm = 10 # capacitance (uF)
taum = Rm*Cm # time constant (msec)
tauref = 4 # refractory period (msec)
Vth = 1 # spike threshold (V)
Vspike = 0.5 # spike delta (V)

# Input stimulus
I = 1.5 # input current (A)

# Iterate over each time step.
for i, t in enumerate(time):
   if t > trest:
       Vm[i] = Vm[i-1] + (-Vm[i-1] + I*Rm) / taum * dt
   if Vm[i] >= Vth:
       Vm[i] += Vspike
trest = t + tauref

# Plot.
plt.plot(time, Vm)
plt.title('Leaky Integrate-and-Fire Example')
plt.ylabel('Membrane Potential (V)')
plt.xlabel('Time (msec)')
plt.ylim([0,2])
plt.show()
