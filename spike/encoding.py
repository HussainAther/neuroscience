import numpy as np
import matplotlib.pyplot as plt

"""
Encoding models address the question of how the firing rate of a neuron is related to the stimulus that is presented.



A biological example would be where we stimulate a neuron with a particular light intensity.
If the firing rate depends linearly on this light intensity then the first property tells us
that doubling (or halving) the light intensity should double (or halve) the firing rate of the neuron.
The second property says that adding two light intensities together should produce a firing rate that
is equal to the sum of the firing rates to each individual lights. This assumption of linearity does
not hold really well in reality as neural firing rates cannot increase indefinitely.
"""

# Definition of the stimulus
dt = 0.01                # Time step [s]
T  = 2                   # Total duration of the simulation
t = np.arange(0, 2, dt)  # Time points

f = 0.5 # Frequency
A = 1   # Amplitude
stim = A*np.sin(2*np.pi*f*t) + np.random.normal(0.0, 0.05, len(t))

# Delayed amplification
response = np.zeros(len(stim))
delay = 50
amplification = 1.5

for i in range(delay,len(stim)):
    response[i] = amplification*stim[i-delay]

plt.plot(t,stim)
plt.plot(t,response)
plt.xlabel('Time (s)')
plt.legend(['Stimulus','response'])
