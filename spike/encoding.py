import numpy as np
import matplotlib.pyplot as plt
from numpy import zeros

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
dt = 0.01 # Time step [s]
T  = 2  # Total duration of the simulation
t = np.arange(0, 2, dt) # Time points

f = 0.5 # Frequency
A = 1 # Amplitude
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

# Running average
N= 10
response = np.zeros(len(stim))

for i in range(N, len(stim)):
        response[i] += np.sum(stim[(i-N):i])/N

plt.plot(t, stim)
plt.plot(t, response)
plt.xlabel('Time (s)')
plt.legend(['Stimulus','response'])

# Leaky average
response = np.zeros(len(stim))
alpha = 30

la_filter = np.exp(-alpha*np.arange(0,N))
la_filter = la_filter/np.sum(la_filter)

for i in range(N, len(stim)):
    for k in range(N):
        response[i] += stim[i-k]*la_filter[k]

plt.plot(t, stim)
plt.plot(t, response)
plt.xlabel('Time (s)')
plt.legend(['Stimulus','response'])

# Retinal ganglion cells
[X,Y] = np.meshgrid(np.arange(-2,2,0.1),np.arange(-2,2,0.1))
B = 1.3
sigma_c = 0.3
sigma_s = 0.4

rgc = (1/(2*np.pi*sigma_c**2))*np.exp(- (X**2 + Y**2)/(2*sigma_c**2)) - (B/(2*np.pi*sigma_s**2))*np.exp(- (X**2 + Y**2)/(2*sigma_s**2))

plt.subplot(1,2,1)
plt.imshow(rgc)
plt.axhline(20,color='g')
plt.title('On center ganglion cell')
plt.axis('off')
plt.subplot(1,2,2)
plt.plot(rgc[20,:])
plt.axis('off')
plt.axhline(0,color='k')
plt.title('Receptive field cross section')

# Primal visual cortex single cells
[X,Y] = np.meshgrid(np.arange(-2,2,0.1),np.arange(-2,2,0.1))

f     = 0.6     # Spatial frequency
theta = 0.0     # Orientation
phi   = np.pi/2 # Phase offset
sigma = 0.3     # Standard deviation of the gaussian kernel
gamma = 0.6     # Aspect ratio

X_rot = X*np.cos(theta) + Y*np.sin(theta)
Y_rot = -X*np.sin(theta) + Y*np.cos(theta)

g = np.exp( - (X_rot**2 + (gamma**2)*(Y_rot**2))/(2*sigma**2))*np.cos(2*np.pi*f*X_rot + phi)
plt.subplot(1,2,1)
plt.title('Simple cell RF')
plt.axis('off')
plt.imshow(g)
plt.subplot(1,2,2)
plt.plot(g[20,:])
plt.title('RF cross-section')
plt.axis('off')
plt.axhline(0,color='k')
