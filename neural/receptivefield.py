import numpy as np
import matplotlib.pyplot as plt
from IPython.disply import Image

"""
In 1952 Hubel and Wiesel made an interesting discovery. When they inserted thin electrodes into cat's visual cortex and presented small spots of light on a screen in front of the cat, some of the cells responded by firing spikes only when the spot was at a particular positions.
The same cell would be inhibited by spots located in different areas of visual field or it would not respond at all.
The responses of such neurons were fully characterised by the inhibitory and excitatory areas in the visual field. They called such neurons simple cells.

They also found another type of cells that manifested more complex responses. These cells, called complex cells, would respond only to a visual object of a particular shape and orientation no matter where it was positioned in the visual field. If the image was slightly rotated
or its shape was changed, the response was supressed.
"""
def receptive_field(XX, YY, theta, width, height):
    """
    Generate a receptive field.
    """
    x_r = XX * np.cos(theta) + YY * np.sin(theta)
    y_r = -XX * np.sin(theta) + YY * np.cos(theta)
    gaussian = np.exp(-np.pi*(y_r**2/height**2)) * (np.abs(x_r) < width)
    F = 1./width * 0.75
    complex_grating = np.exp(1j*2*np.pi*F*x_r)
    return gaussian * complex_grating.real

"""
Receptive field is the set of the stimuli that evoke responses in single neurons. In case of visual receptive fields these are patches of light impaling on retina that are most effecitve in driving a neuron. It is customary to represent visual receptive field as images showing the preffered stimuli.

We will now generate a sample receptive field using numpy. We start with a periodic pattern of light and dark regions along one axis. To represent it as an image we need a two-dimensional grid of $x$ and $y$ coordinates — one pair per pixel — which we store in two numpy arrays XX and YY each containing
only one of the coordinates. This is best illustrated by the following schematics:
"""

xmin, xmax, dx = -20, 20, 0.04
ymin, ymax, dy = -20, 20, 0.04

YY, XX = np.mgrid[xmin:xmax:dx, ymin:ymax:dy]

spatial_freq = 0.1
rf = np.sin(2 * np.pi * spatial_freq * XX)

rf_height = 10
rf_width = 4
rf_angle = 0.
rf_ampl = 1. #Hz

rf_angle = np.pi/3.

"""
Use rotated system of coordinates.
"""
W = rf_ampl * receptive_field(XX, YY, rf_angle, rf_width, rf_height)

rf = np.cos(2 * np.pi * spatial_freq * x_r)
plt.matshow(rf)
plt.gray()

"""
To make the receptive field spatially-limited we mask it with a Gaussian along the stripes and with a rectangle across the stripes.
"""

rf_height = 6.
rf_width = 7.5

mask_x = (np.abs(x_r) < rf_width)
mask_y = np.exp(-(y_r**2/rf_height**2))

rf *= mask_x * mask_y

plt.matshow(rf)


plt.figure(figsize=(6,2.5))
plt.subplot(121)
plt.imshow(W, extent=(xmin, xmax, xmin, xmax), aspect='equal')
plt.xlabel('position (deg)')

plt.axhline(0, ls='--', color='k')
plt.gray()
cbar = plt.colorbar()
plt.subplot(122)
plt.plot(XX[500,:], W[500, :])
plt.xlabel('position (deg)')
cbar.set_label('receptive field')
plt.tight_layout()

"""
Response of a simple cell to the receptive field.
"""

def bar(XX, YY, length, width, angle):
    x_r = XX * np.cos(angle) + YY * np.sin(angle)
    y_r = -XX * np.sin(angle) + YY * np.cos(angle)
    return  1 * (np.abs(x_r) < length / 2) * (np.abs(y_r) < width / 2)
