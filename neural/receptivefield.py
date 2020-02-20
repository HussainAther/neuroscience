import numpy as np
import matplotlib.pyplot as plt

from IPython.disply import Image

'''
In 1952 Hubel and Wiesel made an interesting discovery. When they inserted thin electrodes into a cat's visual cortex and presented small spots of light on a screen in front of the cat, some of the cells responded by firing spikes only when the spot was at a particular positions.
The same cell would be inhibited by spots located in different areas of visual field or it would not respond at all.
The responses of such neurons were fully characterised by the inhibitory and excitatory areas in the visual field. They called such neurons simple cells.

They also found another type of cells that manifested more complex responses. These cells, called complex cells, would respond only to a visual object of a particular shape and orientation no matter where it was positioned in the visual field. If the image was slightly rotated
or its shape was changed, the response was supressed.
'''

def receptive_field(XX, YY, theta, width, height):
    '''
    Generate a receptive field.
    '''
    x_r = XX * np.cos(theta) + YY * np.sin(theta)
    y_r = -XX * np.sin(theta) + YY * np.cos(theta)
    gaussian = np.exp(-np.pi*(y_r**2/height**2)) * (np.abs(x_r) < width)
    F = 1./width * 0.75
    complex_grating = np.exp(1j*2*np.pi*F*x_r)
    return gaussian * complex_grating.real

'''
Receptive field is the set of the stimuli that evoke responses in single neurons. In case of visual receptive fields these are patches of light impaling on retina that are most effecitve in driving a neuron. It is customary to represent visual receptive field as images showing the preffered stimuli.

We will now generate a sample receptive field using numpy. We start with a periodic pattern of light and dark regions along one axis. To represent it as an image we need a two-dimensional grid of $x$ and $y$ coordinates — one pair per pixel — which we store in two numpy arrays XX and YY each containing
only one of the coordinates. This is best illustrated by the following schematics:
'''

xmin, xmax, dx = -20, 20, 0.04
ymin, ymax, dy = -20, 20, 0.04

YY, XX = np.mgrid[xmin:xmax:dx, ymin:ymax:dy]

spatial_freq = 0.1

# receptive field
rf = np.sin(2 * np.pi * spatial_freq * XX)
rf_height = 10
rf_width = 4
rf_angle = 0
rf_ampl = 1 #Hz
rf_angle = np.pi/3.

'''
Use rotated system of coordinates.
'''

W = rf_ampl * receptive_field(XX, YY, rf_angle, rf_width, rf_height)

rf = np.cos(2 * np.pi * spatial_freq * x_r)
plt.matshow(rf)
plt.gray()

'''
To make the receptive field spatially-limited we mask it with a Gaussian along the stripes and with a rectangle across the stripes.
'''

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

'''
Response of a simple cell to the receptive field.
'''

def bar(XX, YY, length, width, angle):
    x_r = XX * np.cos(angle) + YY * np.sin(angle)
    y_r = -XX * np.sin(angle) + YY * np.cos(angle)
    return  1 * (np.abs(x_r) < length / 2) * (np.abs(y_r) < width / 2)


bar_length = 10.
bar_width = 3.
bar_angle = np.pi/2.
bar_lumosity = 1
b = bar(XX, YY, bar_length, bar_width, bar_angle)
plt.imshow(b, extent=(xmin, xmax, ymin, ymax))
plt.gray()

def rectify(x):
    return x * (x > 0)

def simple_cell_model(receptive_field, baseline_rate, dxy):
    '''
    Return how a simple cell will respond to the receptive field
    '''
    dx, dy = dxy
    def _response(stim):
        r = np.sum(stim * receptive_field) * dx * dy + baseline_rate
        return rectify(r)
    return _response


baseline_rate = 20
simple_cell = simple_cell_model(W, 20, (dx, dy))

r = simple_cell(b)
print('Output firing rate: {:.1f} Hz'.format(r))
print('Change from baseline: {:.1f} Hz'.format(r - baseline_rate))

plt.figure(figsize=(3, 3))
plt.imshow(W, extent=(xmin, xmax, ymin, ymax))
b1 = bar(XX, YY, bar_length, bar_width, bar_angle)
b2 = bar(XX - bar_width, YY, bar_length, bar_width, bar_angle)
plt.contour(XX[0,:], YY[:,0], b1)
plt.contour(XX[0,:], YY[:,0], b2)
plt.xlabel('azimuth (deg)')
plt.ylabel('longlitude (deg)')

'''
The two bars are superimposed on the receptive field using white lines. In order to plot them we used contour plot from matplotlib library, which shows contour lines of constant levels in a two dimensional map.

Let us compare the responses of this simple neuron to presentations of each of the bar separately (Bar 1 and 2), sum of the responses, and the response to both bars presented together (Bar 1+2). Since we are only
interested in change from the baseline firing rate, we first subtract the baseline.
'''

r1 = simple_cell(b1) - baseline_rate
r2 = simple_cell(b2) - baseline_rate
r1and2 = simple_cell(b1 + b2) - baseline_rate
plt.bar(range(4), [r1, r2, r1+r2, r1and2])
plt.xticks(np.arange(4) + 0.5, ['Bar 1', 'Bar 2', 'sum of\n responses', 'Bar 1 + 2'])
plt.ylabel('change of firing rate (Hz)')
