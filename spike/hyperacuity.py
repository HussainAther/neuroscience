import numpy as np
import scipy as sp

"""
In some cases of studying the reliability of neural computation,
we find ourselves at the extremes of what we can compute and determine
about the neural functions themselves. What is the scale of spatial precision?
We know the eye samples the world with a discrete lattice of photoreceptors.
In the human fovea this lattice spacing on the retina corresponds to an angular spacing of ~0.01 
degrees in the visual world.

In an acuity task, we distinguish one point source from two point sources of the 
same total intensity. We look at the intensity pattern on the retina from a single
point source.
"""

def I0(x):
    """
    This function represents a single intensity pattern we observe on the retina.
    """
    return np.sin(x)

# for some distance l that separates the point source from the retina
# and x as the distance between the points on the retina itself.
x = 5
l = 20
I = (1/2) * (I0(x-l/2) + I0(x+ l/2) # we add each half of the intensity into a total intensity.

"""
For small displacements l, the differnece between the patterns from one point sources
and the pattern from two points sources can be found by expanding I2 in a Taylor series
and keeping only the leading term.
"""

def d2I0dx(x):
    """
    Return the second derivative of I0 with respect to x.
    """
    return -np.sin(x)

deltaI = (1/8) * (l**2) * d2I0dx # pattern of light on retina

"""
If photon shot noise and other noise in the photoreceptors can be summarized by an 
effective spatial white noise added to the image, thne the signal to noise ration
for discrimination between two signals is an integral of the squared intensity
differences over the entire image.
"""

N0 = 10 # effective noise level
a = 100 # image size
integrand = (d**2) * x * (deltaI)**2
snr = (1/N0) * sp.integrate.quad(integrand, 0, a)
