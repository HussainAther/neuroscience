import numpy as np

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
I = (1/2)*(I0(x-l/2) + I0(x+ l/2) # we add each half of the intensity into a total intensity.
