import numpy as np
import scipy.integrate as integrate

"""
Since pyramidal cell are aligned in parallel, they form a dipole layer of thickness d
when they are synchronized within a cortical molecule. We identify these molecules with
teh anatomical columns in order to collective dipole field potential (as local field potential
generated by a mass of apprxoiamtely 10k pyramidal cells).

The differential of current is proportional to the infinitesimal area in cylindrical coordinates.
"""

def j(x):
    """
    Current density.
    """
    return 3 * x

dI = j * dA

"""
Where the current density j is assumed to be constant scalar within one column.
The differential of the potential dphi at a distance z perpendicular to a cortical column
with radius R cotributed by the current dI is
"""

sigma = 5.670e10-8 # Stefan-Boltzmann constant
d = 4 # for some distance

dphi = (1/4*np.pi*sigma) * j * d *( x - xp) / (x - xp)**3 * dA

"""
in which xp (x') varies across the area of the module. We can integrate this integrand over dA to
determine the potential. By taking advantage of symmetric of the cylindrical area, we can re-write
this as
"""

def phi(z, d, R):
    """
    Potential as dependent upon the distance z for a cylinder of radius R
    """
    constant = (2* np.pi * j(d)) / (4 * np.pi * sigma)
    integrated = integrate.quad(lambda x: x / (x**2 + z**2)**(3/2), 0, R)
    return constant * integrated

"""
Or, performing the integration, we get:
"""

R = 5 # some value idk
z = 4 # another value lol
d = 3

phi = (j(d) * (2*sigma)) * (1 - z/(R**2 + z**2)**(1/2))

"""
Electroencephalograms

We can extend these queations of the summed potential resulting fomr the synchronized
synaptci activity of all pyramidal neurons in a column in a distance z
above the cortical ggray matter. By integrating over a larger domain of cortical tissue,
we can get an estimator of the electrocorticogram (ECoG).

We need to consider the potential generated by a poiny source I at a distance -h from the interface
in the semi-space G1. Another source is I + I''. In the semi-space G2 for the EEG measurement,
we can solve it with:
"""
# G2 is an islate so we can set
sigma2 = 0
sigma1 = sigma

phi = (1/(2*np.pi*(sigma1 + sigma2))) * I/np.sqrt(r**2 + (z + h)**2)
