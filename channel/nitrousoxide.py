import numpy as np

from scipy.optimize import curve_fit

"""
Nitrous oxide (NO) diffusion from an irregular 3D structure. One such equation, known 
as the Crank-Nicholson scheme, is recommended for diffusive problems in one space dimension. 
"""

def diffusion(d2udx2, d2udy2):
    """
    With diffusion coefficient D, measure diffusion over time with using
    a two-dimensional diffusion equation.    
    """
    D = 1 
    dudt = D*(d2udx2+d2udy2
    return dudt

d2udx2 = np.linspace(-10, 10, 101)
d2udy2 = np.linspace(-10, 0, 51)
dudt = diffusion(d2udx2, d2udy2) 

best_vals, covar = curve_fit(diffusion, d2udx2, d2udy2)
