import numpy as np

"""
Dynamic causal modeling for inferences about functional integration
for fMRI time series analysis. Bilinear model of the neurodynamics and 
extended balloon model for hemodynamics. Describe using a multivariate
differential equation.
"""

# zt is the time derivative of neuronal activity with fixed connectivity
# A, modulation of connectivity u, external inputs B and C. Matrix A is the
# fixed or average coupling among nodes without the exogeneous input u, matrix B
# are the change in latent coupling for each input.
zt = (A + np.dot(u, B) * z + C*u

"""
BOLD response y is the BOLD signal convolution of inputs
"""
