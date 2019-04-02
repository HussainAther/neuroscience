import numpy as np

"""
In image processing, the Gabor filter is a linear filter for texture analysis. 
It detects frequency content in the image in specific directions in a localized
region around the point or region of analysis.

We can investigate the extent to which it's similar to perception of the human visual system
with biologically plausible deep learning techniques. If we were to fix hidden layer weights
of a neural network with random Gabor filters or train them with unsupervised methods (such as 
principal component analysis), we can determine whether unsupervised learning leads to better
performanance than networks that are completely connected within themselves. 
"""

def impulse(x, y):
    """
    Input signal has a real and imaginary component representing orthogonal directions. We can
    use a complex number that results fomr the convolution of the Fourier transform of the harmonic
    (sinusoidal) function. For some two-dimensional input x and y in those directions. 
    """
    return np.convolve(x, y)

def gabor(sigma, theta, Lambda, psi, gamma):
    """
    Gabor filters with different frequencies to return features from an image in two dimensions.
    sigma, theta, Lambda, psi, and gamma are the resulting variables of the complex number 
    of the Fourier transform. In this transformation, Lambda is the sine factor wavelength, 
    theta is the normal orientation to the parallel stripes of a Gabor function,
    psi is the phase offset, sigma is the standard deviation of the Gaussian envelope, and
    gamma is the spatial aspect ratio of the ellipticity of the Gabor function support.
    """
    sigmax = sigma
    sigmay = float(sigma) / gamma 
    nstds = 3 # number of standard deviation sigma
