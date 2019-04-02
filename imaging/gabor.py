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

def impulse(x, Nh):
    """
    Input signal has a real and imaginary component representing orthogonal directions. We can
    use a complex number that results fomr the convolution of the Fourier transform of the harmonic
    (sinusoidal) function. For some x input over range Nh. 
    """
    c = y*np.exp(-1j*2*n*np.pi*time/period) # convolution
    f = np.array([2*c(i)*np.exp(1j*2*i*np.pi*x/period) for i in range(1,Nh+1)]) # Fourier transform
    return f.sum

def gabor(sigma, theta, Lambda, psi, gamma):
    """
    Gabor filters with different frequencies to return features from an image in two dimensions.
    """
     
    Gcos = 
