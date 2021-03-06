import numpy as np

"""
In image processing, the Gabor (gabor) filter is a linear filter for texture analysis. 
It detects frequency content in the image in specific directions in a localized
region around the point or region of analysis.

We can investigate the extent to which it"s similar to perception of the human visual system
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
    xmax = max(abs(nstds * sigma_x * np.cos(theta)), abs(nstds * sigma_y * np.sin(theta)))
    xmax = np.ceil(max(1, xmax)) # ceiling function returns the smallest integer i such that i >= x.
    ymax = max(abs(nstds * sigma_x * np.sin(theta)), abs(nstds * sigma_y * np.cos(theta)))
    ymax = np.ceil(max(1, ymax))
    xmin = -xmax # flip the sign to get the minimum in the x direction
    ymin = -ymax # same
    (y, x) = np.meshgrid(np.arange(ymin, ymax + 1), np.arange(xmin, xmax + 1))

    # Rotation in accordance with our theta variable to extract features in both directions. 
    x_theta = x * np.cos(theta) + y * np.sin(theta)
    y_theta = -x * np.sin(theta) + y * np.cos(theta)

    # Gabor function for each dimensino
    gb = np.exp(-.5 * (x_theta ** 2 / sigma_x ** 2 + y_theta ** 2 / sigma_y ** 2)) * np.cos(2 * np.pi / Lambda * x_theta + psi)
   
    return gb

"""
Gabor fatures can identify the script of a word in a multilingual document. Different frequencies and rotatiosn of Gabor functions
let us localize and extract text-only regions from complex document images. We can also use it in optical character recognition,
iris recognition, and fingerprint recognition.

Dennis Gabor also introduced Gaboe wavelets that use complex functions to serve as a basis for Fourier
transforms in information theory application. They"re non-orthongonal, so it"s difficult to create efficient
decomposition into the basis.
"""

def wavelet(x):
    """
    We use a 1-D Gaussian modulated by a complex exponential for our Gabor wavelet.
    x0 is the center of the x distribution, k0 is the rate of modulation, and alpha is the rate of the exponential
    drop-off as the Gabor function becoms exponentially suppressed. Use input x for the position of the wavelet in one
    dimension. 
    """
    x0 = 0
    k0 = 1
    alpha = 10
    exp1 = -(x-x0)**2 # first exponential term
    exp2 = -np.imag(k0*(x-x0) # second exponential term
    return np.exp(exp1) * np.exp(exp2)
