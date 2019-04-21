import numpy as np
import math
import scipy.stats as stats

from scipy.misc import factorial

"""
Estimate static stimulus values on the basis of spike-counting firing rates.
We estimate such a stimulus from the sequence of firing times of the spikes that the stimulus evokes.

If r is the spike-count firing rate of neuron a, we can obtain an estimate of the wind direction
on any given trial from the direction of the vector v_pop.

v_pop = the sum of (r/r_max)*c_a from a=1 to 4.

If we assume each neuron in a popultion fires independently, the firing-rate probability
for the population is the product of the individual probabilities:
"""

def popProb(problist, T, r):
    """
    Population probability of firing. problist is the lsit of probabilities
    for each nueron. T is the time period. r is the firing rate.
    """
    result = 1
    for prob in problist:
        factor = exp(-prob*T)
        factor *= (prob*T)**(r*T)
        factor / factorial(r*T)
        result *= factor
    return result

"""
The Cramér-Rao bound sets a limit of the variance of any estimate s_est according to:

(sigma_est)^2 >= (1 + b'_est(s))^2 / I_F(s)

in which b'_est(s) is the derivative of b_est(s), the bias on all stimulus variables. I_F(s) is the
Fisher information, a measure of encoding accuracy. Through this, the Fisher infomration limits the
accuracy with which any decoding scheme can extract an estimate of an encoded quantity.
"""

def fisherInformation(xvals, sigmavals):
    """
    Fisher information. xvals is the list of x positions at which the 
    spikes are measured, sigmavals is the list of the x position uncertainties.
    """
    npar = len(xvals)
    F = np.zeros([npar,npar])
    for x,sigma in zip(xvals,sigmavals):
        for i in range(npar):
            if i==0:
                dfdpi = x
            else:
                dfdpi = 1
        for j in range(npar):
            if j==0:
                dfdpj = x
            else:
                dfdpj = 1
            F[i,j] += sigma**-2*dfdpi*dfdpj
    return np.mat(F).I

"""
We define the theoretical decoding performance (DPth) of a linear classifier
derived from Averbeck and Lee's "Effects of Noise Correlations on Information
Encoding and Decoding" (2006). We start with DPth = phi(d'/2) where phi(x) is the 
cumulative normal function and d' = sqrt(deltaf^T (sigma^(-1)*(deltaf))) as the
signal-to-noise ratio generalized for a population of neurons. deltaf is the vector
joining the means of the population responses in the two stimulus conditions and 
sigma is the stimulus-invariant noise covariance matrix of the neuronal population.
We re-write the equation by rotating the original N-dimensional neural response
space along the eigenvectors of the covariance matrix such that we get:

DPth = phi((1/2)* |deltaf| * sqrt(summation from i to N of cos^2(thetai/sigmai^2)) 

in which sigmai^2 is the ith eigenvalue of the covariance matrix. 
"""
def vangle(v1, v2):
    """
    Return the angle between two vectors.
    """
    return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def dp(data):
     """
     For 1D or 2D data array of vectors, calculate the decoding performance.
     It's a method of finding the amount of informaiton encoded in a neuronal
     population responses. 
     """
     data = np.matrix(data) # convert to numpy matrix
     if data.shape[0] == 1: # if the matrix is one-dimensional
          xvals = range(len(data.shape[0])) # set the xvalues for the covariance matrix as the range of integers over the length of data
     elif data.shape[0] == 2: # if the matrix is two-dimensional
          xvals = data[0] # set the xvalues as the first array
          data = data[1] # set the data equal to the response variables
     else:
          raise ValueError
     c = np.cov(xvals, data) # get the covariance matrix of the neuronal population
     xvalsmean = mean(xvals) # mean of the xvalues
     datamean = mean(data) # mean of the data
     deltaf_abs = abs(datamean - xvalsmean) # absolute value of deltaf, connecting the two means
     deltaf_ang = np.angle(datamean - xvalsmean) # angle given by deltaf 
     eval, evec = np.linalg.eig(c) # eigenvalues and eigenvectors of the covariance matrix
     angles = []
     for i in range(len(evec)): # for each eigenvector 
         angles.append(vangle((i[0], i[1]), deltaf_ang)) # angles between each eigenvector and the deltaf angle
     summed = 0
     for i in range(len(eval)): # for each eigenvalue  
           num = np.cos(angles[i])**2 # numerator: take the cosine of each angle and square it
           den = eval[i]**2 # denominator: square each eigenvalue 
           summed += num/den # get the summation by adding each one up
     x = .5 * deltaf_abs * np.sqrt(summed) 
     return stats.norm.cdf(x) # return the cumulative normal distribution function of f.

