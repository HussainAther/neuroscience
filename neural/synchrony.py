import numpy as np

from scipy.stats.contingency import margins

"""
Linear measures of neuronal signal synchrony estimate the synchrony between two or sometimes
more continuous time series of brain activitiy which yield low values for independent time
series and high values for correlated time series. 

We can look at both cross-correlation and coherence measures of this synchrony.
"""

def crosscorrelation(tau, x, n):
    """
    For some time lag tau derived from normalized signals x and y of length N and with zero mean, we
    can calculate the unit variance as linear cross correlation and use that as a measure of synchronization.
    The absolute value is symmetric in the x and y directions and reaches a maxixmum of 1 for complete lag 
    synchronization.
    """
    if tau < 0:
        tau = -tau # flip the sign so we're looking at time lag in the same direction every time we run this function
    summ = 0
    for i in range(N-tau):
        summ += x[i+tau] * y[n] # account for time lag tau
    return (1/(N-tau)) * summ

def crossspec(w):
    """
    Return cross-spectrum of some signal w by performing Discrete Fourier transform as the x direction and its inverse for the y.
    """
    return np.fft.fft(w)*np.fft.ifft(w) # multiply the discrete fourier transform of the signal by its complex conjugate.

def coherence(w):
    """
    We can quantify linear correlations in the frequency domain with the cross spectrum for some signal w.
    """
    cross = crossspec(w)
    """
    If we normalize the amplitdue of the power spectrum of both the x and y systems, the we can calculate the 
    cross-spectrum as the coherence function:
    """
    num = abs(crossspec(w))**2 # numerator
    den = abs(np.fft.fft(w)*np.fft.fft(w)) * abs(np.fft.ifft(w) * np.fft.ifft(w)) 
    return num / den

"""
If our information is nonlinear, we can use a few more complicated techniques to measure synchrony.
Between a nonlinear relationship of X and Y, we can use uncertainty of a probability distribution
based on Shannon entropy.
"""

def jointprob(a):
    """
    Given a N-D array a, return the joint probability distribution. For use in the uncertainty function.
    """
    N = a.shape[0]
    jointProbs, edges = np.histogramdd(a, bins=N)
    return JointProbs / jointProbs.sum()

def uncertainty(j):
    """
    Use the Shannon entropy of a marginal distribution and Shannon entropy of the joint distribution
    to calculate nonlinear dependence using our uncertainty of a probability distribution. The joint
    distribution j between X and Y must be known as a 2-D array.
    """
    x, y = margins(j) # x and y describe the joint probability margins (marginal distribution)
    Hx = 0 # Shannon entropy of the x marginal distribution
    for i in x:
        Hx += i*np.log(i)
    Hx = -Hx # flip the sign for entropy
    Hy = 0 # mutatis mutandis for y
    for i in y:
        Hy += i.np.log(i)
    Hy = -Hy
    Hxy = 0 # Shannon entropy of the joint distribution
    for i in x:
        for j in y:
            Hxy += (i*j)*np.log(i*j)
    Hxy = -Hxy
    return Hx + Hy - Hxy # This can be a measure of mutual information from a joint probability distribution

"""
We can use Granger causality using the same principle from transfer entropy to test whether the prediction of a 
signal that relies only on its own past (univariate model) can be improved by incorporating past
information from the other signal (bivariate). We simply compare the univariate model to bivariate.
"""

def lrcoef(x,y):
    """
    For variable arrays x and y, estimate the coefficients of linear regression.
    """
    n = np.size(x) # size
  
    m_x, m_y = np.mean(x), np.mean(y) # get the means

    SS_xy = np.sum(y*x) - n*m_y*m_x # cross derivation
    SS_xx = np.sum(x*x) - n*m_x*m_x
  
    b_1 = SS_xy / SS_xx # regression coefficients
    b_0 = m_y - b_1*m_x
  
    return(b_0, b_1)

def uni(x, y):
    """
    Univariate models uses the information from within the x and y data to predict.
    """ 
    xn = 0 # output univariate model for x
    yn = 0 # for y 
    (ux, uy) = (np.std(x), np.std(y)) # standard error for x and y
    for i in range(len(x)):
        ax, ay = lrcoef(range(i), x[:i]) # model parameters for x of the linear regression model
        xn += ax*x[len(x)-i] + ux # univariate formula
    for i in range(len(y)): # mutatis mutandis for y
        ax, ay = lrcoef(range(i), y[:i])  
        yn += ay*y[len(y)-i] + uy 
    return xn, yn

def biv(x, y):
    """
    Bivariate model to predict. Compare with past prediction of the other signal to improve (if there is an improvement).
    """
    xn1 = 0 # same as for uni function but we split our summation into two components. One for each lr coefficient
    yn1 = 0
    xn2 = 0
    yn2 = 0 
    uxy = np.std(x) - np.std(y) # compare the two methods of standard error 
    uyx = np.std(y) - np.std(x)  
    for i in range(len(x)):
        axy, bxy = lrcoef(x[:i], y[:i]) # model parameters for x and y of the linear regression model
        xn1 += axy*x[len(x)-i] # bivariate formula for y
        xn2 += bxy*y[len(x)-i] + uxy
    for i in range(len(y)):
        ayx, byx = lrcoef(y[:i], x[:i]) # model parameters for y and x of the linear regression model. we flip the direction
        yn1 += ayx*x[len(x)-i] # bivariate formula for y
        yn2 += byx*y[len(x)-i] + uyx
    return (xn1+xn2), (yn1+yn2)

def granger(x, y):
    """
    Use linear regression of stochastic processes. Popular in economics. Only recently examined for potential in neuroscience.  
    """
    return biv(x, y) - uni(x,y) # Return the differences between the bivariate models of x and y with the univariate models of x and y.
