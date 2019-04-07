import numpy as np
import scipy as sp

"""
We use Markov renewal processes to generalize Markov jump processes. Renewal processes are stochastic
point processes to describe a sequence of events in time. In the narrow sense, they presuppose stationary 
input and are defined by the fact that the state of the system, and the probability of generating the
next event depends only on the "age" (t-t_hat) of the system. 
"""

def P0(s):
    """
    Calculate the probability that the next event occurs at time t + s given that the last event time  
    was time t. We use the Laplace transform as the Fourier transform of gamma(x)f(x)exp(b) in which
    gamma is a step function that gets rid of the negative part of the integral and exp(bx) is the
    real part of the complex exponential. In this case, we take the complex number form of the interval
    s and integrate over that to return the renewal probability. 
    """
    com = np.complex(1, s) # complex number form of the imaginary component s
    a = 2
    dt = .5 # some time differential
    f = lambda u: 1/(u-a)
    return (f(com)).real*np.cos(com.imag*dt) * 2*np.exp(com.real*dt)/np.pi 

def survivor(t, t_hat):
    """
    The survivor function defines the probability that the neuron stays quiescent between t_hat and t.
    It's the probability a neuron "survives" from t_hat to t without firing.
    """
    x0 = lambda x: P0(x) 
    return 1 - sp.integrate.quad(x0, t, t_hat)

def delta(n):
    """
    Dirac delta function.
    """
    if n == 0:
        return 1
    return 0

def normauto 
