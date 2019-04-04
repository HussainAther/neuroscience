import math
import numpy as np
import random
import decimal
import scipy.linalg
import numpy.random as nrand
import matplotlib.pyplot as plt

"""
It's like walking when you're drunk. Don't try this at home.
"""

class Brownian:
    """
    For storing objects and functions related to Brownian (brownian) motion.
    """
    def __init__(self, param):
        self.param = param

    def param(self, sigma, delta, time):
        param.all_sigma = sigma # uncertainty assocaited with each measurement according to the Gaussian distributino
        param.all_delta = delta # difference in movement at each step
        param.all_time = time # time over which we simulate Brownian motion

def brownian(x0, n, dt, delta):
    """
    Use Wiener prcoess from mathematics to simulate Brownian motion.
    """
    x0 = np.asarray(x0) # initialize array using numpy
    r = norm.rvs(size=x0.shape + (n,), scale=delta*sqrt(dt))# normalize using Gaussian
    out = np.empty(r.shape) # output result
    np.cumsum(r, axis=-1, out=out) # use cumulative sum as an approximation of Brownian motion
    out += np.expand_dims(x0, axis=-1) # initial condition
    return out

def brownian_motion_log_returns(param):
    """
    This method returns a Wiener process. The Wiener process is also called Brownian motion. For more information
    about the Wiener process check out the Wikipedia page: http://en.wikipedia.org/wiki/Wiener_process
    Param is hte model parameters object. Return the Brownian motion log values. 
    """
    sqrt_delta_sigma = math.sqrt(param.all_delta) * param.all_sigma
    return nrand.normal(loc=0, scale=sqrt_delta_sigma, size=param.all_time)


def brownian_motion_levels(param):
    """
    Returns a price sequence whose returns evolve according to a brownian motion
    """
    return convert_to_prices(param, brownian_motion_log_returns(param))

"""
In dynamic equilibrium, speed of the particles must be equal to mu*m*g. We can
calculate the density of Brownian particles rho at point x at time t to
satisfy the diffusion equation:

∂ρ/∂t = D * ∂^2ρ/∂x^2

in which D is the mass diffusivity defined as a function of Δ, the probability
of displacement.

D = integral from -inf to +inf of ((Δ^2)/2)(2*τ)•φ(Δ) dΔ

We solve this equation for N particles starting at the origin (t=0)

ρ(x, t) = (N/sqrt(4πDt))exp(-x^2 / 4Dt)
"""
D = 0.0016 # Diffusion coefficient for water in mm^2/s
N = 100
t = 1
x = 1

# in one dimension
rho = (N/np.sqrt(4*np.pi*D*t)) * np.exp((-x**2)/(4*D*t))

"""
As stated above, we use the Wiener process to describe a continuous-time stochastic process,
a type of Lévy process (explained below) that we find in describing the activity of a population of neurons.
They're used to study neural coding and can form generative models of brain imaging data.
Used in decision making in the brain can be fitted to behavioral data and used as regressors
in computational fMRI.
"""

delta = 0.25 # Wiener process parameter.
T = 30.0 # time
N = 500 # number of steps, related to time
dt = T/N # from T and N, deduce step size for integration
x = np.empty((2,N+1)) # initialize some values for x.
x[:, 0] = 0.0

brownian(x[:,0], N, dt, delta, out=x[:,1:])

"""
Lévy characterisation involves using a theorem that gives a necessry and sufficient
condition for a continuous R^n-valued stochastic process X to actually be n-dimensional
Brownian motion.
Let X = (X1...Xn) be a continuous stoachstic process on a probability space taking values in R^n.
The following two statements are equivalent:
(1.) X is a Brownian motion w. r. t. P, i.e., the law of X w.r.t. P is the same as the law of an n-dimensional
Brownian motion, i.e., the push-forward measure X*(P) is a classical Wiener measure on C_o([0, +inf]; R^n)
(2.) both (a.) X is a martingale w.r.t. P (and its own natural filtration) and (b.) for all 1 <= i, j <= n,
X_i(t)X_j(t) - δ_ijt is a martingale w.r.t. P (and its own natural filtration) where δ_ij dentoes the Kronecker delta.
"""
