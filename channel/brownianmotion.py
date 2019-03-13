import math
import numpy
import random
import decimal
import scipy.linalg
import numpy.random as nrand
import matplotlib.pyplot as plt

"""
It's like walking when you're drunk. Don't try this at home.
"""

def brownian_motion_log_returns(param):
    """
    This method returns a Wiener process. The Wiener process is also called Brownian motion. For more information
    about the Wiener process check out the Wikipedia page: http://en.wikipedia.org/wiki/Wiener_process
    :param param: the model parameters object
    :return: brownian motion log returns
    """
    sqrt_delta_sigma = math.sqrt(param.all_delta) * param.all_sigma
    return nrand.normal(loc=0, scale=sqrt_delta_sigma, size=param.all_time)


def brownian_motion_levels(param):
    """
    Returns a price sequence whose returns evolve according to a brownian motion
    :param param: model parameters object
    :return: returns a price sequence which follows a brownian motion
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
"""

