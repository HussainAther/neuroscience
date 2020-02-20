import matplotlib.pyplot as plt
import numpy as np

"""
Passive membrane subject to current pulse solved with
Backward Euler method.
"""

def bepswI(dt, Tfin):
    """
    Bakcwards Euler passive membrane 
    """
    VCl = -68 # mV source voltage
    A = 4*np.pi*1e-6 # cm^2 patch area
    Cm = 1 # micro F/cm^2 capacitance per area
    gCl = .3 # mS/cm^2 conductance
    tau = Cm/gCl # ms time constant	
