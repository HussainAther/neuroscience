import matplotlib.pyplot as plt
import numpy as np

from scipy.integrate import odeint

np.random.seed(1000)

# Start and end time (in milliseconds)
tmin = 0.0
tmax = 50.0

gK = 36.0 # Average potassium channel conductance per unit area (mS/cm^2)
gNa = 120.0 # Average sodoum channel conductance per unit area (mS/cm^2)
gL = 0.3 # Average leak channel conductance per unit area (mS/cm^2)
Cm = 1.0 # Membrane capacitance per unit area (uF/cm^2)
VK = -12.0 # Potassium potential (mV)
VNa = 115.0 # Sodium potential (mV)
Vl = 10.613 # Leak potential (mV)

# Time values
T = np.linspace(tmin, tmax, 10000)

# Potassium ion-channel rate functions

def alpha_n(Vm):
    '''
    Alpha rate constant for membrane voltage.
    ''' 
    return (0.01 * (10.0 - Vm)) / (np.exp(1.0 - (0.1 * Vm)) - 1.0)

def beta_n(Vm):
    '''
    Beta rate constant for membrane voltage.
    ''' 
    return 0.125 * np.exp(-Vm / 80.0)

# Sodium ion-channel rate functions
# Alpha and beta are given for the m-th and h-th ion channels.
# n, m, and h are dimensionless quantitites between 0 and 1 for 
# potassium channel acitvation, sodium channel activation,
# and sodium channel inactivation, respectively.

def alpha_m(Vm):
    return (0.1 * (25.0 - Vm)) / (np.exp(2.5 - (0.1 * Vm)) - 1.0)

def beta_m(Vm):
    return 4.0 * np.exp(-Vm / 18.0)

def alpha_h(Vm):
    return 0.07 * np.exp(-Vm / 20.0)

def beta_h(Vm):
    return 1.0 / (np.exp(3.0 - (0.1 * Vm)) + 1.0)

# We can calcualte the individual values n, m, and h.

def n_inf(Vm=0.0):
    return alpha_n(Vm) / (alpha_n(Vm) + beta_n(Vm))

def m_inf(Vm=0.0):
    return alpha_m(Vm) / (alpha_m(Vm) + beta_m(Vm))

def h_inf(Vm=0.0):
    return alpha_h(Vm) / (alpha_h(Vm) + beta_h(Vm))

def Id(t):
    '''
    Input stimulus current.
    '''
    if 0.0 < t < 1.0:
        return 150.0
    elif 10.0 < t < 11.0:
        return 50.0
    return 0.0

def compute_derivatives(y, t0):
    '''
    Use the definition of a derivative to compute the derivates at each point.
    '''
    dy = np.zeros((4,))
    
    Vm = y[0]
    n = y[1]
    m = y[2]
    h = y[3]
    
    # dVm/dt
    GK = (gK / Cm) * np.power(n, 4.0)
    GNa = (gNa / Cm) * np.power(m, 3.0) * h
    GL = gL / Cm
    
    dy[0] = (Id(t0) / Cm) - (GK * (Vm - VK)) - (GNa * (Vm - VNa)) - (GL * (Vm - Vl))
    
    # dn/dt
    dy[1] = (alpha_n(Vm) * (1.0 - n)) - (beta_n(Vm) * n)
    
    # dm/dt
    dy[2] = (alpha_m(Vm) * (1.0 - m)) - (beta_m(Vm) * m)
    
    # dh/dt
    dy[3] = (alpha_h(Vm) * (1.0 - h)) - (beta_h(Vm) * h)
    
    return dy

# State (Vm, n, m, h)
Y = np.array([0.0, n_inf(), m_inf(), h_inf()])

# Solve ODE system using scipy
# Vy = (Vm[t0:tmax], n[t0:tmax], m[t0:tmax], h[t0:tmax])
Vy = odeint(compute_derivatives, Y, T)
