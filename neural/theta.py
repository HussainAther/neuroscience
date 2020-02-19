import matplotlib.pyplot as plt
import numpy as np

"""
The theta model, or Ermentrout–Kopell canonical model, is a biological 
neuron model originally developed to model neurons in the animal Aplysia, 
and later used in various fields of computational neuroscience. The model 
is particularly well suited to describe neuron bursting, which are rapid 
oscillations in the membrane potential of a neuron interrupted by periods of relatively little oscillation.
"""

def I(t):
    """
    Adjust current input over time t.
    """
    return 0 

def dthetadt(theta, t):
    """
    dtheta/dt = 1 - cos(theta) + (1 + cos(theta)) I(t)
    """
    return 1 - np.cos(theta) + (1 + np.cos(theta)) * I(t)

"""
Ermentrout and Kopell shows the "oscillator death" with the model:

dtheta1/dt = omega_1 + sin(theta1)cos(theta2)

dtheta2/dt = omega_2 + sin(theta2)cos(theta1)
"""

def stabplot(theta1, theta2):
    """
    Plot the stability diagram in omega1, omega2 parameter space.
    """
    fig = plt.figure(figsize=(14, 18))
    gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1, 1, 2])
    ax1 = fig.add_subplot(gs[0, 1]) 
    omega1, omega2 = np.mgrid[100, 100]
    dtheta1dt = omega1 + np.sin(theta1)*np.cos(theta2)
    dtheta2dt = omega2 + np.sin(theta2)*np.cos(theta1)
    strm = ax1.streamplot(omega1, omega2, dtheta1dt, dtheta2dt, color=U, linewidth=np.array(5*np.random.random_sample((100, 100))**2 + 1), cmap="winter", density=10, minlength=0.001, maxlength = 0.07, arrowstyle="fancy", integration_direction="forward", start_points = seed_points.T)
    fig.colorbar(strm.lines)
    plt.tight_layout()
    plt.show()
