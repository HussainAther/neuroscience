import brian2 as b2
import matplotlib.pyplot as plt
import numpy as np

from neurodynex.phase_plane_analysis import fitzhugh_nagumo

%{
FitzHugh-Nagumo (fitzhugh nagumo) model
}%

fitzhugh_nagumo.plot_flow()

fixed_point = fitzhugh_nagumo.get_fixed_point()
print("fixed_point: {}".format(fixed_point))

plt.figure()
trajectory = fitzhugh_nagumo.get_trajectory()
plt.plot(trajectory[0], trajectory[1])
