from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

"""
In neuroscience, we use kinetmatic  transformations to determine sensory-motor questions
of how the nervous system can choose a single configuartion of body parts (such as arms and legs) in space.
We can use mathematical methods of transformation and processing input values themselves to model
kinematics transformations. We can study how the central nervous system if the predictable
patterns of variability in limb movement can be attributed to these mathematical coordinate transformations.

We define the states of a neuron using ordinary differential equations.
"""

def minjerk(H1x,H1y,H2x,H2y,t,n):
    """
    Given hand initial position H1x,H1y, final position H2x,H2y and movement duration t,
    and the total number of desired sampled points n,
    Calculates the hand path H over time T that satisfies minimum-jerk.
    
    Flash, Tamar, and Neville Hogan. "The coordination of arm
        movements: an experimentally confirmed mathematical model." The
        journal of Neuroscience 5, no. 7 (1985): 1688-1703.
    
    """
    T = linspace(0,t,n)
    Hx = np.zeros(n)
    Hy = np.zeros(n)
    for i in range(n):
        tau = T[i]/t
        Hx[i] = H1x + ((H1x-H2x)*(15*(tau**4) - (6*tau**5) - (10*tau**3)))
        Hy[i] = H1y + ((H1y-H2y)*(15*(tau**4) - (6*tau**5) - (10*tau**3)))
    return T,Hx,Hy

"""
Sample 9 locations equally spaced in shoulder-joint elbow joint space. For each location,
use gaussian random noise to determine how the joint would move.
"""

l1 = 0.34
l2 = 0.46
angs = np.array([30.0,60.0,90.0]) * np.pi/180
plt.figure(figsize=(5,10))
for i in range(3):
  for j in range(3):
    a1 = angs[i]
    a2 = angs[j]
    plt.subplot(2,1,1)
    plt.plot(a1*180/pi,a2*180/pi,"r+")
    ex,ey,hx,hy = joints_to_hand(a1,a2,l1,l2)
    plt.subplot(2,1,2)
    plt.plot(hx,hy,"r+")
    for k in range(20):
      a1n = a1 + randn()*(sqrt(3)*pi/180)
      a2n = a2 + randn()*(sqrt(3)*pi/180)
      plt.subplot(2,1,1)
      plt.plot(a1n*180/pi,a2n*180/pi,"b.")
      ex,ey,hx,hy = joints_to_hand(a1n,a2n,l1,l2)
      plt.subplot(2,1,2)
      plt.plot(hx,hy,"b.")

plt.subplot(2,1,1)
plt.axis("equal")
plt.xlabel("SHOULDER ANGLE (deg)")
plt.ylabel("ELBOW ANGLE (deg)")
plt.subplot(2,1,2)
plt.axis("equal")
plt.xlabel("HAND POSITION X (m)")
plt.ylabel("HAND POSITION Y (m)")
