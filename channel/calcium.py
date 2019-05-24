from PyDSTool import *
from Ch10_3_model import *

"""
Calcium parabolic bursting
"""

gentype = "radau" # Radau C-based integrator

# parameter arguments for PDSTool of making a Hodgkin-Huxley neuron
par_args = {"gr": 26, "gc": 3.25, "gx": 1.93,
            "tau_r": 5.6, "tau_c": 125, "tau_x": 300, "tau_v": 0.97,
                "Ena": 0.48, "Er": -0.95, "Ec": -0.95, "Ex": 1.4,
                "Iapp": 1.5, "As": 0, "f": 50, "An": 0.0}

ic_args = {"v": -0.82, "r": 0.4, "x": 0.04, "c": 0}

HH = makeHHneuron("HH_Ca", par_args, ic_args, const_I=True,
                  gentype=gentype)

HH.set(tdata=[0,1500])

def test_I(I):
    """
    Test our function to simulate calcium parabolic bursting.
    """
    HH.set(pars={"Iapp": I})
    traj = HH.compute("test")
    pts = traj.sample()
    plt.figure(1)
    plt.clf()
    plt.plot(pts["t"], pts["v"], "k")
    plt.plot(pts["t"], pts["r"], "r")
    plt.plot(pts["t"], pts["x"], "g")
    plt.plot(pts["t"], pts["c"], "y")
    plt.xlim(-5,max(pts["t"]))
    return traj, pts

traj, pts = test_I(1.6)

plt.show()

"""
Calcium-driven adaptation in neurons.
"""

def response(I, r):
    """
    Return response for spike generation rate r and external inupts I.
    """
