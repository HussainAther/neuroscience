import numpy as np
import matplotlib.pyplot as plt

def compute(stim, rho, steps):
    """
    Compute the spike-triggered average (STA) from a stimulus and spike-train.
    stim is the stimulus time-series, rho is the spike-train time series, and
    steps is the number of timesteps. 
    """
    sta = np.zeros(steps)
    times = rho[steps:].nonzero()[0] + steps
    spikes = len(steps)
    for spike in times:
        sta += stim[spike-steps:spike]
    sta /= spikes
    return sta
