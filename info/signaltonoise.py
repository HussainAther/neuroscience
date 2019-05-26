import numpy as np
import matplotlib.pyplot as plt
import nitime.utils as utils
import nitime.timeseries as ts
import nitime.viz as viz

"""
Calculate signal-to-noise (signal to noise snr) ratio. Where SNR(omega) is the ratio of the signal power and the noise power at the 
frequency band centered on \omega.This equation holde true for a Gaussian channel and is an upper bound for all other cases.
The signal power is estimated as the power of the mean response to repeated presentations of the same signal and the noise power 
is calculated as the average of the power in the deviation from this average in each trial
"""
