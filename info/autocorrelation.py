import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""
We use the autocorrelation of a time series function to look for repeating patterns
or serial correlation (correlation between signal at a given time and a later time).
"""

def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size // 2:]
