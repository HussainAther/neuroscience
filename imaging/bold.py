import os
import matplotlib.pyplot as plt
import nitime
import nitime.timeseries as ts
import nitime.analysis as nta
import nitime.viz as viz

from matplotlib.mlab import csv2rec

"""
Extracting the average time-series from one signal, time-locked to the occurence of some type of event in 
another signal is a very typical operation in the analysis of time-series from neuroscience experiments. 
Therefore, we have an additional example of this kind of analysis in Auditory processing in grasshoppers

The following example is taken from an fMRI experiment in which a subject was viewing a motion stimulus, 
while fMRI BOLD was recorded. The time-series in this data set were extracted from motion-sensitive voxels 
near area MT (a region containing motion-sensitive cells) in this subject’s brain. 6 different kinds of trials 
could occur in this experiment (designating different directions and locations of motion). The following example 
shows the extraction of the time-dependent responses of the voxels in this region to the different stimuli.
"""

TR = 2. # interval
len_et = 15 # number of samples 

# Load data
data_path = os.path.join(nitime.__path__[0], 'data')
data = csv2rec(os.path.join(data_path, 'event_related_fmri.csv'))
