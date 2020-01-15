import glob
import matplotlib.pyplot as plt
import nilearn
import numpy as np
import os
import pandas as pd

from nilearn import datasets, plotting, input_data, signal # datasets: for fetching atlas
from nilearn.input_data import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
from nistats.reporting import plot_design_matrix
from nistats.design_matrix import make_design_matrix

