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

"""
Building a pipeline and tutorial for task fMRI analysis 
in nistats and functional connectivity analysis in nilearn.
"""

main_folder = "/home/finc/Dropbox/GitHub/nilearn_task_networks/"
top_dir = "/media/finc/Elements1/extracted_HCPWM/"

fmri_list = glob.glob(top_dir + "1*/**/tfMRI_WM_RL.nii.gz", recursive=True) # if you want all tfMRI_WM_RL.nii.gz files hanging from topdir
confounds_list = glob.glob(top_dir + "1*/**/Movement_Regressors.txt", recursive=True)

fmri_list = sorted(fmri_list)
confounds_list = sorted(confounds_list)

print(f"Number of subjects: {len(fmri_list)}")
