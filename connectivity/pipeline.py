import glob
import matplotlib.pyplot as plt
import nilearn
import numpy as np
import os
import pandas as pd

from nilearn import datasets, plotting, input_data, signal # datasets: for fetching atlas
from nilearn.input_data import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure, sym_matrix_to_vec, vec_to_sym_matrix
from nistats.reporting import plot_design_matrix
from nistats.design_matrix import make_design_matrix
from scipy import stats
from sklearn import neighbors

"""
Building a pipeline and tutorial for task fMRI analysis 
in nistats and functional connectivity analysis in nilearn.
"""

def calculate_lsn_edges(A, labels):
    """
    Function calculates number of edges between and within predefined large-scale networks (LSNs).
    The function takes binary symetrical adjacency matrix, module assignment of each ROI and calculate number of edges between and within each
    large-scale network.
    Parameters
    ------------
    array: N x N binary ajdacency matrix
    array: N-length vector with module assignment for each node
    Returns
    ------------
    array: M x M matrix with number of edges between each module
    """
    columns = np.unique(labels)
    lsn_matrix = np.zeros((len(labels), len(columns)))
    lsn_edges = np.zeros((len(columns), len(columns)))

    for col in range(len(columns)):
        module = columns[col, ]
        for row in range(len(labels)):
            if (labels[row, ] == module):
                lsn_matrix[row, col] = 1
            else:
                lsn_matrix[row, col] = 0

    lsn_edges = lsn_matrix.T @ A @ lsn_matrix
    return lsn_edges

main_folder = "/home/finc/Dropbox/GitHub/nilearn_task_networks/"
top_dir = "/media/finc/Elements1/extracted_HCPWM/"

fmri_list = glob.glob(top_dir + "1*/**/tfMRI_WM_RL.nii.gz", recursive=True) # if you want all tfMRI_WM_RL.nii.gz files hanging from topdir
confounds_list = glob.glob(top_dir + "1*/**/Movement_Regressors.txt", recursive=True)

fmri_list = sorted(fmri_list)
confounds_list = sorted(confounds_list)

print(f"Number of subjects: {len(fmri_list)}")

# Load Power ROIs coordinates.

power = datasets.fetch_coords_power_2011()
power_coords = np.vstack((power.rois["x"], power.rois["y"], power.rois["z"])).T

# Create masker file.

power_spheres = input_data.NiftiSpheresMasker(
    seeds=power_coords, 
    smoothing_fwhm=6, 
    radius=5,
    detrend=True, 
    standardize=True,
    low_pass=0.08, 
    high_pass=0.009,
    t_r=0.72
)

parcellation = power_spheres

# Extract timeseries. 
# Create empty timeseries to store array.
timeseries_all = np.zeros((len(fmri_list), 405, 264)) 

for sub in range(len(fmri_list)):
             
    # Load confound file.
    confounds = np.loadtxt(confounds_list[sub])
        
    # Extract timeseries.
    timeseries = parcellation.fit_transform(fmri_list[sub], confounds = confounds)
    
    # Save timeseries to array.
    timeseries_all[sub, :, :] = timeseries
    
timeseries_all.shape

# Plot.
plt.figure(figsize = (15, 5))

_ = plt.plot(timeseries_all[0,:,:])

# Generate design matrix.

t_r = 0.72
n_scans = 405

onsets_dir = "/home/finc/Dropbox/GitHub/nilearn_task_networks/support/onsets_HCP.csv"

events = pd.read_csv(onsets_dir)
events

frame_times = np.arange(n_scans) * t_r
frame_times

design_matrix = make_design_matrix(frame_times, events, hrf_model = None)
design_matrix = design_matrix.reset_index()

plt.plot(design_matrix["0back"])
plt.plot(design_matrix["2back"])
plt.legend()

# Calculate correlation matrices.
conditions  = ["0back", "2back"]
sub_n = timeseries_all.shape[0]
rois_n = timeseries_all.shape[2]

correlation_matrices = np.zeros((sub_n, len(conditions), rois_n, rois_n))

for sub in range(sub_n):
    for i, cond in enumerate(conditions):
    
        task = timeseries_all[sub, design_matrix[cond].astype("bool").astype("bool"), :]
        
        correlation_measure = ConnectivityMeasure(kind="correlation")
    
        fc = correlation_measure.fit_transform([task])[0]
        np.fill_diagonal(fc, 0)
        
        correlation_matrices[sub, i, :, :] = fc

correlation_matrices.shape

# Calculate edgewise distances.
zero_back = sym_matrix_to_vec(correlation_matrices[:,0,:,:], discard_diagonal = True)
two_back = sym_matrix_to_vec(correlation_matrices[:,1,:,:], discard_diagonal = True)

stat, pvalues = stats.ttest_rel(zero_back, two_back)

import statsmodels.stats.multitest as ssm

_, pvals_corrected, _, _ = ssm.multipletests(pvalues, alpha = 0.05, method = 'fdr_bh')

pvals_corrected_thr = np.zeros((len(pvals_corrected)))

pvals = np.array([0 if p >= 0.05 else 1 for p in pvals_corrected])
sum(pvals)

wei_vector = stat * pvals 
diag = np.zeros((264))

matrix_wei = vec_to_sym_matrix(wei_vector, diagonal = diag)
matrix_bin = vec_to_sym_matrix(pvals, diagonal = diag)
plotting.plot_matrix(matrix_wei)

