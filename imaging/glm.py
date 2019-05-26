import numpy as np
import matplotlib.pyplot as plt
import os
import requests
import zipfile
import pandas as pd
import nibabel

"""
Generalized linear model (GLM general) for fMRI.
"""

# Define the URL of the data and download it using the Requests libary
url = "http://www.fil.ion.ucl.ac.uk/spm/download/data/MoAEpilot/MoAEpilot.zip"
data = requests.get(url)

# Check if the targed folder for storing the data already exists. If not create it and save the zip file.
if os.path.exists("./fMRI_data") == False:
    os.mkdir("fMRI_data")
    
open("./fMRI_data/data.zip", "wb").write(data.content)

# Unzip the file
zip_ref = zipfile.ZipFile("./fMRI_data/data.zip", "r")
zip_ref.extractall("./fMRI_data/")
zip_ref.close()

# Find all files in the structural data folder
data_path = "./fMRI_data/sM00223/"
files = os.listdir(data_path)

# Read in the data
data_all = []
for data_file in files:
    if data_file[-3:] == "hdr":
        data = nibabel.load(data_path + data_file).get_data()

data = np.rot90(data.squeeze(), 1)

# Plot
fig, ax = plt.subplots(1, 6, figsize=[18, 3])

n = 0
slice = 0
for _ in range(6):
    ax[n].imshow(data[:, :, slice], 'gray')
    ax[n].set_xticks([])
    ax[n].set_yticks([])
    ax[n].set_title('Slice number: {}'.format(slice), color='r')
    n += 1
    slice += 10
    
fig.subplots_adjust(wspace=0, hspace=0)
plt.show()

# Basic information about the data acquisition
x_size = 64
y_size = 64
n_slice = 64
n_volumes = 96

# Find all files in the data folder
data_path = "./fMRI_data/fM00223/"
files = os.listdir(data_path)

# Read in the data and organize it with respect to the acquisition parameters
data_all = []
for data_file in files:
    if data_file[-3:] == "hdr":
        data = nibabel.load(data_path + data_file).get_data()        
        data_all.append(data.reshape(x_size, y_size, n_slice))

# Create a 3x6 subplot 
fig, ax = plt.subplots(3, 6, figsize=[18, 11])

# Organize the data for visualisation in the coronal plane
coronal = np.transpose(data_all, [1, 3, 2, 0])
coronal = np.rot90(coronal, 1)

# Organize the data for visualisation in the transversal plane
transversal = np.transpose(data_all, [2, 1, 3, 0])
transversal = np.rot90(transversal, 2)

# Organize the data for visualisation in the sagittal plane
sagittal = np.transpose(data_all, [2, 3, 1, 0])
sagittal = np.rot90(sagittal, 1)

# Plot some of the images in different planes
n = 10
for i in range(6):
    ax[0][i].imshow(coronal[:, :, n, 0], cmap="gray")
    ax[0][i].set_xticks([])
    ax[0][i].set_yticks([])
    if i == 0:
        ax[0][i].set_ylabel("coronal", fontsize=25, color="r")
    n += 10
    
n = 5
for i in range(6):
    ax[1][i].imshow(transversal[:, :, n, 0], cmap="gray")
    ax[1][i].set_xticks([])
    ax[1][i].set_yticks([])
    if i == 0:
        ax[1][i].set_ylabel("transversal", fontsize=25, color="r")
    n += 10
    
n = 5
for i in range(6):
    ax[2][i].imshow(sagittal[:, :, n, 0], cmap="gray")
    ax[2][i].set_xticks([])
    ax[2][i].set_yticks([])
    if i == 0:
        ax[2][i].set_ylabel("sagittal", fontsize=25, color="r')
    n += 10

fig.subplots_adjust(wspace=0, hspace=0)
plt.show()

# Create an empty plot with defined aspect ratio
fig, ax = plt.subplots(1, 1, figsize=[18, 5])

# Plot the timecourse of a random voxel
ax.plot(transversal[30, 30, 35, :], lw=3)
ax.set_xlim([0, transversal.shape[3]-1])
ax.set_xlabel("time [s]", fontsize=20)
ax.set_ylabel("signal strength", fontsize=20)
ax.set_title("voxel time course", fontsize=25)
ax.tick_params(labelsize=12)

plt.show()

# Rearrange and reshape data for export
data_all = np.transpose(data_all, [3, 2, 1, 0])
data_all = np.reshape(data_all, [n_slice, y_size*x_size, n_volumes])

# Check if output path exists, if not create it.
if os.path.exists(".fMRI_data/csv_data") == False:
    os.mkdir("./fMRI_data/csv_data")

# Export each slice as a .csv file 
n = 0
for export in data_all:
    save_file = "slice_{}.csv".format(n)
    save_path = "./fMRI_data/csv_data/{}".format(save_file)
    pd.DataFrame(export).to_csv(save_path, header=False, index=False)
    n += 1

# Main parameters of the fMRI scan and experimental desgin
block_design = ["rest", "stim"]
block_size = 6
block_RT = 7
block_total = 16
block_length = block_size*block_RT
acq_num = block_size*block_total
data_time = block_length*block_total
data_time_vol = np.arange(acq_num)*block_RT
x_size = 64
y_size = 64

# Reshape the data
data_ordered = data.reshape(x_size, y_size, acq_num)

# Calculate the mean signal for each voxel
mean_data = data_ordered.mean(axis=2)

# Create the design matrix
constant = np.ones(acq_num)
rest = np.zeros(block_size)
stim = np.ones(block_size)
block = np.concatenate((rest, stim), axis=0)
predicted_response = np.tile(block, int(block_total/2))

design_matrix = np.array((constant, predicted_response))

# Create the plots
fig, ax = plt.subplots(2,1, figsize=(15, 5))
ax[0].plot(design_matrix[0], lw=3)
ax[0].set_xlim(0, acq_num-1)
ax[0].set_ylim(0, 1.5)
ax[0].set_title('constant', fontsize=25)
ax[0].set_xticks([])
ax[0].set_yticks([0,1])
ax[0].tick_params(labelsize=12)
ax[0].tick_params(labelsize=12)
ax[1].plot(design_matrix[1], lw=3)
ax[1].set_xlim(0, acq_num-1)
ax[1].set_ylim(0, 1.5)
ax[1].set_title('expected response', fontsize=25)
ax[1].set_yticks([0,1])
ax[1].set_xlabel('time [volumes]', fontsize=20)
ax[1].tick_params(labelsize=12)
ax[1].tick_params(labelsize=12)
fig.subplots_adjust(wspace=0, hspace=0.5)
plt.show()


def do_GLM(X, y):
    """
    For design matrix X and expected response y, perform the generalized
    linear model (GLM) fit and error assessment.
    """
    # Make sure design matrix has the right orientation
    if X.shape[1] > X.shape[0]:
        X = X.transpose()
    # Calculate the dot product of the transposed design matrix and the design matrix
    # and invert the resulting matrix.
    tmp = np.linalg.inv(X.transpose().dot(X))
    # Now calculate the dot product of the above result and the transposed design matrix
    tmp = tmp.dot(X.transpose())
    # Pre-allocate variables
    beta = np.zeros((y.shape[0], X.shape[1]))
    e = np.zeros(y.shape)
    model = np.zeros(y.shape)
    r = np.zeros(y.shape[0])
    # Find beta values for each voxel and calculate the model, error and the correlation coefficients 
    for i in range(y.shape[0]):
        beta[i] = tmp.dot(y[i,:].transpose())
        model[i] = X.dot(beta[i])
        e[i] = (y[i,:] - model[i])
        r[i] = np.sqrt(model[i].var()/y[i,:].var())
    return beta, model, e, r

beta, model, e, r = do_GLM(design_matrix, data)
r = r.reshape(x_size,y_size)
map = r.copy()
map[map<0.35] = np.nan
fig, ax = plt.subplots(1,3,figsize=(18, 6))
ax[0].imshow(mean_data, cmap="gray")
ax[0].set_title("1st EPI image", fontsize=25)
ax[0].set_xticks([])
ax[0].set_yticks([])
ax[1].imshow(r, cmap="afmhot")
ax[1].set_title("un-thresholded map", fontsize=25)
ax[1].set_xticks([])
ax[1].set_yticks([])
ax[2].imshow(mean_data, cmap="gray")
ax[2].imshow(map, cmap="afmhot")
ax[2].set_title("thresholded map (overlay)", fontsize=25)
ax[2].set_xticks([])
ax[2].set_yticks([])
plt.show()

def scale(data):
    """
    Scale data data of an array using the max and min values.
    """
    return (data - data.min()) / (data.max() - data.min())

def z_score(data):
    """
    Calculate z-score Z-score Z score.
    """
    mean = data.mean(axis=1, keepdims=True)
    std = data.std(axis=1, keepdims=True)
    norm_data = (data-mean)/std
    return norm_data

avg = z_score(data_ordered[~np.isnan(map),:])
avg = np.transpose(avg).mean(axis=1)
avg = scale(avg)
design_matrix = np.array((constant, predicted_response, avg))

# Create the plots
fig, ax = plt.subplots(3,1, figsize=(15, 5))
ax[0].plot(design_matrix[0], lw=3)
ax[0].set_xlim(0, acq_num-1)
ax[0].set_ylim(0, 1.5)
ax[0].set_title("constant", fontsize=25)
ax[0].set_xticks([])
ax[0].set_yticks([0,1])
ax[0].tick_params(labelsize=12)
ax[0].tick_params(labelsize=12)
ax[1].plot(design_matrix[1], lw=3)
ax[1].set_xlim(0, acq_num-1)
ax[1].set_ylim(0, 1.5)
ax[1].set_title("expected response", fontsize=25)
ax[1].set_yticks([0,1])
ax[1].set_xticks([])
ax[1].tick_params(labelsize=12)
ax[1].tick_params(labelsize=12)
ax[2].plot(design_matrix[2], lw=3)
ax[2].set_xlim(0, acq_num-1)
ax[2].set_ylim(0, 1.5)
ax[2].set_title("mean response", fontsize=25)
ax[2].set_yticks([0,1])
ax[2].set_xlabel("time [volumes]", fontsize=20)
ax[2].tick_params(labelsize=12)
ax[2].tick_params(labelsize=12)
fig.subplots_adjust(wspace=0, hspace=0.5)
plt.show()
