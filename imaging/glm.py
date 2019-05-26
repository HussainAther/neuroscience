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
