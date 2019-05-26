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

