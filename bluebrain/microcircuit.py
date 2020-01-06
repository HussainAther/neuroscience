import os
import urllib # Download files from the web
import neurom # Analyse / view morphologies
import neurom.viewer
import zipfile # Extract zip files

"""
The Neocortical Microcircuit (microcircuit) Collaboration Portal (NMC Portal, 
at [https://bbp.epfl.ch/nmc-portal](https://bbp.epfl.ch/nmc-portal)) 
provides an online public resource of the Blue Brain Project's 
first release of a digital reconstruction of the microcircuitry 
of juvenile Rat somatosensory cortex, access to experimental 
data sets used in the reconstruction, and the resulting single cell models.
"""

# Download neocortical layer 5 thick tufted pyramidal cell model.
urllib.urlretrieve("https://bbp.epfl.ch/nmc-portal/documents/10184/1921755/L5_TTPC2_cADpyr232_1.zip/a058fc9c-6c67-417b-a65b-20742902ccbb','L5_TTPC2_cADpyr232_1.zip")

# Extract the zip file.
with zipfile.ZipFile("L5_TTPC2_cADpyr232_1.zip", "r") as zip_file:
    zip_file.extractall(".")
