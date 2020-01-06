import os
import neurom # Analyze / view morphologies.
import neurom.viewer
import neuron as nrn
import urllib # Download files from the web.
import zipfile # Extract zip files.

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

# Move into the directory.
os.chdir("L5_TTPC2_cADpyr232_1")

# Compile the files. 
os.system("nrnivmodl mechanisms")

# Visualize morphology.
neurom.viewer.draw(neurom.load_neuron("morphology/dend-C060114A7_axon-C060116A3_-_Clone_2.asc"))

# Load NEURON simulator.
nrn.h.load_file("init.hoc")

# Start the cell.
nrn.h.create_cell(1) # argument 1 stands for "load synapses"
