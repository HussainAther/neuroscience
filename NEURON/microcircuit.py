import matplotlib.pyplot as plt
import neurom # Analyze / view morphologies.
import neurom.viewer
import neuron as nrn
import numpy as np
import os
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
cell = nrn.h.cell
soma = cell.soma[0]

# Inject step current.
# Mention source of amplitude
stimulus = nrn.h.IClamp(0.5, sec=soma)
stimulus.dur = 400 # ms
stimulus.delay = 100  # ms     
stimulus.amp = 0.691907 # nA

# Get current amplitudes.
with open("current_amps.dat") as current_file:
    current_content = current_file.read()

print("File content: ", current_content)
holding_current, step1_current, step2_current, step3_current = [float(x) for x in current_content.split()]

print("Holding current: %f nA" % holding_current)
print("Step 1: %f nA" % step1_current)
print("Step 2: %f nA" % step2_current)
print("Step 3: %f nA" % step3_current)

# Activate recording of activity.
nrn.h.create_recording() 

# Run the simuation.
nrn.h.tstop = 600 # ms
nrn.h.dt = 0.05 # 
nrn.h.run()

nrn.h.save_recording()

time = nrn.h.time
voltage = nrn.h.voltage

def plot_tv(time_array, voltage_array, show=True, label=None, constants=[]):
    """
    Plot time and voltage.
    """
    plt.plot(time_array, voltage_array, label=label)
    for constant in constants:
        plt.plot(time_array, constant*np.ones(len(time_array)))
    plt.xlabel("Time (ms)")
    plt.ylabel("Membrane voltage (mV)")
    if show:
        plt.show()
    
plot_tv(time, voltage)

# Define in-vivo like stimulus with Poisson process firing neurons and
# presynaptic morphology type (m-type) specific rates.
with open("synapses/mtype_map.tsv") as mtype_map_file:
    mtype_map_content = mtype_map_file.read()
    
mtype_map = {}
for line in mtype_map_content.split("\n")[:-1]:
    n, mtype = line.split()
    mtype_map[mtype] = int(n)
    
print(mtype_map)

def init_synapses(enabled_mtypes=[]):
    """
    Enable all the synapses that are projected onto this cell 
    from mtype listed in enabled_mtypes.
    """
    enabled_mtype_ints = [mtype_map[mtype] for mtype in enabled_mtypes]
    
    for i in range(0, int(cell.synapses.n_of_mtypes)): # Loop over all the m-type
        if i in enabled_mtype_ints: # Enable synapses.
            # The [were_]active_pre_mtypes is a NEURON vector
            # (it uses the .x syntax to access the elements).
            # When the value in the vector is 1 all the presynaptic neurons
            # of a particular m-types are active (and inactive when it is 0).
            cell.synapses.were_active_pre_mtypes.x[i]= 0
            cell.synapses.active_pre_mtypes.x[i] = 1        
        else: # Disable synapses.
            cell.synapses.were_active_pre_mtypes.x[i]= 1
            cell.synapses.active_pre_mtypes.x[i] = 0
    cell.synapses.update_synapses(nrn.h.synapse_plot) # Update the synapses.

