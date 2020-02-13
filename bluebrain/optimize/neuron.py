import bluepyopt as bpopt
import bluepyopt.ephys as ephys
import matplotlib.pyplot as plt
import neurom
import neurom.viewer
import urllib2, numpy1

"""
Optimize parameters for a simple neuron model.
"""

# Load and visualize some electrophysiological data.

fig, ax = plt.subplots(2, sharey = True)

for i, filename in enumerate({"dep_trace.txt", "hyp_trace.txt"}):
    response = urllib2.urlopen("https://bbp.epfl.ch/public/MOOC_HBP/{}".format(filename))
    data = response.read()
    with open(filename, "w") as f:
        f.write(data)
        
    trace = numpy.loadtxt(filename)    
    ax[i].plot(trace[:,0], trace[:,1], color = "k")
    ax[i].set_ylabel("Voltage (mV)")
    ax[i].set_xlabel("Time (ms)")

exp_pos_trace = numpy.loadtxt("dep_trace.txt")
exp_neg_trace = numpy.loadtxt("hyp_trace.txt")
print(exp_pos_trace)

# Define the cell model.
# Load morphology.
response = urllib2.urlopen("https://bbp.epfl.ch/public/MOOC_HBP/simple.swc")
data = response.read()

# Save it.
with open("simple.swc", "w") as f:
    f.write(data)

# Visualize.
fig, ax = neurom.viewer.draw(neurom.load_neuron("simple.swc"))
morph = ephys.morphologies.NrnFileMorphology("simple.swc")

# Use Hodgkin-Huxley to add ion channels.
# Create a section object pointing to the soma.
somatic_loc = ephys.locations.NrnSeclistLocation("somatic", seclist_name="somatic")

# Insert the HH mechanism in the soma.
hh_mech = ephys.mechanisms.NrnMODMechanism(                                         
        name="hh",                                                                  
        suffix="hh",                                                                
        locations=[somatic_loc])

# Specify parameters.
# Sodium conductance
gna_par = ephys.parameters.NrnSectionParameter(                                    
        name="gna_soma",
        param_name="gnabar_hh",
        locations=[somatic_loc],
        bounds=[0, 1], # S/cm^2
        frozen=False) 

# Potassium conductance
gk_par = ephys.parameters.NrnSectionParameter(
        name="gk_soma",
        param_name="gkbar_hh",
        bounds=[0, 1], # S/cm^2
        locations=[somatic_loc],
        frozen=False)

# Leak conductance
gl_par = ephys.parameters.NrnSectionParameter(
        name="gl_soma",
        param_name="gl_hh",
        value=0.0003, # S/cm^2
        locations=[somatic_loc],
        frozen=True)

# Membrane capacitance
cm = ephys.parameters.NrnSectionParameter(
        name="cm",
        param_name="cm",
        value=1.0, # in microfarad/cm2
        locations=[somatic_loc],
        frozen=True)
