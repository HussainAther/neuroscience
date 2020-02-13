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
