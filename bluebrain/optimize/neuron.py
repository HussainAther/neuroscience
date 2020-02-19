import bluepyopt as bpopt
import bluepyopt.ephys as ephys
import efel
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

# Create cell model.
simple_cell = ephys.models.CellModel(
        name="simple_cell",
        morph=morph,
        mechs=[hh_mech],
        params=[gna_par, gk_par, gl_par, cm])

# Print description.
print(simple_cell)

# Set up evaluator.
soma_loc = ephys.locations.NrnSeclistCompLocation(
        name="soma",
        seclist_name="somatic",
        sec_index=0,
        comp_x=0.5)

# Protocol
sweep_protocols = []

for protocol_name, amp, dur in [("neg_step", -0.01, 3000), ("pos_step", 0.015, 200)]:
    
    stim = ephys.stimuli.NrnSquarePulse(
                step_amplitude=amp,
                step_delay=250,
                step_duration=dur,
                location=soma_loc,
                total_duration=500+dur)
    
    rec = ephys.recordings.CompRecording(
            name="%s.soma.v" % protocol_name,
            location=soma_loc,
            variable="v")
    
    protocol = ephys.protocols.SweepProtocol(protocol_name, [stim], [rec])
    sweep_protocols.append(protocol)

twostep_protocol = ephys.protocols.SequenceProtocol("twostep", protocols=sweep_protocols)

# Instantiate the simulator.
nrn = ephys.simulators.NrnSimulator()

# Set values for the free parameters (they have to be within bounds).
default_params = {"gna_soma": 0.07, "gk_soma": 0.03}

# Run the simulation.
responses = twostep_protocol.run(cell_model=simple_cell, param_values=default_params, sim=nrn)

# Plot.
def plot_responses(responses):
    fig1, ax = plt.subplots(len(responses), 2, sharey = True, sharex = "row")
    ax[0,0].plot(exp_neg_trace[:, 0], exp_neg_trace[:,1], color = "k")            
    ax[0,1].plot(responses["neg_step.soma.v"]["time"], responses["neg_step.soma.v"]["voltage"])
    ax[0,0].set_title("Experiment")
    ax[0,1].set_title("Model")
    ax[1,0].plot(exp_pos_trace[:, 0], exp_pos_trace[:,1], color = "k")            
    ax[1,1].plot(responses["pos_step.soma.v"]["time"], responses["pos_step.soma.v"]["voltage"])
    fig1.tight_layout()

plot_responses(responses)

# Set values for the free parameters.
default_params = {"gna_soma": 0.15, "gk_soma": 0.03}

# Run the simulation.
responses = twostep_protocol.run(cell_model=simple_cell, param_values=default_params, sim=nrn)

# Run the simulation and plot the responses.
plot_responses(responses)

# Get the time and voltage array from the experimental data.
time = exp_pos_trace[:,0]
voltage = exp_pos_trace[:,1]

# Define the trace dictionary for efel.
trace = {"T": time, "V": voltage, "stim_start": [250], "stim_end": [450]}

# Extract features.
feature_values = efel.getFeatureValues([trace], ["Spikecount"])[0]

print("Number of spikes in the experimental trace: %s" % feature_values["Spikecount"])

feat_means = {"neg_step": {"Spikecount": 0}, "pos_step": {"Spikecount": 10}}

features = []
objectives = []

for protocol in sweep_protocols:   
    stim_start = protocol.stimuli[0].step_delay
    stim_end = stim_start + protocol.stimuli[0].step_duration

    for efel_feature_name, mean in feat_means[protocol.name].items():
        feature_name = "%s.%s" % (protocol.name, efel_feature_name)
        
        feature = ephys.efeatures.eFELFeature(
                    feature_name,
                    efel_feature_name=efel_feature_name,
                    recording_names={"": "%s.soma.v" % protocol.name},
                    stim_start=stim_start,
                    stim_end=stim_end,
                    exp_mean=mean,
                    exp_std=0.05 * abs(mean) if mean != 0 else 1)
        
        features.append(feature)
        objective = ephys.objectives.SingletonObjective(
            feature_name,
            feature)
        objectives.append(objective)

# Score function
score_calc = ephys.objectivescalculators.ObjectivesCalculator(objectives)

# Combine them together.
cell_evaluator = ephys.evaluators.CellEvaluator(
        cell_model = simple_cell,
        param_names = ["gna_soma", "gk_soma"],
        fitness_protocols = {twostep_protocol.name: twostep_protocol},
        fitness_calculator = score_calc,
        sim = nrn)

# Get score.
objectives = cell_evaluator.evaluate_with_dicts(default_params)
print("Scores:")
print(objectives)

def plot_objectives(objectives):
    ytick_pos = [x + 0.5 for x in range(len(objectives))]
    obj_val = objectives.values()
    obj_keys = objectives.keys()
    fig, ax = plt.subplots(figsize = (6,2.5))
    ax.barh(ytick_pos,
              obj_val,
              height=0.5,
              align="center",
              color="blue",
              alpha=0.5)
    ax.axvline(3, color = "gray", ls = "--")
    ax.set_yticks(ytick_pos)
    ax.set_ylabel("Score name")
    ax.set_xlabel("Score value (# STD)")
    ax.set_yticklabels(obj_keys)
    fig.tight_layout()
plot_objectives(objectives)  

# Run optimization.
optimization = bpopt.deapext.optimisations.IBEADEAPOptimisation(                              
        evaluator=cell_evaluator,                                                
        offspring_size = 10)    

final_pop, hall_of_fame, logs, hist = optimisation.run(max_ngen=13)

# Optimization results
fig, axs = plt.subplots(1,2, figsize = (5,4))
      
ax0 = axs[0].twinx()
ax0.set_ylabel("gk_soma (S/cm$^2$)")
for param_set in hall_of_fame:
    axs[0].plot(0, param_set[0], "o") 
    axs[0].set_ylabel("gna_soma (S/cm$^2$)")
    ax0.plot(1, param_set[1], "o")
    axs[0].set_title("Best individuals")

ax1 = axs[1].twinx()
axs[1].set_ylabel("gna_soma (S/cm$^2$)")
ax1.set_ylabel("gk_soma (S/cm$^2$)")
for param_set in hist.genealogy_history.values()[:100]:
        axs[1].plot(0,param_set[0], "o", color = "gray")
        ax1.plot(1,param_set[1], "o", color = "gray")
        axs[1].set_title("Other individuals")
      
for ax in axs:
    ax.set_xlim([-1,2])
    ax.set_xticks([])
    
fig.tight_layout()
best_ind = hall_of_fame[0]
print("Best individual: {} ".format(best_ind))
best_ind_dict = cell_evaluator.param_dict(best_ind)
print(best_ind_dict)
responses = twostep_protocol.run(cell_model=simple_cell, param_values=best_ind_dict, sim=nrn)
plot_responses(responses)
objectives = score_calc.calculate_scores(responses)
print("Score: ", objectives)
plot_objectives(objectives)

gen_numbers = logs.select("gen")
min_fitness = logs.select("min")
max_fitness = logs.select("max")

plt.plot(gen_numbers, min_fitness, label="min fitness", color = "k")
plt.xlabel("generation #")
plt.ylabel("score (# std)")
plt.legend()
plt.xlim(min(gen_numbers) - 1, max(gen_numbers) + 1) 
plt.ylim(0.9*min(min_fitness), 1.1 * max(min_fitness));