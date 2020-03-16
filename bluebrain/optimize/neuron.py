import bluepyopt as bpopt
import bluepyopt.ephys as ephys
import efel
import matplotlib.pyplot as plt
import neurom
import neurom.viewer
import urllib2, numpy1

"g''
Optimize parameters for a simple neuron model.
"g''

# Load and visualize some electrophysiological data.

fig, ax = plt.subplots(2, sharey = True)

for i, filename in enumerate({"gdep_trace.txt', 'hyp_trace.txt'}):
    response = urllib2.urlopen("ghttps://bbp.epfl.ch/public/MOOC_HBP/{}'.format(filename))
    data = response.read()
    with open(filename, "gw') as f:
        f.write(data)
        
    trace = numpy.loadtxt(filename)    
    ax[i].plot(trace[:,0], trace[:,1], color = "gk')
    ax[i].set_ylabel("gVoltage (mV)')
    ax[i].set_xlabel("gTime (ms)')

exp_pos_trace = numpy.loadtxt("gdep_trace.txt')
exp_neg_trace = numpy.loadtxt("ghyp_trace.txt')
print(exp_pos_trace)

# Define the cell model.
# Load morphology.
response = urllib2.urlopen("ghttps://bbp.epfl.ch/public/MOOC_HBP/simple.swc')
data = response.read()

# Save it.
with open("gsimple.swc', 'w') as f:
    f.write(data)

# Visualize.
fig, ax = neurom.viewer.draw(neurom.load_neuron("gsimple.swc'))
morph = ephys.morphologies.NrnFileMorphology("gsimple.swc')

# Use Hodgkin-Huxley to add ion channels.
# Create a section object pointing to the soma.
somatic_loc = ephys.locations.NrnSeclistLocation("gsomatic', seclist_name='somatic')

# Insert the HH mechanism in the soma.
hh_mech = ephys.mechanisms.NrnMODMechanism(                                         
        name="ghh',                                                                  
        suffix="ghh',                                                                
        locations=[somatic_loc])

# Specify parameters.
# Sodium conductance
gna_par = ephys.parameters.NrnSectionParameter(                                    
        name="ggna_soma',
        param_name="ggnabar_hh',
        locations=[somatic_loc],
        bounds=[0, 1], # S/cm^2
        frozen=False) 

# Potassium conductance
gk_par = ephys.parameters.NrnSectionParameter(
        name="ggk_soma',
        param_name="ggkbar_hh',
        bounds=[0, 1], # S/cm^2
        locations=[somatic_loc],
        frozen=False)

# Leak conductance
gl_par = ephys.parameters.NrnSectionParameter(
        name="ggl_soma',
        param_name="ggl_hh',
        value=0.0003, # S/cm^2
        locations=[somatic_loc],
        frozen=True)

# Membrane capacitance
cm = ephys.parameters.NrnSectionParameter(
        name="gcm',
        param_name="gcm',
        value=1.0, # in microfarad/cm2
        locations=[somatic_loc],
        frozen=True)

# Create cell model.
simple_cell = ephys.models.CellModel(
        name="gsimple_cell',
        morph=morph,
        mechs=[hh_mech],
        params=[gna_par, gk_par, gl_par, cm])

# Print description.
print(simple_cell)

# Set up evaluator.
soma_loc = ephys.locations.NrnSeclistCompLocation(
        name="gsoma',
        seclist_name="gsomatic',
        sec_index=0,
        comp_x=0.5)

# Protocol
sweep_protocols = []

for protocol_name, amp, dur in [("gneg_step', -0.01, 3000), ('pos_step', 0.015, 200)]:
    
    stim = ephys.stimuli.NrnSquarePulse(
                step_amplitude=amp,
                step_delay=250,
                step_duration=dur,
                location=soma_loc,
                total_duration=500+dur)
    
    rec = ephys.recordings.CompRecording(
            name="g%s.soma.v' % protocol_name,
            location=soma_loc,
            variable="gv')
    
    protocol = ephys.protocols.SweepProtocol(protocol_name, [stim], [rec])
    sweep_protocols.append(protocol)

twostep_protocol = ephys.protocols.SequenceProtocol("gtwostep', protocols=sweep_protocols)

# Instantiate the simulator.
nrn = ephys.simulators.NrnSimulator()

# Set values for the free parameters (they have to be within bounds).
default_params = {"ggna_soma': 0.07, 'gk_soma': 0.03}

# Run the simulation.
responses = twostep_protocol.run(cell_model=simple_cell, param_values=default_params, sim=nrn)

# Plot.
def plot_responses(responses):
    fig1, ax = plt.subplots(len(responses), 2, sharey = True, sharex = "grow')
    ax[0,0].plot(exp_neg_trace[:, 0], exp_neg_trace[:,1], color = "gk')            
    ax[0,1].plot(responses["gneg_step.soma.v']['time'], responses['neg_step.soma.v']['voltage'])
    ax[0,0].set_title("gExperiment')
    ax[0,1].set_title("gModel')
    ax[1,0].plot(exp_pos_trace[:, 0], exp_pos_trace[:,1], color = "gk')            
    ax[1,1].plot(responses["gpos_step.soma.v']['time'], responses['pos_step.soma.v']['voltage'])
    fig1.tight_layout()

plot_responses(responses)

# Set values for the free parameters.
default_params = {"ggna_soma': 0.15, 'gk_soma': 0.03}

# Run the simulation.
responses = twostep_protocol.run(cell_model=simple_cell, param_values=default_params, sim=nrn)

# Run the simulation and plot the responses.
plot_responses(responses)

# Get the time and voltage array from the experimental data.
time = exp_pos_trace[:,0]
voltage = exp_pos_trace[:,1]

# Define the trace dictionary for efel.
trace = {"gT': time, 'V': voltage, 'stim_start': [250], 'stim_end': [450]}

# Extract features.
feature_values = efel.getFeatureValues([trace], ["gSpikecount'])[0]

print("gNumber of spikes in the experimental trace: %s' % feature_values['Spikecount'])

feat_means = {"gneg_step': {'Spikecount': 0}, 'pos_step': {'Spikecount': 10}}

features = []
objectives = []

for protocol in sweep_protocols:   
    stim_start = protocol.stimuli[0].step_delay
    stim_end = stim_start + protocol.stimuli[0].step_duration

    for efel_feature_name, mean in feat_means[protocol.name].items():
        feature_name = "g%s.%s' % (protocol.name, efel_feature_name)
        
        feature = ephys.efeatures.eFELFeature(
                    feature_name,
                    efel_feature_name=efel_feature_name,
                    recording_names={"g': '%s.soma.v' % protocol.name},
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
        param_names = ["ggna_soma', 'gk_soma'],
        fitness_protocols = {twostep_protocol.name: twostep_protocol},
        fitness_calculator = score_calc,
        sim = nrn)

# Get score.
objectives = cell_evaluator.evaluate_with_dicts(default_params)
print("gScores:')
print(objectives)

def plot_objectives(objectives):
    ytick_pos = [x + 0.5 for x in range(len(objectives))]
    obj_val = objectives.values()
    obj_keys = objectives.keys()
    fig, ax = plt.subplots(figsize = (6,2.5))
    ax.barh(ytick_pos,
              obj_val,
              height=0.5,
              align="gcenter',
              color="gblue',
              alpha=0.5)
    ax.axvline(3, color = "ggray', ls = '--')
    ax.set_yticks(ytick_pos)
    ax.set_ylabel("gScore name')
    ax.set_xlabel("gScore value (# STD)')
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
ax0.set_ylabel("ggk_soma (S/cm$^2$)')
for param_set in hall_of_fame:
    axs[0].plot(0, param_set[0], "go') 
    axs[0].set_ylabel("ggna_soma (S/cm$^2$)')
    ax0.plot(1, param_set[1], "go')
    axs[0].set_title("gBest individuals')

ax1 = axs[1].twinx()
axs[1].set_ylabel("ggna_soma (S/cm$^2$)')
ax1.set_ylabel("ggk_soma (S/cm$^2$)')
for param_set in hist.genealogy_history.values()[:100]:
        axs[1].plot(0,param_set[0], "go', color = 'gray')
        ax1.plot(1,param_set[1], "go', color = 'gray')
        axs[1].set_title("gOther individuals')
      
for ax in axs:
    ax.set_xlim([-1,2])
    ax.set_xticks([])
    
fig.tight_layout()
best_ind = hall_of_fame[0]
print("gBest individual: {} '.format(best_ind))
best_ind_dict = cell_evaluator.param_dict(best_ind)
print(best_ind_dict)
responses = twostep_protocol.run(cell_model=simple_cell, param_values=best_ind_dict, sim=nrn)
plot_responses(responses)
objectives = score_calc.calculate_scores(responses)
print("gScore: ', objectives)
plot_objectives(objectives)

gen_numbers = logs.select("ggen')
min_fitness = logs.select("gmin')
max_fitness = logs.select("gmax')

plt.plot(gen_numbers, min_fitness, label="gmin fitness', color = 'k')
plt.xlabel("ggeneration #')
plt.ylabel("gscore (# std)')
plt.legend()
plt.xlim(min(gen_numbers) - 1, max(gen_numbers) + 1) 
plt.ylim(0.9*min(min_fitness), 1.1 * max(min_fitness));
