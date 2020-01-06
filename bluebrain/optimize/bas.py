import bluepyopt as bpop
import bluepyopt.ephys as ephys
import efel
import matplotlib.pyplot as plt
import neurom
import neurom.view
import numpy

"""
Optimize parameters for a ball-and-stick model (ball and stick optimization)
"""

morph_swc_string = """
1 1 0.0 0.0 -10.0 10.0 -1                                                        
2 1 0.0 0.0 0.0 10.0 1                                                           
3 1 0.0 0.0 10.0 10.0 2                                                          
4 3 0.0 10.0 0.0 2.0 2                                                           
5 3 0.0 110.0 0.0 2.0 4
"""

with open("ballandstick.swc", "w") as swc_file:
    swc_file.write(morph_swc_string)

# Draw the ball-and-stick model.
fig, ax = neurom.viewer.draw(neurom.load_neuron("ballandstick.swc"))

# Get the morphology.
morph = ephys.morphologies.NrnFileMorphology("ballandstick.swc")

# Get the morphology by location.
somatic_loc = ephys.locations.NrnSeclistLocation("somatic", seclist_name="somatic")
dend_loc = ephys.locations.NrnSeclistLocation("basal", seclist_name="basal")

cm = ephys.
parameters.NrnSectionParameter(
        name="cm",
        param_name="cm",
        value=1.0, # in microfarad/cm2
        locations=[somatic_loc, dend_loc],
        frozen=True)

# Fix leak conductance dendrite.
gl_dend = ephys.parameters.NrnSectionParameter(
        name="gl_dend",
        param_name="gl_hh",
        value=1e-5,
        locations=[dend_loc],
        frozen=True)

# Disable Na and K.
gnabar_dend = ephys.parameters.NrnSectionParameter(                                    
        name="gnabar_hh_dend",
        param_name="gnabar_hh",
        locations=[dend_loc],
        value=0,
        frozen=True)
     
gkbar_dend = ephys.parameters.NrnSectionParameter(
        name="gkbar_hh_dend",
        param_name="gkbar_hh",
        value=0,
        locations=[dend_loc],
        frozen=True)

# Set boundaries.
gnabar_soma = ephys.parameters.NrnSectionParameter(                                    
        name="gnabar_soma",
        param_name="gnabar_hh",
        locations=[somatic_loc],
        bounds=[0.0, 1.0],
        frozen=False)     
gkbar_soma = ephys.parameters.NrnSectionParameter(
        name="gkbar_soma",
        param_name="gkbar_hh",
        bounds=[0.0, 1.0],
        locations=[somatic_loc],
        frozen=False)

# Compute the cell model.
ballandstick_cell = ephys.models.CellModel(
        name="simple_cell",
        morph=morph,
        mechs=[hh_mech],
        params=[cm, gnabar_dend, gkbar_dend, gl_dend, gnabar_soma, gkbar_soma])  

print(ballandstick_cell)

# Create protocols.
soma_loc = ephys.locations.NrnSeclistCompLocation(
        name="soma",
        seclist_name="somatic",
        sec_index=0,
        comp_x=0.5)

# Current clamp with square current inject
sweep_protocols = []
for protocol_name, amplitude in [("step1", 0.1), ("step2", 0.5)]:
    stim = ephys.stimuli.NrnSquarePulse(
                step_amplitude=amplitude,
                step_delay=100,
                step_duration=50,
                location=soma_loc,
                total_duration=200)
    rec = ephys.recordings.CompRecording(
            name="%s.soma.v" % protocol_name,
            location=soma_loc,
            variable="v")
    protocol = ephys.protocols.SweepProtocol(protocol_name, [stim], [rec])
    sweep_protocols.append(protocol)
twostep_protocol = ephys.protocols.SequenceProtocol("twostep", protocols=sweep_protocols)

# Run protocol.
nrn = ephys.simulators.NrnSimulator()

default_params = {"gnabar_soma": 0.25, "gkbar_soma": 0.1}
responses = twostep_protocol.run(cell_model=ballandstick_cell, param_values=default_params, sim=nrn)

# Plot.
def plot_responses(responses):
    """
    Plot the response traces to the current.
    """
    plt.subplot(2,1,1)
    plt.plot(responses["step1.soma.v"]["time"], responses["step1.soma.v"]["voltage"], label="step1")
    plt.legend()
    plt.subplot(2,1,2)
    plt.plot(responses["step2.soma.v"]["time"], responses["step2.soma.v"]["voltage"], label="step2")
    plt.legend()
    plt.tight_layout()

plot_responses(responses)

# Test out other parameters.
other_params = {"gnabar_soma": 0.1, "gkbar_soma": 0.1}
plot_responses(twostep_protocol.run(cell_model=ballandstick_cell, param_values=other_params, sim=nrn))

# Use eFEL (eFeature Extraction Library) to analyze traces.
responses = twostep_protocol.run(cell_model=ballandstick_cell, param_values=default_params, sim=nrn)

# Show where these names come from.
step2_time = responses["step2.soma.v"]["time"]
step2_voltage = responses["step2.soma.v"]["voltage"]

# Define this dictionary.
trace = {"T": step2_time, "V": step2_voltage, "stim_start": [100], "stim_end": [150]}

# Explain AP_width (from where to where is AP_amplitude...
feature_values = efel.getFeatureValues([trace], ["Spikecount", "AP_width", "AP_amplitude"])[0]

# Plot.
plot_responses(responses)
print("Number of spikes in 2nd trace: %s" % feature_values["Spikecount"])
print("Spike widths (ms) in 2nd trace: %s" % feature_values["AP_width"])
print("Spike amplitude (mV) in 2nd trace: %s" % feature_values["AP_amplitude"])

# Extract features.
efel_feature_means = {"step1": {"Spikecount": 4}, "step2": {"Spikecount": 6}}

objectives = []
features = []

# Define the objects of the protocol.
for protocol in sweep_protocols:
    stim_start = protocol.stimuli[0].step_delay
    stim_end = stim_start + protocol.stimuli[0].step_duration
    for efel_feature_name, mean in efel_feature_means[protocol.name].items():
        feature_name = "%s.%s" % (protocol.name, efel_feature_name)
        feature = ephys.efeatures.eFELFeature(
                    feature_name,
                    efel_feature_name=efel_feature_name,
                    recording_names={"": "%s.soma.v" % protocol.name},
                    stim_start=stim_start,
                    stim_end=stim_end,
                    exp_mean=mean,
                    exp_std=0.05 * abs(mean))
        features.append(feature)
        objective = ephys.objectives.SingletonObjective(
            feature_name,
            feature)
        objectives.append(objective)

# Create cell evaluator.
score_calc = ephys.objectivescalculators.ObjectivesCalculator(objectives) 

# Evaluate.
cell_evaluator = ephys.evaluators.CellEvaluator(
        cell_model=ballandstick_cell,
        param_names=["gnabar_soma", "gkbar_soma"],
        fitness_protocols={twostep_protocol.name: twostep_protocol},
        fitness_calculator=score_calc,
        sim=nrn)

print("Scores:", cell_evaluator.evaluate_with_dicts(default_params))

# Run the optimization algorithm.
optimization_algorithm = bpop.deapext.optimisations.IBEADEAPOptimisation(
        evaluator=cell_evaluator,
        offspring_size = 10)

# Run it for 10 generations.
final_pop, hall_of_fame, logs, hist = optimisation_algorithm.run(max_ngen=10)

# Get the best ones (Hall of fame).
for ind in hall_of_fame:
    print("gnabar_soma=%f, gkbar_soma=%f" % tuple(ind))

best_ind = hall_of_fame[0]
print("Best individual:  ", best_ind)
best_ind_dict = cell_evaluator.param_dict(best_ind)
print(best_ind_dict)

# Plot them.
responses = twostep_protocol.run(cell_model=ballandstick_cell, param_values=best_ind_dict, sim=nrn)
print("Score: ", score_calc.calculate_scores(responses))
plot_responses(responses)

# Get the optimization metrics. 
gen_numbers = logs.select("gen")
min_fitness = logs.select("min")
max_fitness = logs.select("max")
plt.plot(gen_numbers, min_fitness, label="min fitness")
plt.xlabel("generation #")
plt.ylabel("score (# std)")
plt.legend()
plt.xlim(min(gen_numbers) - 1, max(gen_numbers) + 1) 
plt.ylim(0.9*min(min_fitness), 1.1 * max(min_fitness)) 
