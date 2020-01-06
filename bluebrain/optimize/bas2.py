import os
os.unsetenv('PYTHONHOME') # Solve an issue with NEURON simulator import
!pip install -q bluepyopt==1.5.12 matplotlib==2.0.2 numpy==1.13.0 2>&1 | grep -v 'SNIMissingWarning\|InsecurePlatformWarning'
!pip install -q --upgrade --pre -i https://bbpteam.epfl.ch/repository/devpi/simple/ single-cell-mooc-client

%matplotlib inline
import matplotlib.pyplot as plt
import bluepyopt as bpop
import bluepyopt.ephys as ephys

morph_swc_string2 = """
1 1 0.0 -10.0 0.0 10.0 -1                                                        
2 1 0.0 0.0 0.0 10.0 1                                                           
3 1 0.0 10.0 0.0 10.0 2                                                          
4 3 0.0 10.0 0.0 2.0 1                                                           
5 3 0.0 110.0 0.0 2.0 4
"""

with open('ballandstick2.swc', 'w') as swc_file:
    swc_file.write(morph_swc_string2)


cm = ephys.parameters.NrnSectionParameter(
        name="cm",
        param_name="cm",
        value=1.0, # in microfarad/cm2
        locations=[somatic_loc, dend_loc],
        frozen=False)

# Set boundaries.
gnabar_soma = ephys.parameters.NrnSectionParameter(                                    
        name="gnabar_soma",
        param_name="gnabar_hh",
        locations=[somatic_loc],
        bounds=[0.05, .125],
        frozen=False)     
gkbar_soma = ephys.parameters.NrnSectionParameter(
        name="gkbar_soma",
        param_name="gkbar_hh",
        bounds=[0.01, .05],
        locations=[somatic_loc],
        frozen=False)

gl_hh_soma = ephys.parameters.NrnSectionParameter(
        name="gl_hh_soma",
        param_name="gl_hh",
        bounds=[1e-4, 5e-4],
        locations=[somatic_loc],
        frozen=False)

hh_mech = ephys.mechanisms.NrnMODMechanism(
        name='hh',
        suffix='hh',
        locations=[somatic_loc])

# Compute the cell model.
ballandstick_cell = ephys.models.CellModel(
        name="simple_cell",
        morph=morph,
        mechs=[hh_mech],
        params=[cm, gnabar_soma, gkbar_soma, gl_hh])  

# Create protocols.
soma_loc = ephys.locations.NrnSeclistCompLocation(
        name="soma",
        seclist_name="somatic",
        sec_index=0,
        comp_x=0.5)

# Current clamp with square current inject
sweep_protocols = []
for protocol_name, amplitude in [("step1", 0.1), ("step2", 0.5), ("step3", -.3)]:
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

# Extract features.
efel_feature_means = {"step1": {"Spikecount": 1}, "step2": {"Spikecount": 5}, "step3": {"steady_state_voltage_stimend": "-100"}}

# Run protocol.
nrn = ephys.simulators.NrnSimulator()
default_params = {"gnabar_soma": 0.05, "gkbar_soma": 0.01}
responses = twostep_protocol.run(cell_model=ballandstick_cell, param_values=default_params, sim=nrn)
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
                    exp_std=0.05 * abs(float(mean)))
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
        
# Run the optimization algorithm.
optimization_algorithm = bpop.deapext.optimisations.IBEADEAPOptimisation(
        evaluator=cell_evaluator,
        offspring_size = 30)

# Run the optimization algorithm.
optimization_algorithm = bpop.deapext.optimisations.IBEADEAPOptimisation(
        evaluator=cell_evaluator,
        offspring_size = 10)

# An example answer of the expected output
best_ind_dict_ex = {'gl_hh': 0.0, 'gnabar_hh': 0.0, 'gkbar_hh': 0.0}
answer = '%f,%f,%f' % (best_ind_dict_ex['gnabar_hh'], best_ind_dict_ex['gkbar_hh'], best_ind_dict_ex['gl_hh'])
print('Answer: %s' % answer)
