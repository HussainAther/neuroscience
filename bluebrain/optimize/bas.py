import bluepyopt as bpop
import bluepyopt.ephys as ephys
import matplotlib.pyplot as plt
import neurom
import neurom.view

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
