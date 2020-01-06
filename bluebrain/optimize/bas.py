import bluepyopt as bpop
import bluepyopt.ephys as ephys
import matplotlib.pyplot as plt
import neurom
import neurom.view

"""
Optimize parameters for a ball-and-stick model (ball and stick)
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

fig, ax = neurom.viewer.draw(neurom.load_neuron("ballandstick.swc"))

morph = ephys.morphologies.NrnFileMorphology("ballandstick.swc")
