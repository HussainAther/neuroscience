import numpy as np

"""
We model Lobula plate tangential cells (LPTC lptc) used in the fly visual ganglia. The lobula plate is
in the posterior part of the lobula complex (the third visual ganglion of the fly eye). We use electrophysiological
recordings to study the large diameter of the intracellular processes. The integrating neurons can be modeled on a 
biophysical level with direct recordings of their membrane properties. 

We measure the response feature of LPTCs with the spatial integration characteristics when enlarging the area
in which the motion stimulus is displayed. The response saturates signficantly not only for motion along the preferred
direction, but also along the nulll direction. The gain control gives us different saturation plateaus.
"""

def V(E, g, gl):
    """
    Given a list of excitatory and inhibitory potentials E (tuple array) and excitatory and inhibitory conductances g 
    (tuple array) and a leak conductance gl, calcualte the membrane potential V.
    """
