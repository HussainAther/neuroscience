import numpy as np
import pandas as pd
import vtkplotter

from meshparty import trimesh_io, trimesh_vtk
from analysisdatalink.datalink_ext import AnalysisDataLinkExt as AnalysisDataLink

"""
Functionality for electron microscopy dataset
"""

dataset_name = "pinky100"
data_version = 175
sqlalchemy_database_uri = ""
dl = AnalysisDataLink(dataset_name=dataset_name,
                      sqlalchemy_database_uri=sqlalchemy_database_uri,
                      materialization_version=data_version,
                      verbose=False)

mesh_folder = "/data/dynamic_brain_workshop/electron_microscopy/2019/"
voxel_resolution = np.array([4,4,40])

# Find all synapses.
all_man_syn = dl.query_synapses("gsynapse_ai_manual_v2")

# Find the postsynaptic cell ID with the most synapses that are manually annotated.
all_man_syn["syn_num"]=all_man_syn.groupby("post_pt_root_id")["id"].transform(len)
cellid = all_man_syn[all_man_syn.syn_num==34]["post_pt_root_id"].values[0]
print(cellid)
