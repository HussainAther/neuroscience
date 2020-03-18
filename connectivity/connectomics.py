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

# Visualize the cell.
mm = trimesh_io.MeshMeta(disk_cache_path="test/test_files")
mesh = mm.mesh(filename ="/data/dynamic_brain_workshop/electron_microscopy/2019/meshes/%d.h5"%cellid)
mesh_poly =trimesh_vtk.trimesh_to_vtk(mesh.vertices,mesh.faces,None)
plt_actor = vtkplotter.Actor(mesh_poly)
vtkplotter.embedWindow(backend="k3d")
vp = vtkplotter.Plotter(bg="b")
myactor = vtkplotter.Actor(plt_actor, c="r")
myactor.GetMapper().Update()
vp+=myactor
vp.show()

# Find 10 largest synapses automatically extracted on this cell.
biggest_synapses = post_synapse_df.sort_values(by=["size"],ascending=False).head(10)
print(biggest_synapses)
