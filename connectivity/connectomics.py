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
post_synapse_df = dl.query_synapses("pni_synapses_i3", post_ids = np.array([cellid]))
biggest_synapses = post_synapse_df.sort_values(by=["size"],ascending=False).head(10)
print(biggest_synapses)

# Visualize.
mm = trimesh_io.MeshMeta(disk_cache_path="test/test_files")
mesh = mm.mesh(filename ="/data/dynamic_brain_workshop/electron_microscopy/2019/meshes/%d.h5"%cellid)
mesh_poly =trimesh_vtk.trimesh_to_vtk(mesh.vertices,mesh.faces,None)
plt_actor = vtkplotter.Actor(mesh_poly)
syn_pts = np.vstack(biggest_synapses["ctr_pt_position"].values) * voxel_resolution
syn_sizes = biggest_synapses["size"]
syn_actors = trimesh_vtk.point_cloud_actor(syn_pts, size=syn_sizes.values)
vtkplotter.embedWindow(backend="k3d")
vp = vtkplotter.Plotter(bg="b")
myactor = vtkplotter.Actor(plt_actor, c="r")
myactor.GetMapper().Update()
mysynactor = vtkplotter.Actor(syn_actors, c="b")
mysynactor.GetMapper().Update()
vp+= mysynactor
vp.show()

# Find the axon initial segment (AIS) of a pyramidal cell.
mesh_folder = "/data/dynamic_brain_workshop/electron_microscopy/2019/meshes/"

# Query the AIS bounds table and find the AIS points with the pyramidal cell that
# have the functional ID of 1. 
nrn_ind = 1
manual_ais_df = dl.query_cell_ids("manual_ais")
nrn_id = manual_ais_df[manual_ais_df["func_id"]==nrn_ind]["pt_root_id"].values[0]
fname  = "{}.h5".format(nrn_id)

# Visualize.
print(mesh_folder +fname)
mm = trimesh_io.MeshMeta()
mesh = mm.mesh(filename=mesh_folder + "{}.h5".format(nrn_id))
mesh_poly =trimesh_vtk.trimesh_to_vtk(mesh.vertices,mesh.faces,None)
