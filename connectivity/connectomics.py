import numpy as np
import pandas as pd
import scipy as sp
import vtkplotter

from matplotlib import cm
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
ais_bounds_df = dl.query_cell_ids("ais_bounds_v3")
nrn_bounds_df = ais_bounds_df[ais_bounds_df["pt_root_id"] == nrn_id]
voxel_resolution = np.array([4,4,40])
ais_bounds_pts = np.vstack(nrn_bounds_df["pt_position"].values) * voxel_resolution
pts_actor = trimesh_vtk.point_cloud_actor(ais_bounds_pts, size=1200, color=(0.2, 0.8, 0.2))
print(ais_bounds_pts)
vtkplotter.embedWindow(backend="k3d")
vp = vtkplotter.Plotter(bg="b")
mesh_poly_actor = vtkplotter.Actor(mesh_poly)
mesh_poly_actor.GetMapper().Update()
vp+=mesh_poly_actor
plot_actor1 = vtkplotter.Actor(pts_actor,c="r")
plot_actor1.GetMapper().Update()
vp+=plot_actor1
vp.show()

# Find the vertex indices of the mesh closest to the AIS bounds points.
ds, mesh_inds = mesh.kdtree.query(ais_bounds_pts)
dist_from_pts = np.linalg.norm(mesh.vertices[mesh_inds] - ais_bounds_pts, axis=1)
print(dist_from_pts)
mesh_pt_actor = trimesh_vtk.point_cloud_actor(mesh.vertices[mesh_inds],
                                              size=1200, color=(0.9, 0, 0.8))

vtkplotter.embedWindow(backend="k3d")
vp = vtkplotter.Plotter(bg="b")
plot_actor = vtkplotter.Actor(mesh_pt_actor,c="b")
plot_actor.GetMapper().Update()
vp+=plot_actor
plot_actor1 = vtkplotter.Actor(pts_actor,c="r")
plot_actor1.GetMapper().Update()
vp+=plot_actor1
vp.show()

# Find the vertex indices of the mesh "between" the top and bottom parts.
ds = sp.sparse.csgraph.dijkstra(mesh.csgraph, indices=mesh_inds)
d_padding = 2000
ais_len = ds[0,mesh_inds[1]]
is_ais = ds.sum(axis=0) < ais_len + d_padding

# Visualize the AIS in the whole neuron mesh.
ais_mesh = mesh.apply_mask(is_ais)
rest_of_mesh = mesh.apply_mask(~is_ais)
ama_poly =trimesh_vtk.trimesh_to_vtk(ais_mesh.vertices,ais_mesh.faces,None)
rma_poly =trimesh_vtk.trimesh_to_vtk(rest_of_mesh.vertices,rest_of_mesh.faces,None)
vtkplotter.embedWindow(backend="k3d")
vp = vtkplotter.Plotter(bg="b")
ama_poly_actor = vtkplotter.Actor(ama_poly,c="r")
ama_poly_actor.GetMapper().Update()
rma_poly_actor = vtkplotter.Actor(rma_poly,c="b")
rma_poly_actor.GetMapper().Update()
vp+=ama_poly_actor
vp+=rma_poly_actor
vp.show()

# Get an array of synapse locations onto the same cell you looked at above.
# Filter only those onto mesh points corresponding to the AIS.
syn_df = dl.query_synapses("pni_synapses_i3", post_ids=[int(nrn_id)])
syn_locs = np.vstack(syn_df["ctr_pt_position"].values) * voxel_resolution
_, syn_mesh_inds = mesh.kdtree.query(syn_locs)
ais_syn_df = syn_df[is_ais[syn_mesh_inds]]
ais_syn_df

# Count how many synapses per presynaptic object there are.
pre_ids, cts = np.unique(ais_syn_df["pre_pt_root_id"], return_counts=True)
count_df = pd.DataFrame(data={"root_id":pre_ids, "synapses":cts})
count_df

# Visualize.
ais_pts = np.vstack(ais_syn_df["ctr_pt_position"].values) * voxel_resolution
ais_pre_ids = ais_syn_df["pre_pt_root_id"].values
c_mapping = {}
for ii, oid in enumerate(pre_ids):
    c_mapping[oid] = cm.tab20.colors[ii]

s_mapping = {}
for ii, oid in enumerate(pre_ids):
    s_mapping[oid]= cts[ii]
    
pt_colors = [c_mapping[oid] for oid in ais_pre_ids]
pt_sizes = [100 * s_mapping[oid] for oid in ais_pre_ids]
syn_actors = trimesh_vtk.point_cloud_actor(ais_pts, color=pt_colors, size=pt_sizes)
vtkplotter.embedWindow(backend="k3d")
vp = vtkplotter.Plotter(bg="b")
plot_actor = vtkplotter.Actor(syn_actors,c="b")
plot_actor.GetMapper().Update()
vp+=plot_actor
vp.show()

# Plot.
fig, ax = plt.subplots(1, 1, figsize=(3, 3))
ax.hist(inh_dist, density=True, alpha=0.5, label='Inhibitory')
ax.hist(ex_dist, density=True, alpha=0.5, label='Excitatory')
ax.set_xlabel('Num. inhibitory input cells')
ax.set_ylabel('Frequency')
ax.legend(loc=0, frameon=False)

sns.despine(fig)
fig.tight_layout()
