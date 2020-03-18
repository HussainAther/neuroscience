import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import vtk
import vtkplotter as vtkp

from analysisdatalink import datalink_ext as de
from analysisdatalink.datalink_ext import AnalysisDataLinkExt as AnalysisDataLink
from annotationframeworkclient import infoservice
from cloudvolume import CloudVolume
from matplotlib import colors
from meshparty import skeletonize, trimesh_io
from meshparty.trimesh_vtk import trimesh_to_vtk
from meshparty import mesh_filters, trimesh_vtk
from pykdtree.kdtree import KDTree
from scipy.spatial import cKDTree
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

"""
Inhibitory cell connectomics
"""

# Initialize.
dataset_name = "pinky100"
data_version = 175
sqlalchemy_database_uri = ""
dl = AnalysisDataLink(dataset_name=dataset_name,
                      sqlalchemy_database_uri=sqlalchemy_database_uri,
                      materialization_version=data_version,
                      verbose=False)
mesh_folder = "/data/dynamic_brain_workshop/electron_microscopy/2019/meshes/"
voxel_resolution = np.array([4,4,40])

# Get tables from the dataset.
print(dl.sqlalchemy_engine.table_names())

# Get soma valence v2 information.
cell_types_df = dl.query_cell_types("soma_valence_v2")
cell_types_df.head()

# Get cell type.
print(cell_types_df["cell_type"].unique())

# How many inhibitory cells and excitatory cells are labeled?
inh_ids = cell_types_df.loc[cell_types_df["cell_type"] == "i"]["pt_root_id"]
exc_ids = cell_types_df.loc[cell_types_df["cell_type"] == "e"]["pt_root_id"]
all_ids = pd.concat((inh_ids, exc_ids))
print(str(len(inh_ids))+" inhibitory")
print(str(len(exc_ids))+" excitatory")
print(str(len(all_ids))+" neuron ids")

inh_ids = inh_ids.unique()
exc_ids = exc_ids.unique()
all_ids = all_ids.unique()
print(str(len(inh_ids))+" unique inhibitory ids")
print(str(len(exc_ids))+" unique excitatory ids")
print(str(len(all_ids))+" unique neuron ids")

# Get synapse dataframe and synapse information.
synapse_df = dl.query_synapses("pni_synapses_i3", pre_ids=inh_ids, post_ids=all_ids)
synapse_df.head()

# Make an inhibitory-to-all connectivity matrix using the total synapse size between
# two neurons as the connection strength.
Ni = len(inh_ids)
Ne = len(exc_ids)
J = np.zeros((Ne+Ni, Ni)) # post, pre
for j, pre in enumerate(inh_ids):
    this_pre = synapse_df.loc[synapse_df["pre_pt_root_id"] == pre]
    for i, post in enumerate(all_ids):
        this_pre_post = this_pre.loc[this_pre["post_pt_root_id"] == post]
        J[i, j] = this_pre_post["size"].sum()

# View the matrix.
fig, ax = plt.subplots(1, 1, figsize=(3, 3))
ax.imshow(J, cmap="gray_r", aspect=0.1)
ax.set_xlabel("Presynaptic")
ax.set_ylabel("Postsynaptic")

# Compute frequency distributions of connection number from excitatory neurons
# to all inhibitory postsynaptic neurons (e-i connections).
A = np.minimum(J, np.ones(J.shape))
ex_dist = np.sum(A[Ni:, :Ni], axis=1)
