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
