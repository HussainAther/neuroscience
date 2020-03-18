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
sqlalchemy_database_uri = "postgres://postgres:welcometothematrix@swdb-em-db.crjvviai1xxh.us-west-2.rds.amazonaws.com"
dl = AnalysisDataLink(dataset_name=dataset_name,
                      sqlalchemy_database_uri=sqlalchemy_database_uri,
                      materialization_version=data_version,
                      verbose=False)
