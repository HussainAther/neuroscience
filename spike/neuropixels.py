import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
from allensdk.brain_observatory.ecephys import ecephys_session

cache = EcephysProjectCache.fixed(manifest=manifest_path)

"""
Neuropixels dataset from the Allen Institute
"""

# Load data.
import platform
platstring = platform.platform()

if "Darwin" in platstring:
    # OS X 
    data_root = "/Volumes/Brain2019/"
elif "Windows"  in platstring:
    # Windows (replace with the drive letter of USB drive)
    data_root = "E:/"
elif ("amzn1" in platstring):
    # then on AWS
    data_root = "/data/"
else:
    # then your own linux platform
    # EDIT location where you mounted hard drive
    data_root = "/media/$USERNAME/Brain2019/"

manifest_path = os.path.join(data_root, "dynamic-brain-workshop/visual_coding_neuropixels/2019/manifest.json")

# Get information.
sessions = cache.get_sessions()
