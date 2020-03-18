import allensdk.brain_observatory.behavior.swdb.utilities as tools
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import platform
import seaborn as sns

from allensdk.brain_observatory.behavior.swdb import behavior_project_cache as bpc

sns.set_context("notebook", font_scale=1.5, rc={"lines.markeredgewidth": 2})
sns.set_style("white")
sns.set_palette("deep")

platstring = platform.platform()

if "Darwin" in platstring:
    # macOS
    data_root = "/Volumes/Brain2019/"
elif "Windows"  in platstring:
    # Windows (replace with the drive letter of the hard drive)
    data_root = "E:/"
elif ("amzn1" in platstring):
    # then on AWS
    data_root = "/data/"
else:
    # then linux (default here is for Ubuntu - insert your username; your distribution may differ)
    data_root = "/media/$USERNAME/Brain2019"
    
cache_path = os.path.join(data_root, "dynamic-brain-workshop/visual_behavior/2019")

"""
Visual behavior
"""

# Load data.
cache = bpc.BehaviorProjectCache(cache_path)

experiments = cache.experiment_table
selected_experiments = experiments[(experiments.full_genotype=="Vip-IRES-Cre/wt;Ai148(TIT2L-GC6f-ICL-tTA2)/wt") & (experiments.imaging_depth==175) & (experiments.stage_name=="OPHYS_4_images_B")]
selected_experiments
