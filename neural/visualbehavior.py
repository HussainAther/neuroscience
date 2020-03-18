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

# Print dataset dimensions.
print("targeted structures:", experiments.targeted_structure.unique())
print("\ncre_lines:", experiments.full_genotype.unique())
print("\nstage_types:", experiments.stage_name.unique())

# Random active behavior experiment
active_experiments = experiments[experiments.passive_session==False]
experiment_id = active_experiments.ophys_experiment_id.sample(1).values[0]

# Metadata
row = experiments[experiments["ophys_experiment_id"] == experiment_id]
print(row.targeted_structure.values[0])
print(row.imaging_depth.values[0])
print(row.full_genotype.values[0])
print(row.stage_name.values[0])

# Container ID
container_id = experiments[experiments.ophys_experiment_id==experiment_id]["container_id"].values[0]
experiments.groupby("container_id").get_group(container_id)
