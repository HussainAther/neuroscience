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
sessions.head()

# Recordings across visual areas
cache = EcephysProjectCache.fixed(manifest=manifest_path)
sessions = cache.get_sessions()
vis_areas = ["VISp","VISl","VISal","VISrl","VISam","VISpm"]
vis_session_list = []
for session_id in sessions.index:
    session_areas = sessions.structure_acronyms.loc[session_id]
    vis_areas_in_session = [area for area in vis_areas if area in session_areas]
    if len(vis_areas_in_session)==6:
        vis_session_list.append(session_id)
print(vis_session_list)
session_id = vis_session_list[0]
session = cache.get_session_data(session_id)

# Visualize.
for i, area in enumerate(vis_areas):
    vis_unit_list = session.units[session.units.structure_acronym==area].index
    for i,unit in enumerate(vis_unit_list[:20]):
        spike_times =  session.spike_times[unit]                       
        plt.plot(spike_times, np.repeat(i, len(spike_times)), "|", color="gray") 
    plt.title(area) 
    plt.xlim(0,300)
    plt.show()

# Inter-spike interval distribution (ISI distribution)
unit_id = session.units.index.values[0]
unit_spikes = session.spike_times[unit_id]
isi = np.diff(unit_spikes)
isi = isi*1000 # convert to ms
print(isi[:20])

# Visualize.
fig,ax = plt.subplots(1,2,figsize=(10,4),sharey=True)
ax[0].hist(isi,bins=200,range=(0,200))
ax[1].hist(isi,bins=20,range=(0,20))
plt.ylabel("Count")
plt.xlabel("Inter-spike interval (ms)")
plt.show()
