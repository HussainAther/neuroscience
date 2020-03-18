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

# Show in chart.
snr_df = session.units.sort_values(by=["snr"], ascending=False)
snr_df.head()

# Get spike times for 50 units with highest SNR.
unit_list = snr_df.index.values[:50]

# Figure setup
fig,ax = plt.subplots(5,10,figsize=(15,7),sharex=True)
ax = ax.ravel()

# Plot ISI distribution for each unit.
for i,unit in enumerate(unit_list):
    unit_spikes = session.spike_times[unit]
    if len(unit_spikes) > 3000:
        isi = np.diff(unit_spikes)
        ax[i].hist(isi,bins=50,range=(0,0.3))
        ax[i].set_title(str(unit_list[i]),fontsize=7)

plt.xlim(0,0.3)
for i in ax:
    i.axis("off")

# Explore responses to natural scenes.
stim_table = session.get_presentations_for_stimulus("natural_scenes")
stim_ids = stim_table.index.values
frames = session.get_stimulus_parameter_values(stimulus_presentation_ids=stim_ids, drop_nulls=False)["frame"]
frames = np.sort(frames)
frames

# Calculate spikes in 5 ms bins.
stim_presentation_ids = stim_table[stim_table.frame==46].index.values
bin_width = 0.005
duration = stim_table.duration.mean()
pre_time = -duration
post_time = 2*duration
bins = np.arange(pre_time, post_time+bin_width, bin_width)   
histograms = session.presentationwise_spike_counts(
    bin_edges=bins,
    stimulus_presentation_ids=stim_presentation_ids,
    unit_ids=None
)
mean_histograms = histograms.mean(dim="stimulus_presentation_id")
rates = mean_histograms/bin_width

# Plot the peristimulus time histogram (PSTH).
def plot_psth(unit_id, rates, ax=None, title=None):
    #Default params
    if not ax:
        fig,ax = plt.subplots(1,1,figsize=(6,3))
    rates.loc[{"unit_id":unit_id}].plot(ax=ax)
    ax.axvspan(0, duration, color="gray", alpha=0.1)
    ax.set_ylabel("Firing rate (spikes/second)")
    ax.set_xlabel("Time (s)")
    ax.set_xlim(pre_time,post_time)
    if ax:
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("")

unit_id = 914686471 #unit_list[0]
plot_psth(unit_id, rates)

unit_id = 914686471
fig,ax = plt.subplots(15, 8, figsize=(18,20), sharex=True, sharey=True)
ax = ax.ravel()

# Calculate histograms for all presentations at once - this may take a minute.
histograms = session.presentationwise_spike_counts(
    bin_edges=bins,
    stimulus_presentation_ids=stim_table.index.values,
    unit_ids=None
    )

for i,frame in enumerate(frames):
    stim_presentation_ids = stim_table[stim_table.frame==frame].index.values
    # Select the histograms for this frame and average.
    frame_histograms = histograms.loc[{"stimulus_presentation_id":stim_presentation_ids}]
    mean_histograms = frame_histograms.mean(dim="stimulus_presentation_id")
    rates = mean_histograms/bin_width
    
    plot_psth(unit_id, rates, ax=ax[i], title=frame)
plt.tight_layout()

# Plot the image that shows the largest response for this unit.
plt.imshow(cache.get_natural_scene_template(105), cmap="gray")

# Responses to natural scenes
spike_stats = session.conditionwise_spike_statistics(stimulus_presentation_ids=stim_ids, unit_ids=unit_list)
spike_stats.head()

unit_id = 914686471
response_mean = np.empty((len(frames)))
response_sem = np.empty((len(frames)))
for i,frame in enumerate(frames):
    stim_id = stim_table[stim_table.frame==frame].stimulus_condition_id.iloc[0]
    response_mean[i] = spike_stats.loc[(unit_id, stim_id)].spike_mean
    response_sem[i] = spike_stats.loc[(unit_id, stim_id)].spike_sem
plt.errorbar(range(118), response_mean[1:], yerr=response_sem[1:], fmt="o")
plt.axhspan(response_mean[0]+response_sem[0], response_mean[0]-response_sem[0], color="gray", alpha=0.3)
plt.axhline(y=response_mean[0], color="gray", ls="--")

# Responses to natural scenes
response_norm = response_mean - response_mean[0]
N = float(len(response_norm))
ls = ((1-(1/N) * ((np.power(response_norm.sum(),2)) / (np.power(response_norm,2).sum()))) / (1-(1/N)))
print(ls)
