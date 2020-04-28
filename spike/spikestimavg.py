# Spike-triggered stimulus average

spikes = mat["rho"]
stimVals = mat["stim"]
n_timeSteps = 150
spikeTimes = spikes[n_timeSteps:].nonzero()[0] + n_timeSteps
sepInterval = np.array([2,20,40,60,80,100])
sepInterval = (sepInterval/2).astype("int")
print("Time intervals:", sepInterval, "steps")
colors = ["b", "g", "r", "c", "m", "y"]
handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in colors]
labels = [str(int(x*2)) + " ms" for x in sepInterval]
MAX_SINGLE_RESPONSE = 29.472907029165032
diff_dict = dict()

j = 0
plt.figure(figsize=(15,5))
for interval in sepInterval:
    spikePairs = []
    i = 0
    while i < len(spikeTimes) - 1:
        if spikeTimes[i+1] - spikeTimes[i] == interval:
            spikePairs.append(spikeTimes[i+1]) 
        i += 1
    spikePairs = np.array(spikePairs)
    sta = np.zeros((n_timeSteps,))
    n_spikes = len(spikePairs)
    for sp_time in spikePairs:
        sta += stimVals[sp_time-n_timeSteps:sp_time]
    sta /= n_spikes
    diff_dict[interval*2] = abs(max(sta) - 2*MAX_SINGLE_RESPONSE)
    plt.plot(sta, c=colors[j], alpha=0.6)
    j += 1

plt.title("Spike-Triggered Average for Spike Pair")
plt.xlabel("Time Period (2ms)")
plt.ylabel("Stimulus Value")
plt.xticks(np.arange(0, 160, 10))
plt.legend(handles, labels)
plt.show()

plt.figure(figsize=(15,5))
plt.scatter(diff_dict.keys(), diff_dict.values())
plt.xlabel("Time Interval, ms")
plt.ylabel("Absolute diffrence")
plt.title("Absolute diffrence between two STA and the sum of two single STA by Time Interval")
plt.show()
