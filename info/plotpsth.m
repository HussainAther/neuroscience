function psth,bins= bin_for_psth(raw_dta, sampling_rate, threshold, trial_count, trial_length, bin_size)
    % Locate evens above threshold in raw data and generate PSTH (peristimulus time histogram)
    % from multi-trial recording. Trials should be contiguous.
    % For input data raw_data, sampling rate in Hz, sampling_rtae, threshold for events
    % threshold, number of contiguous trials trial_count, length of each trial trial_length,
    % and size of each PSTH bin bin_size, output the psth (count in each bin) and bins 
    % (center position of ech bin, in seconds relative to trial start). All time units 
    % are in seconds.
    events = raw_data > threshold;
    % only positive threshold crossings (not sustained activity above threshold)
    events = diff(events) == 1;
    events = [0 events];
    events = reshape(events, trial_count, trial_length*sampling_rate);
    % events should be MxN for M = trial and N = sample
    summed_events = sum(events);
    % summmed_events should be the sum of events at each sample relative to
    % the start of the trial
    max_event_count = max(summed_events);
    for count = 0:max_event_count-1
        above_count = find(summed_events > count);
        event_offsets = [event_offsets above_count];
    end
    % event_offsets should be the offset in sample counts where events occur
    event_times = event_offsets / sampling_rate;
    [psth, bins] = hist(event_times, trial_length/bin_size);
end
