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
