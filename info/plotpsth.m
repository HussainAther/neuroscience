function psth,bins= bin_for_psth(raw_dta, sampling_rate, threshold, trial_count, trial_length, bin_size)
    % Locate evens above threshold in raw data and generate PSTH (peristimulus time histogram)
    % from multi-trial recording. Trials should be contiguous.
