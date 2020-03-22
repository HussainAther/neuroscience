function out = align_waveforms(x)
    % Determine alignments for a set of waves
    % relative to the initial waveform using
    % cross correlation.
    % Takes x as an MxN matrix of waveforms where x(m,:) is the 
    % nth waveform and outputs a vector of length M where out(m) is
    % the offset relative to the first wave.
