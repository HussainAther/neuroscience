function [spikeTimesMS] = spikeExtractorForVoltageMatrix(timesMS, voltageMatrixUV);
%
%   You need to edit one line of this function to make it work with your
%   own spike detector.
%   You may also want to check that your detector is working on this data.
%
[nRuns nSamples] = size(voltageMatrixUV);

if (nSamples~= length(timesMS)) 
    error('The number of columns of the voltage matrix is not equal to the number of time values');
    return;
end


for run = 1:nRuns
    voltageUV = voltageMatrixUV(run,:);     %% get one row of the voltage matrix
    
    % now you have data just like in project 1:  an list of times and a
    % list of voltages
    
    %% you should uncomment the next line and replace with your function
    %% (and then comment out the text line)
    
    %spikeTimesOnOneRunMS = runMyProject1SpikeDetector(timesMS,voltageUV);    
    spikeTimesOnOneRunMS = spikedetector(timesMS,voltageUV);
    %['YOU NEED TO CALL YOUR OWN FUNCTION IN HERE.']
    
    spikeTimesMS{run} = spikeTimesOnOneRunMS;   %% a list of spike times for this run put in one cell of a cell array

end


   
    
