% Stimulus program for psychophysics experiment
temp = unit8(zeros(500,500,3)); % Create a dark stimulus matrix.
temp1 = cell(10,1); % Create a cell that can hold 10 matrices.
for ii = 1:10 
    temp(200,200,:) = 255; % Insert a fixation point.
    temp(200,240,:) = (ii-1)*10; % Insert a test point 40 pixels right of it
                                   % with brightness range 0 to 90.
    temp1(ii) = temp; % Put the respective modified matrix in cell.
end
h = figure; % Create a figure.
stimulusorder = randperm(200); % Create a random order from 1 to 200
                               % for the 200 trials. This lets you have 
                               % a precisely equal number per condition.
stimulusorder = mod(stimulusorder, 10); % Use the modulus function to create 
                                        % a range from 0 to 9. 20 each.
stimulusorder = stimulusorder + 1; % Now, the range is from 1 to 10.
score = zeros(10,1); % Keep score. How many stimuli were reported seen.
for ii = 1:200 % 200 trials and 20 per condition
    image(temp1{stimulusorder(1, ii)}); % Image the respective matrix.
    ii % Present the current trial.  
    pause;  
    temp2 = get(h, 'CurrentCharacter'); % Get the keypress. "." for present and "," absent
    temp3 = strcmp('.', temp2); % Compare strings. If present, temp3 = 1. Otherwise 0.
    score(stimulusorder(1, ii)) = score(stimulusorder(1, ii)) + temp3; % Add up.
end  
