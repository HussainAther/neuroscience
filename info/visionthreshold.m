% Stimulus program for psychophysics experiment
temp = unit8(zeros(500, 500, 3)); % Create a dark stimulus matrix.
temp1 = cell(10, 1); % Create a cell that can hold 10 matrices.
for ii = 1:10 
    temp(200, 200, :) = 255; % Insert a fixation point.
    temp(200, 240, :) = (ii-1)*10; % Insert a test point 40 pixels right of it
                                   % with brightness range 0 to 90.
    temp1(ii) = temp; % Put the respective modified matrix in cell.
