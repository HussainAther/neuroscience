function plotTuningCurveBar(meanResponse, SEM);
%
%  meanResponses is a list (vector) of mean responses
%  SEM is a list (vectore) of SEM associated with each of the means
%
%  

if (length(meanResponse) ~= length(SEM))
    error('Input vectors must have same length');
    return;
end
 

gray = [0.5 0.5 0.5];

x = 1:length(meanResponse);

figure; hold on;
h = bar(x,meanResponse,0.7);
set(h,'Facecolor',gray);
set(h,'Edgecolor',gray);

for i = 1:length(x);
    hh = plot([x(i) x(i)],[meanResponse(i)-SEM(i) meanResponse(i)+SEM(i)],'k-');
    set(hh,'linewidth',3);
end
xlabel('Condition');


