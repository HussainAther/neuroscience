function [SEM] = computeSEM(responses);

n = length(responses);  %% number of observations

if (n<=1) error('Need at least two responses to compute SEM.'); return; end;

standardDeviation = std(responses);
SEM = standardDeviation/sqrt(n);

