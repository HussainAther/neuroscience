% Trains a independent maximum-entropy model on a set of samples. This is a sub-function 
% invoked by trainModel() and typically should not be invoked directly.
%
% Usage:
%   model_out = trainModelExhaustive(model,samples)
%
% Input:
%   refer to trainModel.m
%
% Output:
%   model_out - trained ME model (normalized and with the entropy stored in model_out.H)
%
% Last update: Ori Maoz 20/11/2016
function model_out = trainModelIndependent(model,raster)


if (~strcmp(model.type,'indep'))
   error('This function can only train independent models'); 
end


% check if our input was a raster (of samples from a distribution) or an actual Boltzmann distribution
if (isstruct(raster))
    % we got a model (that we can compute the marginals directly on)
    use_exact_marginals = true;
    model_base = raster;
    ncells = model_base.ncells;
else
    use_exact_marginals = false;
    [ncells,nsamples] = size(raster);    
    if (nsamples < ncells)    
        warning('Input raster must be of the form (nsamples x ncells), are you sure you did not transpose it?');
    end
end



if (ncells ~= model.ncells)
    error('Number of cells in model does not match number of cells in data');
    
end


% compute the observables which are the mean firing rates. If we are preseted with experimental data,
% just take the empirical count. Otherwise we need to compute it explicitly on the base model
if (use_exact_marginals)
    npatterns = 2^ncells;
    unique_words = logical(de2bi(0:(npatterns-1)))';
    
    % get probabilities of the input patterns
    lp = maxent.getLogProbability(model_base,unique_words);
    lp = lp - log(sum(exp(lp(:)))); % normalize
    empirical_probs = exp(lp);
    
    % get the observables - our goal is to fit the model marginals to these.
    firing_rates = maxent.getWeightedMarginals(unique_words,model,empirical_probs)';

else    
    % mean firing rates
    firing_rates = mean(raster,2);
end


model_out = model;

% we need special handling for cells that never fire or that always fire
idx_always_zero = logical(firing_rates == 0);
idx_always_one = logical(firing_rates == 1);

% assume that our measurements are approximate: if we encountered a marginal of zero, then the actual marginal
% is just below the quantization error i.e. 1/(2*nsamples). Similarly a marginal of one is actually just below one.
firing_rates(idx_always_zero) = 1 / (2*nsamples);
firing_rates(idx_always_one) = 1 - (1 / (2*nsamples));

% compute the partition function: z = (1-p1)(1-p2).....(1-pn)
logz = sum(log(1-firing_rates));


% now compute the entropy
model_out.entropy = sum(-firing_rates .* log2(firing_rates) - (1-firing_rates) .* log2(1-firing_rates));


% compute each of the factors
for i = 1:ncells
    
    % first compute (1-p1)(1-p2)....(pi)(1-p_1+1)...(1-pn)
    vec = ones(ncells,1);
    vec(i) = 2*firing_rates(i);
    lambda = -(sum(log(vec - firing_rates)) - logz);
    model_out.factors(i) = lambda;
end


model_out.z = logz;


end