% Use a composite model.
% Load spiking data of 15 neurons.
load data/example15

% Randomly divide it into a training set and a test set (so we can verify how well we trained).
[ncells,nsamples] = size(spikes15);
idx_train = randperm(nsamples,ceil(nsamples/2));
idx_test = setdiff(1:nsamples,idx_train);
samples_train = spikes15(:,idx_train);
samples_test = spikes15(:,idx_test);

% Create a model with independent factors, k-synchrony and third-order correlations.
% We will do this by initializing 3 separate models and then combining them to a single model.
third_order_correlations = num2cell(nchoosek(1:ncells,3),2);
model_indep = maxent.createModel(ncells,'indep');
model_ksync = maxent.createModel(ncells,'ksync');
model_thirdorder = maxent.createModel(ncells,'highorder',third_order_correlations);
model = maxent.createModel(ncells,'composite',{model_indep,model_ksync,model_thirdorder});

% Train it.
model = maxent.trainModel(model,samples_train,'threshold',1);

% Use the model to predict the frequency of activity patterns.
% We will start by observing all the patterns that repeated at least twice (because a pattern that repeated at least
% once may grossly overrepresent its probability and is not meaningful in this sort of analysis)
limited_empirical_distribution = maxent.getEmpiricalModel(samples_test,'min_count',2);

% Get the model predictions for these patterns.
model_logprobs = maxent.getLogProbability(model,limited_empirical_distribution.words);

% Plot on a log scale.
% Output saved to "doc/composite.png".
figure
plot(limited_empirical_distribution.logprobs,model_logprobs,'bo');
hold on;
minval = min(limited_empirical_distribution.logprobs);
plot([minval 0],[minval 0],'-r');  % identity line
xlabel('empirical pattern log frequency');
ylabel('predicted pattern log frequency');
title(sprintf('Composite model: activity patterns in %d cells',ncells));
