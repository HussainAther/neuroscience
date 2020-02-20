%{
Any neurobiological model should account for the perceived successes of the diffusion
model in describing the characteristics of reaction time in decision processes under psychometric 
testing. One such model was proposed by Shadlen and Newsome (2001) to account
for interactions between visual areas MT and LIP. Area MT is a cortical area sensitive to
visual motion, and area LIP is a cortical area implicated in decision processes. Under this
model, neurons in area LIP function as integrators, accumulating rate information from
neurons in area MT over time. Neurons in area MT, whose preferred directions oppose each
other, inhibit the corresponding LIP neuron of the opposing neuron. In other words, LIP
neurons integrate the difference between sensory cells with opposing preferred direction
sensitivities. For this simulation, you will draw firing rates from a normal distribution of typical rates
for an MT neuron, varying with orientation.
%}

function rate = mtneurons(preferred, stimulus, neurons)
    % assuming index 1 is preferred
    typical_rate_mean = [30 20 15 10 5 5 5 5 5 5 5 5 5 10 15 20];
    typical_rate_stdev = sqrt(typical_rate_mean);
    mean = typical_rate_mean(1+mod(stimulus-preferred));
    stdev = typical_rate_stdev(1+mod(stimulus-preferred));
    rate = normpdf(mean*ones(1,neurons), stdev*ones(1,neurons));
end
