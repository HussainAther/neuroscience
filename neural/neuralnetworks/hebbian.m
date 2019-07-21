function net = hebbRNN_create_model(N, B, I, p, g, dt, tau, varargin)

% net = hebbRNN_create_model(N, B, I, p, g, dt, tau, varargin)
%
% This function initializes a recurrent neural network for later training
% and execution
%
% INPUTS:
%
% N -- the number of recurrent neurons in network
%
% B -- the number of outputs
%
% I -- the number of inputs
%
% p -- the sparseness of the J (connectivity) matrix, (range: 0-1)
%
% g -- the spectral scaling of J
%
% dt - the integration time constant
%
% tau - the time constant of each neuron
%
%
% OPTIONAL INPUTS:
%
% actFun -- the activation function used to tranform activations into
% firing rates
% Default: 'tanh'
%
% netNoiseSigma - the variance of random gaussian noise added at each time
% point
% Default: 0
%
% useBiasNeurons -- whether or not to included neurons that have a fixed
% value
% Default: false
%
% numBiasNeurons -- how many bias neurons to use
% Default: 0
%
% feedback -- whether or not to feed the output of the network back
% Default: false
%
% energyCost -- how much to weight the firing rate of the output units in
% the calculation of total error. This parameter encourages the network to
% find solutions with as low firing rate as possible.
% Default: 0
%
% outputUnitIdent -- user-defined indexing of output units
% Default: 1:B
%
% biasUnitIdent -- user-defined indexing of bias units
% Default: []
%
%
% OUTPUTS:
%
% net -- the network structure

actFunType = 'tanh'; % Default activation function
netNoiseSigma = 0.0; % Default noise-level
useBiasNeurons = false; % Default use of bias neurons
numBiasNeurons = 0; % Default number of bias neurons
feedback = false; % Default use of output feedback
energyCost = 0; % Default weight of cost function
outputUnitIdent = 1:B; % Default identity of output units
biasUnitIdent = []; % Default identity of bias units
