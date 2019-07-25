% Leaky Integrate-and-fire (IAF) neuron
clear; clf;

% parameters
tau_inv=.1; % inverse time constant
uu=0; % initial membrane voltage
I_ext=12; % constant external input
tspan=[0 100] % integration interval
theta=10; % firing threshold
