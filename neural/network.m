%{
Groups of neurons acting together in networks.
}%
N_exc = 80; % excitatory neurons
N_inh = 20; % inhibitory neurons
v_i = 10; % initial velocity
time = 100; % time interval
neuronV = v_i + zeros(N_exc+N_inh,length(time));
