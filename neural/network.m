%{
Groups of neurons acting together in networks.
}%
N_exc = 80; % excitatory neurons
N_inh = 20; % inhibitory neurons
v_i = 10; % initial velocity
time = 100; % time interval
neuronV = v_i + zeros(N_exc+N_inh,length(time));
spikedNeurons = neuronV(:,ti) > volt_thresh;
neuronV(spikedNeurons,ti) = volt_reset;
r_i = volt_rest + input(ti)*R_m;
n2stm = round(N_exc+N_inh)/2;
r_i = [r_i*ones(1,n2stm) volt_rest*ones(1,n2stm)]';

