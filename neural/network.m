%{
Groups of neurons acting together in networks.
}%
N_exc = 80; % excitatory neurons
N_inh = 20; % inhibitory neurons
v_i = 10; % initial velocity
time = 100; % time interval
neuronV = v_i + zeros(N_exc+N_inh,length(time));
ti = 1; % initial time
volt_rest = -70; % resting potential (mV)
volt_thresh = -50; % action potential thresh. (mV)
volt_reset = -75; % post-spike reset voltage

% membrane parameters
R_m = 10; % neuron membrane resistance (MOhm)
tau = 10; % time constant of decay (ms)
volt_thresh = 15;
srate = 10000; % sampling rate in Hz
sim_dur = 1; % stimulus duration in seconds
time = 0:1/srate:sim_dur - 1/srate;
input = zeros(1,length(time)); input(dsearchn(time',.3):dsearchn(time',.7)) = 3;

if neuronV(ti) > volt_thresh
    r_i = volt_rest + input(ti)*R_m;
    neuronV(ti+1) = r_i + (neuronV(ti) - r_i) ...
                    * exp(-1000/srate/tau);
    neuronV(ti) = volt_reset;
    spiketimes = cat(1,spiketimes,ti);
end

spikedNeurons = neuronV(:,ti) > volt_thresh;
neuronV(spikedNeurons,ti) = volt_reset;
r_i = volt_rest + input(ti)*R_m;
n2stm = round(N_exc+N_inh)/2;
r_i = [r_i*ones(1,n2stm) volt_rest*ones(1,n2stm)]';

