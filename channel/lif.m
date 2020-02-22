%leaky integrate and fire with random inputs  (LIF)

%test first constant current injection

%time constant
tau = 20; %ms

%input resistance
R = 10; %MOhms

v_thres = 10; %mV

%integration time step
dt = 0.05; %ms

%simulation time
t_max = 1000;
t_v = 0:dt:t_max-dt;
nt_v = length(t_v);

%compute constants to save time
dttau = dt/tau;
one_dttau = 1-dttau; %forward Euler
m1_dttau = 1/(1+dttau); %backward Euler

%number of simulations
%n_sims = 500;
n_sims = 1; %debugging

%input current
%constant current for debugging purposes
%i_v = 2*ones(size(t_v)); %nA

%to test the effect of a single synaptic input
%i_v = (qex/dt)*[zeros(size([0:dt:10-dt])) 1 zeros(size([10+dt:dt:t_max-dt]))];

%
%
% compute response to current-based excitation only.
%
%

%charge and number of excitatory inputs
qex = 2; %pCb
nex = 500;

%save number of spikes per trial and isi distribution
n_spks_e = zeros(1,n_sims);
isis_e = [];

for j = 1:n_sims
    %EPSP input times uniformly distributed between 0 and t_max
    r_v = sort(t_max*rand(1,nex));

    %find corresponding indices in t_v
    ind_v = ceil(r_v/dt);
    ind_v = ind_v ( (ind_v > 0) & (ind_v <= nt_v) );
    
    %generate pulse current vector
    i_v = zeros(size(t_v));
    i_v(ind_v) = 1;
    i_v = (qex/dt)*i_v;
    i_v_n = dttau * R * i_v;

    v_v = zeros(size(t_v));
    s_v = zeros(size(t_v));

    %simulate LIF response to input
    for i = 2:nt_v
        %forward Euler
        %v_v(i) = one_dttau * v_v(i-1) + i_v_n(i-1);
    
        %backward Euler
        v_v(i) = m1_dttau*(v_v(i-1) + i_v_n(i));
    
        if (v_v(i) >= v_thres)
            v_v(i) = 0;
            s_v(i) = 1;
        end;
    end;
    
    %spike indices
    ind_spk = find(s_v > 0.5);
    
    %compute number of spikes
    n_spks_e(j) = length(ind_spk);
    
    isis_e = [isis_e diff(t_v(ind_spk))]; 
end;

h_f1 = figure;
h_a1 = axes;

line('Parent',h_a1,'XData',t_v,'YData',v_v);


for i = 1:length(ind_spk);
    line('Parent',h_a1,'XData',[t_v(ind_spk(i)) t_v(ind_spk(i))],'YData',[0 30]);
end;

set(h_a1,'YLim',[-5 30]);

%for constant current input, compare simulated firing rate with theoretical f-i curve
%of LIF neuron
%fr_sim = sum(s_v)/t_max; %spk/ms
%fr_th = -1/(tau*log(1-v_thres/(R*i_v(1))));

%
%
% current-based excitation and inhibition
%
%

%number of inhibitory inputs adjusted for similar firing rate as without
%inhibition
qex = 2; %pCb
nex = 660;
qin = 4; %pCb
nin = 100;

%save number of spikes per trial and isi distribution
n_spks_i = zeros(1,n_sims);
isis_i = [];

for j = 1:n_sims
    %EPSP input times uniformly distributed between 0 and t_max
    r_v = sort(t_max*rand(1,nex));

    %find corresponding indices in t_v
    ind_v = ceil(r_v/dt);
    ind_v = ind_v ( (ind_v > 0) & (ind_v <= nt_v) );
    
    %generate pulse current vector
    i_v = zeros(size(t_v));
    i_v(ind_v) = 1;
    i_v = (qex/dt)*i_v;
    i_v_n = dttau * R * i_v;

    %IPSP input times uniformly distributed between 0 and t_max
    r_v = sort(t_max*rand(1,nin));

    %find corresponding indices in t_v
    ind_v = ceil(r_v/dt);
    ind_v = ind_v ( (ind_v > 0) & (ind_v <= nt_v) );
    
    %generate pulse current vector
    i_v = zeros(size(t_v));
    i_v(ind_v) = 1;
    i_v = (qin/dt)*i_v;
    
    %add to current vector
    i_v_n = i_v_n - dttau * R * i_v;
    
    v_v = zeros(size(t_v));
    s_v = zeros(size(t_v));


    %simulate LIF response to input
    for i = 2:nt_v
        %forward Euler
        %v_v(i) = one_dttau * v_v(i-1) + i_v_n(i-1);
    
        %backward Euler
        v_v(i) = m1_dttau*(v_v(i-1) + i_v_n(i));
    
        if (v_v(i) >= v_thres)
            v_v(i) = 0;
            s_v(i) = 1;
        end;
    end;
    
    %spike indices
    ind_spk = find(s_v > 0.5);
    
    %compute number of spikes
    n_spks_i(j) = length(ind_spk);    
    isis_i = [isis_i diff(t_v(ind_spk))]; 
end;

h_f2 = figure;
h_a2 = axes;

line('Parent',h_a2,'XData',t_v,'YData',v_v);

for i = 1:length(ind_spk);
    line('Parent',h_a2,'XData',[t_v(ind_spk(i)) t_v(ind_spk(i))],'YData',[0 30]);
end;

set(h_a2,'YLim',[-5 30]);

%
%
% conductance-based excitation
%
%

%number of inhibitory inputs adjusted for similar firing rate as without
%inhibition
qex = 2; %pCb
nex = 600;

%excitatory synaptic conductance
t_ex = 1; %ms
v_ex = 70; %mV

%the excitatory conductance in normalized to yield the same charge at rest 
%as the excitatory current sysnapse
g_ex = qex/(exp(1)*v_ex); %microS

%alpha function for excitatory input
t_syn = 0:dt:8*t_ex-dt;
g_se = (g_ex/t_ex)*t_syn.*exp(1-t_syn/t_ex);

%save number of spikes per trial and isi distribution
n_spks_eg = zeros(1,n_sims);
isis_eg = [];

for j = 1:n_sims
    %EPSP input times uniformly distributed between 0 and t_max
    r_v = sort(t_max*rand(1,nex));

    %find corresponding indices in t_v
    ind_v = ceil(r_v/dt);
    ind_v = ind_v ( (ind_v > 0) & (ind_v <= nt_v) );
    
    %generate pulse synaptic activation vector
    i_v = zeros(size(t_v));
    i_v(ind_v) = 1;
    %generate conductance vector
    g_v = conv(i_v,g_se);
    g_v = g_v(1:nt_v);
    %execute multiplication in advance to save time
    g_v_n = dttau * R * g_v;
    
    v_v = zeros(size(t_v));
    s_v = zeros(size(t_v));


    %simulate LIF response to input
    for i = 2:nt_v
        %forward Euler
        %v_v(i) = one_dttau * v_v(i-1) + g_v_n(i-1)*(v_ex-v_v(i-1));
    
        %backward Euler
        v_v(i) = (v_v(i-1) + g_v_n(i)*v_ex)/(1 + dttau + g_v_n(i));
    
        if (v_v(i) >= v_thres)
            v_v(i) = 0;
            s_v(i) = 1;
        end;
    end;
    
    %spike indices
    ind_spk = find(s_v > 0.5);
    
    %compute number of spikes
    n_spks_eg(j) = length(ind_spk);    
    isis_eg = [isis_eg diff(t_v(ind_spk))]; 
end;

h_f3 = figure;
h_a3 = axes;
line('Parent',h_a3,'XData',t_v,'YData',v_v);

for i = 1:length(ind_spk);
    line('Parent',h_a3,'XData',[t_v(ind_spk(i)) t_v(ind_spk(i))],'YData',[0 30]);
end;

set(h_a3,'YLim',[-5 30]);

%
%
% conductance-based excitation and inhibition
%
%

%number of inhibitory inputs adjusted for similar firing rate as without
%inhibition
qex = 2; %pCb
nex = 690;
qin = 4; %pCb
nin = 100;

%excitatory synaptic conductance
t_ex = 1; %ms
v_ex = 70; %mV

%the excitatory conductance in normalized to yield the same charge at rest 
%as the excitatory current sysnapse
g_ex = qex/(exp(1)*v_ex); %microS

%alpha function for excitatory input
t_syn_e = 0:dt:8*t_ex-dt;
g_se = (g_ex/t_ex)*t_syn_e.*exp(1-t_syn_e/t_ex);

%excitatory synaptic conductance
t_in = 4; %ms
v_in = 0; %mV

%the inhibitory conductance in normalized 
%as the excitatory synapse
g_in = qin/(exp(1)*v_ex); %microS

%alpha function for inhibitory input
t_syn_i = 0:dt:8*t_in-dt;
g_si = (g_in/t_in)*t_syn_i.*exp(1-t_syn_i/t_in);

%save number of spikes per trial and isi distribution
n_spks_ig = zeros(1,n_sims);
isis_ig = [];

for j = 1:n_sims
    %EPSP input times uniformly distributed between 0 and t_max
    r_v = sort(t_max*rand(1,nex));

    %find corresponding indices in t_v
    ind_v = ceil(r_v/dt);
    ind_v = ind_v ( (ind_v > 0) & (ind_v <= nt_v) );
    
    %generate pulse synaptic activation vector
    i_v = zeros(size(t_v));
    i_v(ind_v) = 1;
    %generate conductance vector
    g_ve = conv(i_v,g_se);
    g_ve = g_ve(1:nt_v);
    %execute multiplication in advance to save time
    g_ve_n = dttau * R * g_ve;

    %IPSP input times uniformly distributed between 0 and t_max
    r_v = sort(t_max*rand(1,nin));

    %find corresponding indices in t_v
    ind_v = ceil(r_v/dt);
    ind_v = ind_v ( (ind_v > 0) & (ind_v <= nt_v) );
    
    %generate pulse synaptic activation vector
    i_v = zeros(size(t_v));
    i_v(ind_v) = 1;
    %generate conductance vector
    g_vi = conv(i_v,g_si);
    g_vi = g_vi(1:nt_v);
    %execute multiplication in advance to save time
    g_vi_n = dttau * R * g_vi;
        
    v_v = zeros(size(t_v));
    s_v = zeros(size(t_v));


    %simulate LIF response to input
    for i = 2:nt_v
        %forward Euler
        %v_v(i) = one_dttau * v_v(i-1) + g_ve_n(i-1)*(v_ex-v_v(i-1)) + g_vi_n(i-1)*(v_in-v_v(i-1));
    
        %backward Euler
        v_v(i) = (v_v(i-1) + g_ve_n(i)*v_ex + g_vi_n(i)*v_in)/(1 + dttau + g_ve_n(i) + g_vi_n(i));
    
        if (v_v(i) >= v_thres)
            v_v(i) = 0;
            s_v(i) = 1;
        end;
    end;
    
    %spike indices
    ind_spk = find(s_v > 0.5);
    
    %compute number of spikes
    n_spks_ig(j) = length(ind_spk);    
    isis_ig = [isis_ig diff(t_v(ind_spk))]; 
end;

h_f4 = figure;
h_a4 = axes;
line('Parent',h_a4,'XData',t_v,'YData',v_v);

for i = 1:length(ind_spk);
    line('Parent',h_a4,'XData',[t_v(ind_spk(i)) t_v(ind_spk(i))],'YData',[0 30]);
end;

xlabel(h_a4,'time (ms)');
ylabel(h_a4,'membrane potential (mV)');

set(h_a4,'YLim',[-5 30]);

%print(handles.figure1,'-depsc','lif_rand_inp.eps');
