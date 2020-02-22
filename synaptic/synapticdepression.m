% Synaptic depression

function [tv, Pv, Vv, Qv, ind_v] = mt_mod7(ap_freq, t_max)
%
% ap_freq: stimulation frequency (Hz)
% t_max: end of simulation (in ms)
% 
%
% returns the time vector corresponding to the simulations (tv)
% the resource vector Pv, with maximum value of 1
% the membrane potential vector Vv (in mV)
% the charge vector Qv 
% a vector of indices of the time of action potentials (ind_v)
%
%basic setup for the Markram-Tsodyks model
%This file generates a stochastic MT model
%

t_rec = 800; %ms
dt = 0.5; %ms

t_frac = dt/t_rec;
den_t = (1 + t_frac);

Pvmax = 1;
Pv0 = Pvmax;
p = 0.55;

if ( isempty(t_max) )
    t_max = 1000;
end;

tv = 0:dt:t_max;
n_tv = length(tv);

Pv = zeros(size(tv));
Pv(1) = Pv0;

if ( isempty(ap_freq) )
    ap_freq = 20; % 1/s
end;

%convert to time interval in ms and round to the closest bin
dt_ap = round((1000/ap_freq)/dt);

%time of first pulse 
t_aps = 20; %in ms
t_ind = round(t_aps/dt);
%indices of pulses
ind_v = t_ind:dt_ap:n_tv;

%leaky IF membrane potential
Vv = zeros(size(tv));
V0 = 0; %rest is at zero
Vv(1) = V0;

tm = 50; %ms
tm_frac = dt/tm;
den_v = (1 + tm_frac);

%synaptic factor
Rin = 250; %MOhm
qmax = 1.5; %pC
syn = p*Rin*qmax/(Pvmax*tm);
ap = (p*qmax/Pvmax); %in pC
Qv = zeros(size(tv));

for i = 2:n_tv
    Pv(i) = (Pv(i-1) + t_frac * Pvmax)/den_t;

    %backward Euler on membrane potential
    Vv(i) = Vv(i-1)/den_v;
    
    if ( ~isempty(find(ind_v == i)) )
        %time for a pulse
        
        %first update V and Q by taking the fraction p of the 
        %available resources
        Vv(i) = Vv(i) + syn*Pv(i)/den_v;
        Qv(i) = ap*Pv(i);
        
        %then update the available resources as 1-p those immediately
        %before the AP
        Pv(i) = (1-p)*Pv(i);

    end;
end;


h_f1 = figure; 
h_a1 = subplot(4,1,1);
h_a2 = subplot(4,1,2);
h_a3 = subplot(4,1,3);
h_a4 = subplot(4,1,4);

plot_intermediate = 1;

%stimulation frequencies
stim_freq = [1 2:2:40];
nstim_freq = length(stim_freq);

freq_q = zeros(1,nstim_freq);
freq_v = zeros(1,nstim_freq);
freq_m = zeros(1,nstim_freq);

for j = 1:nstim_freq
%for j = 6:6
    [tv, Pv, Vv, Qv, ind_v] = mt_mod7(stim_freq(j),5000);

    ind_max = length(tv);

    ind_ss = find(ind_v > 2000);
    nind_ss = length(ind_ss);

    Q_vals = zeros(1,nind_ss);
    V_vals = zeros(1,nind_ss);

    for i = 1:nind_ss-1
        Q_vals(i) = max(Qv(ind_v(ind_ss(i)):ind_v(ind_ss(i+1))));
        V_vals(i) = max(Vv(ind_v(ind_ss(i)):ind_v(ind_ss(i+1))));
    end;

    Q_vals(nind_ss) = max(Qv(ind_v(ind_ss(nind_ss)):ind_max));
    V_vals(nind_ss) = max(Vv(ind_v(ind_ss(nind_ss)):ind_max));

    freq_q(j) = mean(Q_vals);
    freq_v(j) = mean(V_vals);

    freq_m(j) = mean(Vv(ind_ss(1):ind_max));

    if ( (j == 6) & (plot_intermediate == 1) )
        line('Parent',h_a1,'XDAta',tv,'YData',Pv);
        set(h_a1,'XLim',[0 1000]);
        
        for i = 1:length(ind_v)
            line('Parent',h_a2,'XData',[tv(ind_v(i)) tv(ind_v(i))],'YData',[0 1]);
        end;
        set(h_a3,'XLim',[0 1000]);
        
        line('Parent',h_a3,'XDAta',tv,'YData',Vv);
        set(h_a3,'XLim',[0 1000]);
        ylabel(h_a3,'Vm (mV)');
        
        line('Parent',h_a4,'XDAta',tv,'YData',Qv);
        set(h_a4,'XLim',[0 1000]);
        xlabel(h_a4,'time (ms)');
        ylabel(h_a4,'charge (pC)')
    end;
end;

%

h_f2 = figure; 
h_a5 = subplot(3,1,1);
h_a6 = subplot(3,1,2);
h_a7 = subplot(3,1,3);

line('Parent',h_a5,'XData',stim_freq,'YData',freq_q);
ylabel(h_a5,'charge (pC)');
line('Parent',h_a6,'XData',stim_freq,'YData',freq_v);
ylabel(h_a6,'peak Vm (mV)');
line('Parent',h_a7,'XData',stim_freq,'YData',freq_m);
xlabel(h_a7,'stimulation frequency (Hz)');
ylabel(h_a7,'average Vm (mV)');

%theoretical curves
qmax = 1.5; %pC
p = 0.55;
tau_rec = 800; %ms
tau_rec_sec = tau_rec *1e-3; %in seconds

v_e = exp(-1./(stim_freq*tau_rec_sec));

q_ss = (qmax*p)*(1-v_e)./(1-(1-p)*v_e);

line('Parent',h_a5,'XData',stim_freq,'YData',q_ss,...
    'LineStyle','none','Marker','o','MarkerFaceColor','none',...
    'MarkerEdgeColor','r','MarkerSize',4);

Rin = 250; %MOhm
tm = 50; %ms
tm_sec = tm*1e-3; %in seconds

v_e2 = 1./(1-exp(-1./(stim_freq*tm_sec)));

v_ss = (q_ss*Rin/tm).*v_e2;

line('Parent',h_a6,'XData',stim_freq,'YData',v_ss,...
    'LineStyle','none','Marker','o','MarkerFaceColor','none',...
    'MarkerEdgeColor','r','MarkerSize',4);

v_av = q_ss*Rin .* stim_freq*1e-3;
line('Parent',h_a7,'XData',stim_freq,'YData',v_av,...
    'LineStyle','none','Marker','o','MarkerFaceColor','none',...
    'MarkerEdgeColor','r','MarkerSize',4);

lim_val = qmax*Rin/tau_rec;

line('Parent',h_a6,'XData',[0 40],'YData',[lim_val lim_val],...
    'LineStyle','--');

line('Parent',h_a7,'XData',[0 40],'YData',[lim_val lim_val],...
    'LineStyle','--');

%print(h_f1,'-depsc2','mt_mod6.eps');


