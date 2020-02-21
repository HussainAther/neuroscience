% Impact of spontaneous synaptic activity on membrane potential. 
% set up the destexhe equations on a patch
% and solve via hybrid euler. Includes background
% synaptic activity.

function [n_mean s_mean m_isis c_isis] = dest_f2

%compute f-i curve of destexhy2 model
%as well as cv 

dt = 0.05;

currents = [ 4.5:0.5:9 ];
n_currents = length(currents);

n_trials = 10;
n_spk = zeros(n_currents,n_trials);
m_isis = zeros(1,n_currents);
s_isis = zeros(1,n_currents);
n_isis = zeros(1,n_currents);

for j = 1:n_currents
    %info_str = sprintf('processing current %i ...',j);
    %disp(info_str);
    isis = [];
    
    for k = 1:n_trials
        v = destexhy2(dt,10,1010,1020,currents(j));

        %find spike times
        spk = [];
        inds = find(v > 0);
        if ( ~isempty(inds) )
            %at least one spike
            gaps = diff(inds);
            gi = find(gaps>1);
            if ( ~isempty(gi) )
                %at least two spikes
                [valm indm] = max(v(inds(1:gi(1))));
                spk(1) = indm + (inds(1)-1);

                for i =2:length(gi)
                    %if there are 3 or more spikes
                    [valm indm] = max(v(inds(gi(i-1)+1:gi(i))));
                    spk(i) = indm + (inds(gi(i-1)+1)-1);
                end;
        
                %last spike
                [valm indm] = max(v(inds(gi(end)+1:end)));
                spk(end+1) = indm + (inds(gi(end)+1)-1);
            else
                %get the only spike
                [valm indm] = max(v(inds(1:end)));
                spk(1) = indm + (inds(1) - 1);
            end;
        end;
        isis = [isis diff(spk)];
        n_spk(j,k) = length(spk);
    end;
    m_isis(j) = mean(isis);
    s_isis(j) = std(isis);
    n_isis(j) = length(isis);
end;

n_mean = mean(n_spk,2);
s_mean = std(n_spk,0,2);
c_isis = s_isis./m_isis;

%figure; 
%h = axes; 
%line('Parent',h,'XData',currents,'YData',n_mean,'Marker','o','MarkerFaceColor','k');
%for i = 1:n_currents
%    line('Parent',gca,'XData',[currents(i) currents(i)],...
%        'YData',[n_mean(i)-s_mean(i) n_mean(i)+s_mean(i)]); 
%end;
%set(h,'XLim',[currents(1)-1 currents(n_currents)+1]);

%figure; 
%plot(m_isis,s_isis./m_isis,'o');

function v = destexhy(dt,t1,t2,T,Istim)

% maximal conductances, in mS / (cm)^2
%gK = 0;
%gM = 0;
%gNa = 0;
gK = 10; %delayed rectifier            
gM = 0.5; %M current
gNa = 51.6; %HH-type sodium channel        
gl = 0.045; %leak           

% reversal potentials, in mV		
VK = -90;             
VNa = 50;           
Vl = -80;         

%capacitance
C = 1;                  % micro F / (cm)^2

N = floor((T-dt)/dt);
v = zeros(N,1);		% allocate space for v
%t = (1:N)'*dt;

% initial conditions
v(1) = -80.3935;

%steady-state values
m = am(v(1))/(am(v(1)) + bm(v(1)));
h = ah(v(1))/(ah(v(1)) + bh(v(1)));
n = an(v(1))/(an(v(1)) + bn(v(1)));
p = ap(v(1))/(ap(v(1)) + bp(v(1)));

j = 2;			% time counter

while j*dt < T,		% run for T ms

t = j*dt;

    Iapp = 0;		% apply Istim between t1 and t2 
    if (t>t1 & t<t2)
       Iapp = Istim;
    end

    % advance the gating variables, m, h, n, and p using hybrid Euler

    tmp = am(v(j-1));
    m = (m+dt*tmp)/(1+dt*(tmp+bm(v(j-1))));

    tmp = ah(v(j-1));
    h = (h+dt*tmp)/(1+dt*(tmp+bh(v(j-1))));

    tmp = an(v(j-1));
    n = (n+dt*tmp)/(1+dt*(tmp+bn(v(j-1))));

    tmp = ap(v(j-1));
    p = (p+dt*tmp)/(1+dt*(tmp+bp(v(j-1))));

    % compute the conductances of the K, Na, NaP, h & Ca currents

    cIK = gK*(n^4);
    cINa = gNa*m^3*h;
    
    %M current 
    cIM = gM*p;

    % advance the voltage, v

    tmp1 = cINa*VNa+(cIK+cIM)*VK+gl*Vl+Iapp;
    tmp2 = cINa+cIK+cIM+gl;

    %the correct formula is 
    %v(j) = (C*v(j-1) + dt*tmp1)/(C + dt*tmp2)
    %but with C = 1 the following will do.
    v(j) = (v(j-1) + dt*tmp1)/(1 + dt*tmp2);
   
    j = j + 1;

end

return

% steady-state activation, time constant and opening and closing
% functionals from Destexhe and  Pare


%forward rate of delayed rectifier
function val = an(V)
VT = -58;
tmp = (V-VT-15);
val = -0.032*tmp./(exp(-tmp/5)-1);

%backward rate of delayed rectifier
function val = bn(V)
VT = -58;
val = 0.5*exp(-(V-VT-10)/40);

%forward rate of fast sodium current  activation
function val = am(V)
VT = -58; 
tmp = (V-VT-13);
val = -0.32*tmp/(exp(-tmp/4)-1);

%backward rate of fast sodium current activation
function val = bm(V)
VT = -58;
tmp = (V-VT-40);
val = 0.28*tmp/(exp(tmp/5)-1);

%forward rate of fast sodium current inactivation
function val = ah(V)
VT = -58;
VS = -10;
tmp = (V-VT-VS-17);
val = 0.128*exp(-tmp/18);

%backward rate of fast sodium current inactivation
function val = bh(V)
VT = -58;
VS = -10;
tmp = (V-VT-VS-40);
val = 4/(1+exp(-tmp/5));

%forward rate of non-inactivating K current
function val = ap(V)
tmp = (V+30);
val = 0.0001*tmp./(1-exp(-tmp/9));

%backward rate of non-inactivating K current
function val = bp(V)
tmp = (V+30);
val = -0.0001*tmp./(1-exp(tmp/9));

%compute the distribution of the membrane potential
dt = 0.05;
v_stat = destexhy2(dt,10,100,2000,0);

%skip the first 50 ms
n1 = ceil(50/dt);
mv = mean(v_stat(n1:end));
sv = std(v_stat(n1:end));

text_str = sprintf('mean vm: %.3g std: %.3g',mv,sv);
disp(text_str);

function v  = destexhy2(dt,t1,t2,T,Istim)

% maximal conductances, in mS / (cm)^2
%gK = 0;
%gM = 0;
%gNa = 0;
gK = 10; %delayed rectifier            
gM = 0.5; %M current
gNa = 51.6; %HH-type sodium channel        
gl = 0.045; %leak           

% reversal potentials, in mV		
VK = -90;             
VNa = 50;           
Vl = -80;         

%capacitance
C = 1;                  % micro F / (cm)^2

N = floor((T-dt)/dt);
v = zeros(N,1);		% allocate space for v
t = (1:N)'*dt;

%total dendritic area
Ad = 34636e-8; %cm2

%--------------------

%mean excitatory conductance
ge0 = 0.012e-3; %mS 
Ve = 0;

%relaxation time constant and steady-state SD
tau_e = 2.7; %ms
sig_e = 0.003e-3; %mS

%numerical simulation factors
ef_1 = exp(-dt/tau_e);
ef_2 = sqrt(1-exp(-2*dt/tau_e));

%initial value of the random conductance variable
x0 = 0;

%random white noise increment
w_ve = normrnd(0,1,1,N);
u_ve = sig_e*ef_2*w_ve;

%save space for random vector
x_ve = zeros(1,N);

%initial condition
x_ve(1) = x0;

for i =2:N
    %stochastic update
    x_ve(i) = x_ve(i-1)*ef_1 + u_ve(i);
end;

%add steady-state value
x_ve = x_ve + ge0;

cAMPA = x_ve/Ad; %mS/cm2

%-------------------

%mean inhibitory conductance
gi0 = 0.057e-3; %mS 
Vi = -80;

%relaxation time constant and steady-state SD
tau_i = 10.5; %ms
sig_i = 0.0066e-3; %mS

%numerical simulation factors
if_1 = exp(-dt/tau_i);
if_2 = sqrt(1-exp(-2*dt/tau_i));

%initial value of the random conductance variable
x0 = 0;

%random white noise increment
w_vi = normrnd(0,1,1,N);
u_vi = sig_i*if_2*w_vi;

%save space for random vector
x_vi = zeros(1,N);

%initial condition
x_vi(1) = x0;

for i =2:N
    %stochastic update
    x_vi(i) = x_vi(i-1)*if_1 + u_vi(i);
end;

%add steady-state value
x_vi = x_vi + gi0;
cGABA = x_vi/Ad; %mS/cm2

%-------------------

% initial conditions
%v(1) = -80.3935;
v(1) = -69;

%steady-state values
m = am(v(1))/(am(v(1)) + bm(v(1)));
h = ah(v(1))/(ah(v(1)) + bh(v(1)));
n = an(v(1))/(an(v(1)) + bn(v(1)));
p = ap(v(1))/(ap(v(1)) + bp(v(1)));

j = 2;			% time counter

%while j*dt < T,		% run for T ms

for j = 2:N
    tc = t(j);
    
    Iapp = 0;		% apply Istim between t1 and t2 
    if (tc>t1 & tc<t2)
       Iapp = Istim;
    end

    % advance the gating variables, m, h, n, and p using hybrid Euler

    tmp = am(v(j-1));
    m = (m+dt*tmp)/(1+dt*(tmp+bm(v(j-1))));

    tmp = ah(v(j-1));
    h = (h+dt*tmp)/(1+dt*(tmp+bh(v(j-1))));

    tmp = an(v(j-1));
    n = (n+dt*tmp)/(1+dt*(tmp+bn(v(j-1))));

    tmp = ap(v(j-1));
    p = (p+dt*tmp)/(1+dt*(tmp+bp(v(j-1))));

    % compute the conductances of the K, Na, NaP, h & Ca currents

    cIK = gK*(n^4);
    cINa = gNa*m^3*h;
    
    %M current 
    cIM = gM*p;

    % advance the voltage, v

    tmp1 = cINa*VNa+(cIK+cIM)*VK+gl*Vl+cAMPA(j)*Ve+cGABA(j)*Vi+Iapp;
    tmp2 = cINa+cIK+cIM+gl+cAMPA(j)+cGABA(j);

    %the correct formula is 
    %v(j) = (C*v(j-1) + dt*tmp1)/(C + dt*tmp2)
    %but with C = 1 the following will do.
    v(j) = (v(j-1) + dt*tmp1)/(1 + dt*tmp2);
   
end

return

% steady-state activation, time constant and opening and closing
% functionals from Destexhe and  Pare


%forward rate of delayed rectifier
function val = an(V)
VT = -58;
tmp = (V-VT-15);
val = -0.032*tmp./(exp(-tmp/5)-1);

%backward rate of delayed rectifier
function val = bn(V)
VT = -58;
val = 0.5*exp(-(V-VT-10)/40);

%forward rate of fast sodium current  activation
function val = am(V)
VT = -58; 
tmp = (V-VT-13);
val = -0.32*tmp/(exp(-tmp/4)-1);

%backward rate of fast sodium current activation
function val = bm(V)
VT = -58;
tmp = (V-VT-40);
val = 0.28*tmp/(exp(tmp/5)-1);

%forward rate of fast sodium current inactivation
function val = ah(V)
VT = -58;
VS = -10;
tmp = (V-VT-VS-17);
val = 0.128*exp(-tmp/18);

%backward rate of fast sodium current inactivation
function val = bh(V)
VT = -58;
VS = -10;
tmp = (V-VT-VS-40);
val = 4/(1+exp(-tmp/5));

%forward rate of non-inactivating K current
function val = ap(V)
tmp = (V+30);
val = 0.0001*tmp./(1-exp(-tmp/9));

%backward rate of non-inactivating K current
function val = bp(V)
tmp = (V+30);
val = -0.0001*tmp./(1-exp(tmp/9));

%histogram of values
[nv, xoutv] = hist(v_stat(n1:end),50);

h_f1 = figure; 
h_a2 = subplot(2,2,2);
bar(h_a2,xoutv,nv,'r');
set(h_a2,'TickDir','out','XLim',[-80 -60]);
xlabel(h_a2,'Vm (mV)');
ylabel(h_a2,'number of occurrences');

%plot sample noise trace
t = (1:length(v_stat))*dt;
inds = find( (t>500)&(t<1000));
h_a1 = subplot(2,2,1);
line('Parent',h_a1,'XData',t(inds),'YData',v_stat(inds));
set(h_a1,'XLim',[500 1000]);

vn = destexhy(dt,150,350,550,-0.25);
N = length(vn);
t = (1:N)*dt;
h_a3 = subplot(2,2,3);
line('Parent',h_a3,'XData',t,'YData',vn);
xlabel(h_a3,'time (ms)');
ylabel(h_a3,'Vm (mv)');

reps = 1000;
vb = zeros(N,reps);

disp('iterating 1000 times. This may take a while...');
for i = 1:reps
    if ( mod(i,100) == 0 )
        info_str = sprintf('at iteration %i',i);
        disp(info_str);
    end;
    
    vb(:,i) = destexhy2(dt,150,350,550,-0.25);
end;
vba = mean(vb,2);

h_a4 = subplot(2,2,4);
line('Parent',h_a4,'XData',t,'YData',vba);
xlabel(h_a4,'time (ms)');

inds1 = find( (t>40) & (t<140) );
inds2 = find( (t>200) & (t<300) );

vpre = mean(vba(inds1));
vpul = mean(vba(inds2));
vdiff = vpre-vpul;

text_str2 = sprintf('mean vm prepulse: %.2g during: %.2g difference: %.2g',vpre,vpul,vdiff);
disp(text_str2);

vpren = mean(vn(inds1));
vpuln = mean(vn(inds2));
vdiffn = vpren-vpuln;

text_str3 = sprintf('mean vm prepulse: %.2g during: %.2g difference: %.2g',vpren,vpuln,vdiffn);
disp(text_str3);

%print(handles.figure1,'-depsc2','destex_f1.eps');

% Suprathreshold effects of spontaneous synaptic activity.

dt = 0.05;
v = destexhy2(dt,10,1010,1020,5);
N = length(v);
t = (1:N)*dt;

h_f1 = figure; 
h_a1 = subplot(2,2,1);
line('Parent',h_a1,'XData',t,'YData',v)
set(h_a1,'XLim',[200 800],'YLim',[-80 50]);
xlabel(h_a1,'time (ms)');
ylabel(h_a1,'Vm (mV)');

disp('computing 1st panel...');
[n_mean s_mean m_isis c_isis] = dest_f2;
h_a2 = subplot(2,2,2);
line('Parent',h_a2,'XData',m_isis,'YData',c_isis,'Marker','o',...
'LineStyle','none','MarkerFaceColor','k','MarkerSize',4);
set(h_a2,'XLim',[0 1500],'YLim', [0 1.1]);
xlabel(h_a2,'mean ISI (ms)');
ylabel(h_a2,'CV');


disp('computing 2nd panel...');
[currents n_mean n_mean2 n_mean3] = dest_f2b;
h_a3 = subplot(2,2,3);
line('Parent',h_a3,'XData',currents,'YData',n_mean,...
    'Marker','o','MarkerFaceColor','k','MarkerSize',4);
line('Parent',h_a3,'XData',currents,'YData',n_mean2,...
    'Marker','o','MarkerFaceColor','r','MarkerSize',4);
line('Parent',h_a3,'XData',currents,'YData',n_mean3,...
    'Marker','o','MarkerFaceColor','w','MarkerSize',4);
set(h_a3,'YLim',[0 300]);
xlabel(h_a3,'current (muA/cm2)');
ylabel(h_a3,'firing frequency (spk/s)');

disp('computing 3rd panel...');
[currents n_mean n_mean2 n_mean3] = dest_f2c;
h_a4 = subplot(2,2,4);
line('Parent',h_a4,'XData',currents,'YData',n_mean,...
    'Marker','o','MarkerFaceColor','k','MarkerSize',4);
line('Parent',h_a4,'XData',currents,'YData',n_mean2,...
    'Marker','o','MarkerFaceColor','r','MarkerSize',4);
line('Parent',h_a4,'XData',currents,'YData',n_mean3,...
    'Marker','o','MarkerFaceColor','w','MarkerSize',4);
set(h_a4,'YLim',[0 300]);
xlabel(h_a4,'current (muA/cm2)');

%print(handles.figure1,'-depsc2','destex_f2.eps');

function [currents n_mean n_mean2 n_mean3] = dest_f2b
%compute f-i curve of destexhy3 model
%changes peak conductance of inhibitory 
%conductance
dt = 0.05;

currents = [ 4:2:28 ];
n_currents = length(currents);

n_trials = 10;
n_spk = zeros(n_currents,n_trials);

for j = 1:n_currents
    %info_str = sprintf('processing current %i ...',j);
    %disp(info_str);
    isis = [];
    
    for k = 1:n_trials
        v = destexhy3(dt,10,1010,1020,currents(j),[]);

        %find spike times
        spk = [];
        inds = find(v > 0);
        if ( ~isempty(inds) )
            %at least one spike
            gaps = diff(inds);
            gi = find(gaps>1);
            if ( ~isempty(gi) )
                %at least two spikes
                [valm indm] = max(v(inds(1:gi(1))));
                spk(1) = indm + (inds(1)-1);

                for i =2:length(gi)
                    %if there are 3 or more spikes
                    [valm indm] = max(v(inds(gi(i-1)+1:gi(i))));
                    spk(i) = indm + (inds(gi(i-1)+1)-1);
                end;
        
                %last spike
                [valm indm] = max(v(inds(gi(end)+1:end)));
                spk(end+1) = indm + (inds(gi(end)+1)-1);
            else
                %get the only spike
                [valm indm] = max(v(inds(1:end)));
                spk(1) = indm + (inds(1) - 1);
            end;
        end;
        n_spk(j,k) = length(spk);
    end;
end;

n_mean = mean(n_spk,2);
s_mean = std(n_spk,0,2);

% figure; 
% h = axes; 
% line('Parent',h,'XData',currents,'YData',n_mean,'Marker','o','MarkerFaceColor','k');
% for i = 1:n_currents
%     line('Parent',gca,'XData',[currents(i) currents(i)],...
%         'YData',[n_mean(i)-s_mean(i) n_mean(i)+s_mean(i)]); 
% end;
% set(h,'XLim',[currents(1)-1 currents(n_currents)+1]);


n_spk2 = zeros(n_currents,n_trials);
syn_par = [0.012e-3 0.003e-3 2*0.057e-3 0.0066e-3];

for j = 1:n_currents
    %info_str = sprintf('processing current %i ...',j);
    %disp(info_str);
    isis = [];
    
    for k = 1:n_trials
        v = destexhy3(dt,10,1010,1020,currents(j),syn_par);

        %find spike times
        spk = [];
        inds = find(v > 0);
        if ( ~isempty(inds) )
            %at least one spike
            gaps = diff(inds);
            gi = find(gaps>1);
            if ( ~isempty(gi) )
                %at least two spikes
                [valm indm] = max(v(inds(1:gi(1))));
                spk(1) = indm + (inds(1)-1);

                for i =2:length(gi)
                    %if there are 3 or more spikes
                    [valm indm] = max(v(inds(gi(i-1)+1:gi(i))));
                    spk(i) = indm + (inds(gi(i-1)+1)-1);
                end;
        
                %last spike
                [valm indm] = max(v(inds(gi(end)+1:end)));
                spk(end+1) = indm + (inds(gi(end)+1)-1);
            else
                %get the only spike
                [valm indm] = max(v(inds(1:end)));
                spk(1) = indm + (inds(1) - 1);
            end;
        end;
        n_spk2(j,k) = length(spk);
    end;
end;

n_mean2 = mean(n_spk2,2);
s_mean2 = std(n_spk2,0,2);

% line('Parent',h,'XData',currents,'YData',n_mean2,'Marker','o','MarkerFaceColor','r');
% for i = 1:n_currents
%     line('Parent',gca,'XData',[currents(i) currents(i)],...
%         'YData',[n_mean2(i)-s_mean2(i) n_mean2(i)+s_mean2(i)],'Color','r'); 
% end;

n_spk3 = zeros(n_currents,n_trials);
syn_par = [0.012e-3 0.003e-3 0.5*0.057e-3 0.0066e-3];

for j = 1:n_currents
    %info_str = sprintf('processing current %i ...',j);
    %disp(info_str);
    isis = [];
    
    for k = 1:n_trials
        v = destexhy3(dt,10,1010,1020,currents(j),syn_par);

        %find spike times
        spk = [];
        inds = find(v > 0);
        if ( ~isempty(inds) )
            %at least one spike
            gaps = diff(inds);
            gi = find(gaps>1);
            if ( ~isempty(gi) )
                %at least two spikes
                [valm indm] = max(v(inds(1:gi(1))));
                spk(1) = indm + (inds(1)-1);

                for i =2:length(gi)
                    %if there are 3 or more spikes
                    [valm indm] = max(v(inds(gi(i-1)+1:gi(i))));
                    spk(i) = indm + (inds(gi(i-1)+1)-1);
                end;
        
                %last spike
                [valm indm] = max(v(inds(gi(end)+1:end)));
                spk(end+1) = indm + (inds(gi(end)+1)-1);
            else
                %get the only spike
                [valm indm] = max(v(inds(1:end)));
                spk(1) = indm + (inds(1) - 1);
            end;
        end;
        n_spk3(j,k) = length(spk);
    end;
end;

n_mean3 = mean(n_spk3,2);
s_mean3 = std(n_spk3,0,2);

% line('Parent',h,'XData',currents,'YData',n_mean3,'Marker','o','MarkerFaceColor','b');
% for i = 1:n_currents
%     line('Parent',gca,'XData',[currents(i) currents(i)],...
%         'YData',[n_mean3(i)-s_mean3(i) n_mean3(i)+s_mean3(i)],'Color','b'); 
% end;
function [currents n_mean n_mean2 n_mean3] = dest_f2c
%compute f-i curve of destexhy3 model
%changes inhibitory SD
dt = 0.05;

currents = [ 4:2:28 ];
n_currents = length(currents);

syn_par_def = [0.012e-3 0.003e-3 0.057e-3 0.0066e-3];

n_trials = 5;
n_spk = zeros(n_currents,n_trials);

for j = 1:n_currents
    %info_str = sprintf('processing current %i ...',j);
    %disp(info_str);
    isis = [];
    
    for k = 1:n_trials
        v = destexhy3(dt,10,1010,1020,currents(j),[]);

        %find spike times
        spk = [];
        inds = find(v > 0);
        if ( ~isempty(inds) )
            %at least one spike
            gaps = diff(inds);
            gi = find(gaps>1);
            if ( ~isempty(gi) )
                %at least two spikes
                [valm indm] = max(v(inds(1:gi(1))));
                spk(1) = indm + (inds(1)-1);

                for i =2:length(gi)
                    %if there are 3 or more spikes
                    [valm indm] = max(v(inds(gi(i-1)+1:gi(i))));
                    spk(i) = indm + (inds(gi(i-1)+1)-1);
                end;
        
                %last spike
                [valm indm] = max(v(inds(gi(end)+1:end)));
                spk(end+1) = indm + (inds(gi(end)+1)-1);
            else
                %get the only spike
                [valm indm] = max(v(inds(1:end)));
                spk(1) = indm + (inds(1) - 1);
            end;
        end;
        n_spk(j,k) = length(spk);
    end;
end;

n_mean = mean(n_spk,2);
s_mean = std(n_spk,0,2);

% figure; 
% h = axes; 
% line('Parent',h,'XData',currents,'YData',n_mean,'Marker','o','MarkerFaceColor','k');
% for i = 1:n_currents
%     line('Parent',gca,'XData',[currents(i) currents(i)],...
%         'YData',[n_mean(i)-s_mean(i) n_mean(i)+s_mean(i)]); 
% end;
% set(h,'XLim',[currents(1)-1 currents(n_currents)+1]);


n_spk2 = zeros(n_currents,n_trials);
syn_par(1) = 4*syn_par_def(1);
syn_par(3) = 4*syn_par_def(3);
syn_par(2) = 8*syn_par_def(2);
syn_par(4) = 8*syn_par_def(4);

for j = 1:n_currents
    %info_str = sprintf('processing current %i ...',j);
    %disp(info_str);
    isis = [];
    
    for k = 1:n_trials
        v = destexhy3(dt,10,1010,1020,currents(j),syn_par);

        %find spike times
        spk = [];
        inds = find(v > 0);
        if ( ~isempty(inds) )
            %at least one spike
            gaps = diff(inds);
            gi = find(gaps>1);
            if ( ~isempty(gi) )
                %at least two spikes
                [valm indm] = max(v(inds(1:gi(1))));
                spk(1) = indm + (inds(1)-1);

                for i =2:length(gi)
                    %if there are 3 or more spikes
                    [valm indm] = max(v(inds(gi(i-1)+1:gi(i))));
                    spk(i) = indm + (inds(gi(i-1)+1)-1);
                end;
        
                %last spike
                [valm indm] = max(v(inds(gi(end)+1:end)));
                spk(end+1) = indm + (inds(gi(end)+1)-1);
            else
                %get the only spike
                [valm indm] = max(v(inds(1:end)));
                spk(1) = indm + (inds(1) - 1);
            end;
        end;
        n_spk2(j,k) = length(spk);
    end;
end;

n_mean2 = mean(n_spk2,2);
s_mean2 = std(n_spk2,0,2);

% line('Parent',h,'XData',currents,'YData',n_mean2,'Marker','o','MarkerFaceColor','r');
% for i = 1:n_currents
%     line('Parent',gca,'XData',[currents(i) currents(i)],...
%         'YData',[n_mean2(i)-s_mean2(i) n_mean2(i)+s_mean2(i)],'Color','r'); 
% end;

n_spk3 = zeros(n_currents,n_trials);
syn_par(1) = 2.5*syn_par_def(1);
syn_par(3) = 2.5*syn_par_def(3);
syn_par(2) = 5*syn_par_def(2);
syn_par(4) = 5*syn_par_def(4);

for j = 1:n_currents
    %info_str = sprintf('processing current %i ...',j);
    %disp(info_str);
    isis = [];
    
    for k = 1:n_trials
        v = destexhy3(dt,10,1010,1020,currents(j),syn_par);

        %find spike times
        spk = [];
        inds = find(v > 0);
        if ( ~isempty(inds) )
            %at least one spike
            gaps = diff(inds);
            gi = find(gaps>1);
            if ( ~isempty(gi) )
                %at least two spikes
                [valm indm] = max(v(inds(1:gi(1))));
                spk(1) = indm + (inds(1)-1);

                for i =2:length(gi)
                    %if there are 3 or more spikes
                    [valm indm] = max(v(inds(gi(i-1)+1:gi(i))));
                    spk(i) = indm + (inds(gi(i-1)+1)-1);
                end;
        
                %last spike
                [valm indm] = max(v(inds(gi(end)+1:end)));
                spk(end+1) = indm + (inds(gi(end)+1)-1);
            else
                %get the only spike
                [valm indm] = max(v(inds(1:end)));
                spk(1) = indm + (inds(1) - 1);
            end;
        end;
        n_spk3(j,k) = length(spk);
    end;
end;

n_mean3 = mean(n_spk3,2);
s_mean3 = std(n_spk3,0,2);

% line('Parent',h,'XData',currents,'YData',n_mean3,'Marker','o','MarkerFaceColor','b');
% for i = 1:n_currents
%     line('Parent',gca,'XData',[currents(i) currents(i)],...
%         'YData',[n_mean3(i)-s_mean3(i) n_mean3(i)+s_mean3(i)],'Color','b'); 
% end;

