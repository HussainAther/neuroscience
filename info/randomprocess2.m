%{Three sample paths of a homogeneous Poisson process. 
Below each path the jump times are indicated by spikes. 
%}

m_t = 175e-3; %20 ms
t_max = 1; %time we sample process
n_rnd = 3*round(t_max/m_t); %3 times as many samples as needed to be on the safe side
v = exprnd(m_t,1,n_rnd);
v1 = cumsum(v);
v2 = v1(find(v1<t_max));
v3 = [0 v2];
n_v3 = length(v3);
y3 = [0:1:n_v3];

h_f1 = figure; 
h_a1 = subplot(6,1,1);
h_a2 = subplot(6,1,2);

for i = 1:n_v3-1
    line('Parent',h_a1,'XData',[v3(i) v3(i+1)],'YData',[y3(i) y3(i)]);
end;
line('Parent',h_a1,'XData',[v3(n_v3) 1],'YData',[y3(n_v3) y3(n_v3)]);

for i = 2:n_v3
    line('Parent',h_a2,'XData',[v3(i) v3(i)],'YData',[0 0.5]);
end;

set(h_a1,'XLim',[0 1]);
set(h_a2,'XLim',[0 1], 'YLim',[0 1]);

%%
v = exprnd(m_t,1,n_rnd);
v1 = cumsum(v);
v2 = v1(find(v1<t_max));
v3 = [0 v2];
n_v3 = length(v3);
y3 = [0:1:n_v3];

h_a3 = subplot(6,1,3);
h_a4 = subplot(6,1,4);

for i = 1:n_v3-1
    line('Parent',h_a3,'XData',[v3(i) v3(i+1)],'YData',[y3(i) y3(i)]);
end;
line('Parent',h_a3,'XData',[v3(n_v3) 1],'YData',[y3(n_v3) y3(n_v3)]);

for i = 2:n_v3
    line('Parent',h_a4,'XData',[v3(i) v3(i)],'YData',[0 0.5]);
end;

set(h_a3,'XLim',[0 1]);
set(h_a4,'XLim',[0 1],'YLim',[0 1]);

% set(handles.axes2,'XLim',[0 1]);
% set(handles.axes8,'XLim',[0 1]);
% set(handles.axes8,'YLim',[0 1]);

%%

v = exprnd(m_t,1,n_rnd);
v1 = cumsum(v);
v2 = v1(find(v1<t_max));
v3 = [0 v2];
n_v3 = length(v3);
y3 = [0:1:n_v3];

h_a5 = subplot(6,1,5);
h_a6 = subplot(6,1,6);

for i = 1:n_v3-1
    line('Parent',h_a5,'XData',[v3(i) v3(i+1)],'YData',[y3(i) y3(i)]);
end;
line('Parent',h_a5,'XData',[v3(n_v3) 1],'YData',[y3(n_v3) y3(n_v3)]);

for i = 2:n_v3
    line('Parent',h_a6,'XData',[v3(i) v3(i)],'YData',[0 0.5]);
end;

set(h_a5,'XLim',[0 1]);
set(h_a6,'XLim',[0 1],'YLim',[0 1]);
axes(h_a6);
xlabel('time (s)');

%%

h_f2 = figure; 
h_a7 = subplot(3,1,1);
h_a8 = subplot(3,1,2);
h_a9 = subplot(3,1,3);

%generate sample paths of inhomogeneous Poisson process with oscillating
% rate
dt = 0.1e-3; %time step in s

%one second
t_v = 0:dt:1-dt;
f_s = 4; %four cycles per s
y_v = 10*(sin(2*pi*f_s*t_v-pi/2) + 1);
line('Parent',h_a7,'XData',t_v,'YData',y_v);

n_ind = length(y_v);
y_vdt = y_v*dt; %speed up integration step below

%generate 10 ms bins (100*dt bins times 100 gives 1s)
m_vect = zeros(1,100);

%number of spike trains simulated
n_sims = 1000;
for sims = 1:n_sims
    %simulate an inhomogeneous Poisson spike train
    spk_ind = [];
    ind = 1;
    thres = exprnd(1);
    y_int = 0;
    while ( ind < n_ind )
        while ( (y_int < thres) && (ind < n_ind) )
            y_int = y_int + y_vdt(ind);
            ind = ind + 1;
        end;

        if ( y_int >= thres )
            spk_ind(end+1) = ind-1;
            %rounds to the nearest 10 ms bin
            mvind = floor((ind-1)/100)+1;
            m_vect(mvind) = m_vect(mvind)+1;
            thres = exprnd(1);
            y_int = 0;
        end;
    end;

    if (sims <=10)
        for i = 1:length(spk_ind)
            t_spk = t_v(spk_ind(i));
            line('Parent',h_a8,'XData',[t_spk t_spk],'YData',[0.1 0.9]+(sims-1));
        end;
    end;
end;

%average spike number in each bin
m_vect = m_vect/n_sims;
%convert to firing rate, 1/10 ms = 100 spk/s
m_vect = m_vect*100;

%generate a time vector centered at the middle of each
%bin, 5, 15, 25, ...
t_m = (1:100)*10e-3;
t_m = t_m - 5e-3;
line('Parent',h_a9,'XData',t_m,'YData',m_vect);
axes(h_a9);
xlabel('time (s)');
ylabel('instantaneous firing rate (spk/s)');
set(h_a9,'YLim',[0 25]);

%print(handles.figure1,'-depsc','rand_fig3.eps'); 
