%{Simulation of an integrate-and-fire neuron with random, gamma distributed 
threshold of order 10 stimulated with a sinusoidal current.
Plot hte prediction instantaneous rate of a corresponding unmodulated process.
}%
function rand_gamma

n_g = 10;

[h_fig, f_factor] = single_trial(n_g,1);

%time vector and increment
dt = 0.1e-3;
tv_s = [0:dt:1-dt];

%first compute the time transformed variable
f_s = 4; %four cycles per s
fr_max = 20; %in spk/s
y_v = (fr_max/2)*(sin(2*pi*f_s*tv_s-pi/2) + 1);
%plot(tv_s,y_v);

z_v = (fr_max/2)*(tv_s - (1/(2*pi*f_s))*cos(2*pi*f_s*tv_s-pi/2));
%figure; 
%plot(tv_s,z_v,'r');


%rate to compensate for order giving an overall rate of 1
rho = 10; %in spk/s
%mean rate in spk/sec
chi = rho/n_g;

%compute the chi_c function over the transformed time interval
tv_s2 = [0:dt:max(z_v)];

k_v = [1:1:9];
n_k = length(k_v);
zk = exp(2*pi*i*k_v/n_g);
t_mat = repmat(tv_s2,n_k,1);
t_mat2 = zeros(size(t_mat));
for j = 1:n_k
    t_mat2(j,:) = rho*(zk(j)-1)*t_mat(j,:);
end;
rhoc = chi*(1+ zk*exp(t_mat2));
rhoc_r = real(rhoc);

%figure; plot(tv_s2,rhoc_r);

%map true time to transformed time chi_c
alpha_v = interp1(tv_s2,rhoc_r,z_v);
full_rate = y_v.*alpha_v;

%figure; plot(tv_s,full_rate);

figure(h_fig);
subplot(3,1,1);
hold on;
plot(tv_s,full_rate,'r');

%uncomment to print figure
%print(h_fig,'-depsc2','randgamma_p1.eps')


function [h_fig, f_factor] = single_trial(n_g,plot_fig)
%generate sample paths of inhomogeneous gamma processes with oscillating
%rate
dt = 0.1e-3; %time step in s

%one second
t_v = 0:dt:1-dt;
f_s = 4; %four cycles per s
fr_max = 20; %in spk/s
y_v = (fr_max/2)*(sin(2*pi*f_s*t_v-pi/2) + 1);
if ( plot_fig == 1 )
    h_fig = figure; h_axes5 = subplot(3,1,3);  
    line('Parent',h_axes5,'XData',t_v,'YData',y_v);
    h_axes6 = subplot(3,1,2);
    h_axes4 = subplot(3,1,1);
else
    h_fig = [];
end;

n_ind = length(y_v);
y_vdt = y_v*dt; %speed up integration step below

%generate 10 ms bins (100*dt bins times 100 gives 1s)
m_vect = zeros(1,100);

%number of spike trains simulated
n_sims = 1000;

%order of the gamma distribution
%n_g = 10;

spk_count = zeros(1,n_sims);
for sims = 1:n_sims
    %simulate an inhomogeneous Poisson spike train
    if ( (sims <= 10) & (plot_fig == 1) )
        %save spike times for plot below
        spk_ind = [];
    end
    ind = 1;
    thres = gamrnd(n_g,1/n_g); %mean of 1
    y_int = 0;
    while ( ind < n_ind )
        while ( (y_int < thres) && (ind < n_ind) )
            y_int = y_int + y_vdt(ind);
            ind = ind + 1;
        end;

        if ( y_int >= thres )
            %time to spike...
            
            if ( (sims <= 10) & (plot_fig == 1) )
                %save spike time for plotting
                spk_ind(end+1) = ind-1;
            end;
            
            %rounds to the nearest 10 ms bin to compute mean 
            %instantaneous firing rate
            mvind = floor((ind-1)/100)+1;
            m_vect(mvind) = m_vect(mvind)+1;
            
            %add to current trial spike count
            spk_count(1,sims) = spk_count(1,sims)+1;
            
            %get a new random threshold
            thres = gamrnd(n_g,1/n_g);
            
            %reset voltage
            y_int = 0;
        end;
    end;

    %plot the first 10 spike trains
    if ( (sims <=10) & (plot_fig == 1) )
        for i = 1:length(spk_ind)
            t_spk = t_v(spk_ind(i));
            line('Parent',h_axes6,'XData',[t_spk t_spk],'YData',[0.1 0.9]+(sims-1));
        end;
    end;
end;

%average spike number in each bin
m_vect = m_vect/n_sims;
%convert to firing rate, 1/10 ms = 100 spk/s
m_vect = m_vect*100;

%mean and variance of the spike count
m_spk_count = mean(spk_count);
v_spk_count = var(spk_count);
f_factor = v_spk_count/m_spk_count;

if ( plot_fig == 1 )
    %generate a time vector centered at the middle of each
    %bin, 5, 15, 25, ...
    t_m = (1:100)*10e-3;
    t_m = t_m - 5e-3;
    line('Parent',h_axes4,'XData',t_m,'YData',m_vect);
    axes(h_axes5);
    xlabel('time (s)');
    ylabel('instantaneous firing rate (spk/s)');
    set(h_axes5,'YLim',[0 fr_max+5],'TickDir','out');
    set(h_axes4,'YLim',[0 fr_max+5],'TickDir','out');
    set(h_axes6,'TickDir','out');
end;

