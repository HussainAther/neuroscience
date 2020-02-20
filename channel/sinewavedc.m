function  sinewavedc
%Computes the correlation model output for a translating sinewave stimulus with addition of a constant luminance component. 

%set the value of the phase of the cosine wave
phi = 0;

%stimulus contrast
alpha = 0.8;

%spatial wavelength
omega = 0.05; %c/deg

t_start = 0;
t_end = 1; %s
dt = 0.05e-3; %s
t_vect = t_start:dt:t_end;

%grating translation speed
v = 640; %deg/s

%position of the first sampling station
x0 = 0;

lum_vect0 = alpha*cos(2*pi*(omega*x0-omega*v*t_vect+phi));

%distance to the second sampling station
dx = 2; %deg
x1 = x0+dx;
lum_vect1 = alpha*cos(2*pi*(omega*x1-omega*v*t_vect+phi));

%compute responses of the 1st order low-pass filters
tau_s = 50e-3; %50 ms
tau_f = 1e-3; %1 ms
lv0_s_num = lp_filt_fn(lum_vect0,dt,tau_s);
lv0_f_num = lp_filt_fn(lum_vect0,dt,tau_f);
lv1_s_num = lp_filt_fn(lum_vect1,dt,tau_s);
lv1_f_num = lp_filt_fn(lum_vect1,dt,tau_f);

%baseline luminance/contrast value
i_l = 0.1;

%half correlator responses
r_hc1 = (i_l + lv0_s_num).*(i_l + lv1_f_num);
r_hc2 = (i_l + lv1_s_num).*(i_l + lv0_f_num);

%full opponent response
r_hr = r_hc1 - r_hc2;

%compare with the standard correlation model response
r_hc1_std = lv0_s_num.*lv1_f_num;
r_hc2_std = lv1_s_num.*lv0_f_num;

%full opponent response
r_hr_std = r_hc1_std - r_hc2_std;


%steady-state initial time
t_ss = 0.5; 

%compute steady-state oscillation amplitude
sso_amp = max(r_hr(t_vect>=t_ss))-min(r_hr(t_vect>=t_ss));

%compute steady-state signal amplitude
sss_amp = mean(r_hr(t_vect>=t_ss));

%compute signal to noise ratio
snr_hr = sss_amp/sso_amp;

%setup an array of correlation detectors to check that the steady state
%response is again constant
r_hr_mean = zeros(size(t_vect));
r_hr_std_mean = zeros(size(t_vect));

n_corr = 72;
x_vect = (0:n_corr)*dx;
for i = 1:n_corr
    %contrast input for the two detectors
    lum_vect0 = alpha*cos(2*pi*(omega*x_vect(i)-omega*v*t_vect+phi));
    lum_vect1 = alpha*cos(2*pi*(omega*x_vect(i+1)-omega*v*t_vect+phi));
    
    %low pass filtering
    lv0_s_num = lp_filt_fn(lum_vect0,dt,tau_s);
    lv0_f_num = lp_filt_fn(lum_vect0,dt,tau_f);
    lv1_s_num = lp_filt_fn(lum_vect1,dt,tau_s);
    lv1_f_num = lp_filt_fn(lum_vect1,dt,tau_f);
    
    %add DC component and compute half-correlator responses
    r_hc1 = (i_l + lv0_s_num).*(i_l + lv1_f_num);
    r_hc2 = (i_l + lv1_s_num).*(i_l + lv0_f_num);
    
    %full opponent response
    r_hr = r_hc1 - r_hc2;
    
    %add to the mean vector
    r_hr_mean = r_hr_mean + r_hr;

    %Repeat for the standard correlation model
    r_hc1_std = lv0_s_num.*lv1_f_num;
    r_hc2_std = lv1_s_num.*lv0_f_num;
    
    %full opponent response
    r_hr_std = r_hc1_std - r_hc2_std;

    %add to the mean vector
    r_hr_std_mean = r_hr_std_mean + r_hr_std;

end

%compute average
r_hr_mean = r_hr_mean/n_corr;
r_hr_std_mean = r_hr_std_mean/n_corr;

sso_amp_mean = max(r_hr_mean(t_vect>=t_ss))-min(r_hr_mean(t_vect>=t_ss));

%compute the decrease in oscillation amplitude caused by averaging relative
%to baseline
disp(['baseline oscillation amplitude: ' num2str(sso_amp) ', after spatial averaging: ' num2str(sso_amp_mean) ', ratio: ' num2str(sso_amp_mean/sso_amp)]);

%compute a signal to noise ratio, defined as the mean response of the
%detector divided by the oscillation amplitude
sss_amp_mean = mean(r_hr_mean(t_vect>=t_ss));
snr_hr_mean = sss_amp_mean/sso_amp_mean;

disp(['single correlator signal to noise ratio: ' num2str(snr_hr) 'snr after spatial averaging: ' num2str(snr_hr_mean)]);

h_f = figure; 
ha1 = subplot(2,1,1);
hp1 = plot(t_vect,r_hr,'r');
hold on;
hp2 = plot(t_vect, r_hr_std,'k');
legend([hp1 hp2],{'DC correlator response', 'standard model'});
set(ha1,'TickDir','out');
set(ha1,'XLim',[0 1],'YLim',[-0.25 0.45]);

ha2 = subplot(2,1,2);
hp3 = plot(t_vect,r_hr_mean,'k');
legend([hp3],{'DC spatially averaged response'});
%hold on;
%hp4 = plot(t_vect,r_hr_std_mean);
%legend([hp3 hp4],{'DC spatially averaged response', 'standard model spatially averaged'});
set(ha2,'TickDir','out');
xlabel('time (s)');
ylabel('model output');
set(ha2,'XLim',[0 1],'YLim',[-0.25 0.45]);

%print(h_f,'figures/dc_component.eps','-depsc');

end

function lp_vect = lp_filt_fn(lum_vect,dt, tau_lp)
% Numerical implementation of first order low-pass filter for correlator models
%
% Usage:   vect = lp1filt_fn(lum_vec,dt,tau)
%
% where lum_vect is to be filtered,
% tau_lp is the time constant in s,
% dt is time step in s. 

if ( nargin ~= 3 )
    disp('Unexpected number of inputs. Aborting...');
    return;
end

% Model parameters

if ( isempty(tau_lp) )
    %time constants of 1st order low pass (ms)
    tau_lp = 20e-3; %s
end

if ( isempty(dt) )
    dt = 0.05e-3; %s
end

%setup constants for backward euler
lp_const1 =  dt/(tau_lp + dt);
lp_const2 = 1 - lp_const1;

lp_vect = zeros(size(lum_vect));

%initial condition 
lp_vect(1) = lum_vect(1);

%marching algorithm
for i = 2:length(lum_vect)
    lp_vect(i) = lp_const1*lum_vect(i)+lp_const2*lp_vect(i-1);
end
end
