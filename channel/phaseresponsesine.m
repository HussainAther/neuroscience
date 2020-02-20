%%
%function phaseresponsesine
%Plots the steady state response to sinewave stimuli for the correlator output as a function of pattern speed.

alpha = 1;
tau_s = 50e-3; %s
tau_f = 1e-3; %s

%spatial frequency
omega_x1 = 1/8; %c/deg
omega_x2 = 7/8;

%compute the stimuli
dx = 0.05; %spatial resolution
x_vect = 0:dx:10; %deg
c_vect1 = cos(2*pi*omega_x1*x_vect);
c_vect2 = cos(2*pi*omega_x2*x_vect);

%phase shift of 2nd stimulus 
phase = pi;
c_vect2b = cos(2*pi*omega_x2*x_vect + phase);

%mixing both stimuli
alpha3 = 0.5;
beta3 = 0.5;
c_vect3 = alpha3*c_vect1 + beta3*c_vect2;
c_vect3b = alpha3*c_vect1 + beta3*c_vect2b;

alpha4 = 0.8;
beta4 = 0.2;
c_vect4 = alpha4*c_vect1 + beta4*c_vect2;
c_vect4b = alpha4*c_vect1 + beta4*c_vect2b;

%plot the two stimuli
hf2 = figure;
ha2 = axes;
hp4 = plot(x_vect,c_vect3,'k');
hold on;
hp5 = plot(x_vect,c_vect3b,'r');
legend([hp4 hp5],{'equally weighted sum', 'phase shifted sum'});
xlabel('spatial position (deg)');
ylabel('contrast');
set(ha2,'TickDir','out');

%print(hf2,'figures/phase_resp1.eps','-depsc');

hf3 = figure; 
ha3 = axes;
hp6 = plot(x_vect,c_vect4,'k');
hold on;
hp7 = plot(x_vect,c_vect4b,'r');
legend([hp6 hp7],{'unequally weighted sum', 'phase shifted sum'});
xlabel('spatial position (deg)');
ylabel('contrast');
set(ha3,'TickDir','out');

%print(hf3,'figures/phase_resp2.eps','-depsc');

%%
%spacing between two sampling stations
dx_hr = 2; %in deg

%angular speed
v_vect = [0.1:0.1:1000]; %deg/s

%steady-state response for one sinewave
r_ss_th1 = alpha^2*2*pi*omega_x1*(tau_s-tau_f)*sin(2*pi*omega_x1*dx_hr)*v_vect./...
    ((1+(2*pi*omega_x1*tau_s*v_vect).^2).*(1+(2*pi*omega_x1*tau_f*v_vect).^2));

%steady-state response for one sinewave
r_ss_th2 = alpha^2*2*pi*omega_x2*(tau_s-tau_f)*sin(2*pi*omega_x2*dx_hr)*v_vect./...
    ((1+(2*pi*omega_x2*tau_s*v_vect).^2).*(1+(2*pi*omega_x2*tau_f*v_vect).^2));

r_ss_th3 = alpha3^2*r_ss_th1 + beta3^2*r_ss_th2;
r_ss_th4 = alpha4^2*r_ss_th1 + beta4^2*r_ss_th2;

hf5 = figure;
ha5 = axes;
semilogx(v_vect,r_ss_th3,'k');
hold on;
semilogx(v_vect,r_ss_th4,'r');
xlabel('angular speed (deg/s)');
ylabel('model response');
title('steady state response of correlator to sums of translating sinewaves');
set(ha5,'TickDir','out');

%% first cosine wave parameters

%compute the numerical response to the stimuli
%v_num = 20; %deg/s, test case
v_num = [0.5 1 2 5 10 20 50 100 200 500];

%number of temporal samples
n_t = 1000;

%temporal frequency
omega_t1 = v_num*omega_x1; % c/deg

%temporal wavelength
lambda_t1 = 1./omega_t1; %s

%maximal time of simulation; 10 times wavelength
t_max1 = 10*lambda_t1;

%100 points per wavelength
dt1 = t_max1/n_t;

%location of the sampling stations
x_1 = 0; 
x_2 = dx_hr; 

debug = 0; 

%% second cosinewave parameters

%temporal frequency
omega_t2 = v_num*omega_x2; % c/deg

%temporal wavelength
lambda_t2 = 1./omega_t2; %s

%maximal time of simulation; 10 times wavelength
t_max2 = 10*lambda_t2;

%100 points per wavelength
dt2 = t_max2/n_t;

%% third cosine wave

%frequency of the second wave is 7 times higher so temporal wavelength is 7 times shorter for
%the second wave. Hence to use the same length as in the first wave, we
%need 7 times more sampling points. 
n_t2 = 7*n_t;

%compute steady-state value; last seven times two hundred points
ind_frame2 = ((n_t2-(7*200)+1):n_t2);

for i = 1:length(v_num)
    
    %time vector
    t_vect1 = (0:(n_t2-1))*dt2(i);
    
    %compute low pass filtered version of stimulus
    ct_vect1a = cos(2*pi*(omega_x1*x_1 - omega_t1(i)*t_vect1));
    ct_vect2a = cos(2*pi*(omega_x1*x_2 - omega_t1(i)*t_vect1));

    %repeat for second stimulus
    ct_vect1b = cos(2*pi*(omega_x2*x_1 - omega_t2(i)*t_vect1));
    ct_vect2b = cos(2*pi*(omega_x2*x_2 - omega_t2(i)*t_vect1));

    ct_vect1 = alpha3*ct_vect1a + beta3*ct_vect1b;
    ct_vect2 = alpha3*ct_vect2a + beta3*ct_vect2b;
    
    %repeat for phase shifted second stimulus
    ct_vect1bp = cos(2*pi*(omega_x2*x_1 - omega_t2(i)*t_vect1) + phase);
    ct_vect2bp = cos(2*pi*(omega_x2*x_2 - omega_t2(i)*t_vect1) + phase);

    ct_vect1p = alpha3*ct_vect1a + beta3*ct_vect1bp;
    ct_vect2p = alpha3*ct_vect2a + beta3*ct_vect2bp;
    
    if ( debug == 1 )
        hf6 = figure;
        plot(t_vect1,ct_vect1);
        hold on;
        plot(t_vect1,ct_vect2);
    end
    
    %compute the lowpass filtered version of the stimuli
    ls_ct_vect1 = lp1filt_fn(ct_vect1,dt2(i),tau_s);
    lf_ct_vect1 = lp1filt_fn(ct_vect1,dt2(i),tau_f);
    
    ls_ct_vect2 = lp1filt_fn(ct_vect2,dt2(i),tau_s);
    lf_ct_vect2 = lp1filt_fn(ct_vect2,dt2(i),tau_f);
    
    hr_out3 = ls_ct_vect1.*lf_ct_vect2 - ls_ct_vect2.*lf_ct_vect1;
    
    if ( debug == 1 )
        hf7 = figure;
        plot(t_vect1,hr_out3);
    end
    
    hr_ss3(i) = mean(hr_out3(ind_frame2));

    %repeat for phase shifted stimulus
    ls_ct_vect1p = lp1filt_fn(ct_vect1p,dt2(i),tau_s);
    lf_ct_vect1p = lp1filt_fn(ct_vect1p,dt2(i),tau_f);
    
    ls_ct_vect2p = lp1filt_fn(ct_vect2p,dt2(i),tau_s);
    lf_ct_vect2p = lp1filt_fn(ct_vect2p,dt2(i),tau_f);
    
    hr_out3p = ls_ct_vect1p.*lf_ct_vect2p - ls_ct_vect2p.*lf_ct_vect1p;
    hr_ss3p(i) = mean(hr_out3p(ind_frame2));

end

figure(hf5);
plot(v_num,hr_ss3,'xk');
plot(v_num,hr_ss3p,'ok');

%% fourth cosine wave

for i = 1:length(v_num)
    
    %time vector
    t_vect1 = (0:(n_t2-1))*dt2(i);
    
    %compute low pass filtered version of stimulus
    ct_vect1a = cos(2*pi*(omega_x1*x_1 - omega_t1(i)*t_vect1));
    ct_vect2a = cos(2*pi*(omega_x1*x_2 - omega_t1(i)*t_vect1));

    %repeat for second stimulus
    ct_vect1b = cos(2*pi*(omega_x2*x_1 - omega_t2(i)*t_vect1));
    ct_vect2b = cos(2*pi*(omega_x2*x_2 - omega_t2(i)*t_vect1));

    ct_vect1 = alpha4*ct_vect1a + beta4*ct_vect1b;
    ct_vect2 = alpha4*ct_vect2a + beta4*ct_vect2b;
    
    %repeat for phase shifted second stimulus
    ct_vect1bp = cos(2*pi*(omega_x2*x_1 - omega_t2(i)*t_vect1) + phase);
    ct_vect2bp = cos(2*pi*(omega_x2*x_2 - omega_t2(i)*t_vect1) + phase);

    ct_vect1p = alpha4*ct_vect1a + beta4*ct_vect1bp;
    ct_vect2p = alpha4*ct_vect2a + beta4*ct_vect2bp;
    
    if ( debug == 1 )
        hf6 = figure;
        plot(t_vect1,ct_vect1);
        hold on;
        plot(t_vect1,ct_vect2);
    end
    
    %compute the lowpass filtered version of the stimuli
    ls_ct_vect1 = lp1filt_fn(ct_vect1,dt2(i),tau_s);
    lf_ct_vect1 = lp1filt_fn(ct_vect1,dt2(i),tau_f);
    
    ls_ct_vect2 = lp1filt_fn(ct_vect2,dt2(i),tau_s);
    lf_ct_vect2 = lp1filt_fn(ct_vect2,dt2(i),tau_f);
    
    hr_out4 = ls_ct_vect1.*lf_ct_vect2 - ls_ct_vect2.*lf_ct_vect1;
    
    if ( debug == 1 )
        hf7 = figure;
        plot(t_vect1,hr_out4);
    end
    
    hr_ss4(i) = mean(hr_out4(ind_frame2));

    %repeat for phase shifted stimulus
    ls_ct_vect1p = lp1filt_fn(ct_vect1p,dt2(i),tau_s);
    lf_ct_vect1p = lp1filt_fn(ct_vect1p,dt2(i),tau_f);
    
    ls_ct_vect2p = lp1filt_fn(ct_vect2p,dt2(i),tau_s);
    lf_ct_vect2p = lp1filt_fn(ct_vect2p,dt2(i),tau_f);
    
    hr_out4p = ls_ct_vect1p.*lf_ct_vect2p - ls_ct_vect2p.*lf_ct_vect1p;
    hr_ss4p(i) = mean(hr_out4p(ind_frame2));

end

figure(hf5);
plot(v_num,hr_ss4,'xr');
plot(v_num,hr_ss4p,'or');
%print(hf5,'figures/phase_resp3.eps','-depsc');

%end
