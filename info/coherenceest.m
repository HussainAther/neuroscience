% Compute analytically the coherence between an OU (Ornstein–Uhlenbeck process) process, X and Z=X+Y, where Y is an independent white 
% Gaussian noise, up to a cut-off frequency fN. The OU process is assumed to have parameters σOU and τ, 
% while the white noise is assumed to have variance σw2.
% Ornstein–Uhlenbeck process (Ornstein Uhlenbeck)
%relaxation time constant and steady-state SD
tau_ou = 2.7; %ms
sig_ou = 0.55; %mu S

%time step and time vector
dt = 0.5; %ms
fs = 1/(dt*1e-3); %sampling freq
fn = fs/2; %nyquist freq

%this gives 16 segments to average over for the 
%power spectrum
n_seg = 32;
seg_len = 2048;
n_v = 0:(n_seg*seg_len-1);
t_v = n_v*dt;
n_t_v = length(t_v);

%numerical simulation factors
nsf_1 = exp(-dt/tau_ou);
nsf_2 = sqrt(1-exp(-2*dt/tau_ou));

%initial value of the random conductance variable
x0 = 0;

%random white noise increment
w_vou = normrnd(0,1,1,n_t_v);
u_vou = sig_ou*nsf_2*w_vou;

%save space for random vector
x_vou = zeros(1,n_t_v);

%initial condition
x_vou(1) = x0;

for i =2:n_t_v
    %stochastic update
    x_vou(i) = x_vou(i-1)*nsf_1 + u_vou(i);
end;

%white noise
sigma_w = 1;
y_vw = randn(size(x_vou));

z = x_vou + y_vw;

h_f1 = figure; 
h_a1 = subplot(2,1,1);
%show data
line('Parent',h_a1,'XData',t_v(1:250),'YData',z(1:250));
line('Parent',h_a1,'XData',t_v(1:250),'YData',x_vou(1:250),'Color','r');
set(h_a1,'XLim',[0 125],'YLim',[-3 3]);
xlabel(h_a1,'time (ms)');

%% 

%debug
%2048 = window length (hamming)
%1024 = overlap (50 percent)
%2048 = number of samples for fft
%2000 = sampling frequency (0.5 ms = 2kHz)
%[px_w,f_w] = pwelch(x_vou,seg_len,seg_len/2,seg_len,fs,'twosided');
%line('Parent',handles.axes4,'XData',f_w(1:seg_len/2 + 1),'YData',px_w(1:seg_len/2+1),'Color','r');

%[py_w,f_w] = pwelch(y_vw,seg_len,seg_len/2,seg_len,fs,'twosided');
%line('Parent',handles.axes4,'XData',f_w(1:seg_len/2 + 1),'YData',py_w(1:seg_len/2+1));

%[pz_w,f_w] = pwelch(z,seg_len,seg_len/2,seg_len,fs,'twosided');
%line('Parent',handles.axes4,'XData',f_w(1:seg_len/2 + 1),'YData',pz_w(1:seg_len/2+1),'Color','k');

%theoretical power spectrum OU process
f = 0:0.1:1000;
pfou = 2*sig_ou^2*tau_ou*1e-3./(1 + (2*pi*f*tau_ou*1e-3).^2);
%line('Parent',handles.axes4,'XData',f,'YData',pfou,'Color','g');

%theoretical power spectrum white noise
pfw = sigma_w^2/(2*fn)*ones(size(f));

%line('Parent',handles.axes4,'XData',f,'YData',pfw,'Color','b');

%theoretical power spectrum sum of OU and white noise
pfz = pfw + pfou;
%line('Parent',handles.axes4,'XData',f,'YData',pfz,'Color','r');

%top_ylim = 10e-8;
%set(handles.axes4,'XLim',[0 200]); %,'YLim',[0 top_ylim]);

nfft = 2048;
noverlap = 1024;
window = 2048;
[Cxy,F] = mscohere(z,x_vou,window,noverlap,nfft,fs);
h_a2 = subplot(2,1,2);
line('Parent',h_a2,'XData',F,'YData',Cxy);
xlabel(h_a2,'frequency (Hz)');
ylabel(h_a2,'coherence');

%theoretical mean square coherence
cxy_t = pfou./pfz;
line('Parent',h_a2,'XData',f,'YData',cxy_t,'Color','r');

%%
%
% modified from whiten1
%

% Generate an oversampled white noise segment in the frequency domain
% Passes the white noise through a static non-linearity. Computes 
% autocorrelations of the input and outpout, cross correlation and 
% corresponding quantities in the Fourier domain. Analytical results
% are plotted as well.
   
%
% Generate the white noise
%

stddev = 1;
fc = 100; %cut-off frequency Hz
deltat = 1e-3; %sampling step, converts from msec in sec
fs = 1/deltat; %sampling frequency
seg_len = 1024;
n = 32*seg_len; %number of samples
sigma = sqrt((stddev^2*n)/(4*fc*deltat));
 
%typically, a random seed controls the state of the random number 
%generator. To make the following sequence reproducible it should
%be initialized here 

whitef = zeros(1,n);
   
%zero frequency component
ai = randn(1) * sigma;
whitef(1) = complex(ai,0);

%all component except for nyquist
for i=2:n/2
    if ( i/(n*deltat) < fc )
        ai = randn(1) * sigma;
        bi = randn(1) * sigma;
        whitef(i) = complex(ai,bi);
        whitef(n-i+2) = complex(ai,-bi);
    end;
end;

ai = randn(1) * sigma;
whitef(n/2+1) = complex(ai,0);

whitet = fft(whitef);
whitet = whitet/n;
whitet = real(whitet);

%
%autocorrelation function of the input
%

%analytical formula
n_tp = 128;
t_v = (-n_tp:n_tp)* deltat;
ind_z = find(t_v == 0);
auto_th = sin(2*pi*fc*t_v)./(2*pi*fc*t_v);
auto_th(ind_z) = 1;

%
%power spectrum of input
%

%computed from the FT of the autocorrelation
df = 1/(2*n_tp*deltat);
f_auto = (0:2*n_tp-1)*df;
auto_sf = circshift(auto_th(2:end),[0 -(n_tp-1)]);
ps_est1 = real(fft(auto_sf))*deltat;

%
% static non-linearity
%

%gaussian signal
x = -5:0.01:5;
y = normpdf(x,0,1);

l2 = 1/10;
nl2 = erf(x/(sqrt(2)*l2));
nw2 = erf(whitet/(sqrt(2)*l2));

[n,xout] = hist(whitet,100);
dxout = xout(2) - xout(1);
n = n/(sum(n*dxout));

h_f2 = figure;
h_a3 = subplot(1,3,1);
bar(h_a3,xout,n);
line('Parent',h_a3,'XData',x,'YData',y,'Color','r');
set(h_a3,'TickDir','out');

h_a4 = subplot(1,3,2);
line('Parent',h_a4,'XData',x,'YData',nl2);
set(h_a4,'XLim',[-1 1]);

[n,xout] = hist(nw2,100);
dxout = xout(2) - xout(1);
n = n/(sum(n*dxout));

h_a5 = subplot(1,3,3);
bar(h_a5,xout,n);
set(h_a5,'XLim',[-1.1 1.1],'YLim',[0 25],'TickDir','out');

%%
n_pts = 250;
t_v_ex = (1:n_pts)*dt;

h_f3 = figure; 
h_a6 = subplot(2,1,1);
line('Parent',h_a6,'XData',t_v_ex,'YData',whitet(1:n_pts),'Color','r');
line('Parent',h_a6,'XData',t_v_ex,'YData',nw2(1:n_pts),'Color','k');
set(h_a6,'XLim',[0 125],'YLim',[-3 3]);
xlabel(h_a6,'time (ms)');

%%
%
%autocorrelation of output
%
% %analytical formula
auto_th2 = (2/pi)*asin(auto_th/(1 + l2^2));

%
%power spectrum of output
%

%analytical forumula from autocorrelation function
auto_sf2 = circshift(auto_th2(2:end),[0 -(n_tp-1)]);
ps_est2 = real(fft(auto_sf2))*deltat;

%
%cross correlation
%
cross_th2 = (sqrt(2)/sqrt(pi))*(1/sqrt(1 + l2^2))*auto_th;

%
%cross spectrum
%

%value from cross correlation
cross_sf = circshift(cross_th2(2:end),[0 -(n_tp-1)]);
cs_est = (abs(fft(cross_sf))*deltat).^2;

%
%mean squared coherence
%

%seg_len = window 
%seg_len/2 = noverlap
%seg_len = nfft
[Cxy,F] = mscohere(nw2,whitet,seg_len,seg_len/2,seg_len,1000);
 
%theoretical value
coh_th = cs_est./(ps_est1.*ps_est2);
% 
%plot ms coherence
h_a7 = subplot(2,1,2);
line('Parent',h_a7,'XData',F,'YData',Cxy);
line('Parent',h_a7,'XData',f_auto(1:n_tp+1),'YData',coh_th(1:n_tp+1),'Color','r');
set(h_a7,'YLim',[0 1]);
xlabel(h_a7,'frequency (Hz)');

%print(handles.figure1,'-depsc2','coherence_est.eps');

