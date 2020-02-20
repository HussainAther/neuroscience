%{
The autocovariance of Y is C_YY = (2alpha^2/pi) * arcsin(C_XX/sigma_X^2+l^2)
Use this result and Exercise 7 to derive numerically the coherence between X and Y when X is white noise with a cut-off frequency of 100 Hz. Show that this allows you to reproduce the red curve in Figure 18.1D. 
%}

% Coherence estimation
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
n_t_v = length(t_v)

%numerical simulation factors
nsf_1 = exp(-dt/tau_ou);
nsf_2 = sqrt(1-exp(-2*dt/tau_ou));

%initial value of the random conductance variable
x0 = 0;;

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
