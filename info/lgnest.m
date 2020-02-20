% LGN (lateral geniculate nucleus) neuron estimation

fc = 5.5; %Hz
fc1 = 2*pi*fc*1e-3; %converts to kHz and circular frequency

%sampling step in msec
dt = 1;
fs = 1/(dt*1e-3); %sampling frequency
fn = fs/2; %nyquist frequency
nlgn = 512;
t_lgn = (0:nlgn-1)*dt; %in msec
v_lgn1 = exp(-fc1*t_lgn);
v_lgn = t_lgn.*v_lgn1 - (fc1/2) * (t_lgn.^2) .* v_lgn1;

max_fr = 50;
v_lgn = max_fr*(v_lgn./max(v_lgn));

nfft = 1024;
nw = 32*nfft;
sigma_w = 0.25;
w_n = sigma_w*randn(1,nw);

%the matlab definition of the convolution requires the lgn 
%filter to be first and the stimulus to be second. To get
%no phase delay we need to take the full wave form and truncate
%it, as seen by convolving v_lgn with a [1 0 0 ... 0] sequence
f_v = conv(v_lgn,w_n);
f_v = f_v(1:length(w_n));

[gf_lgn,f_lgn] = tfestimate(w_n,f_v,nfft,nfft/2,nfft,fs,'twosided');

gt_lgn = ifft(gf_lgn);

h_f1 = figure; 
h_a1 = subplot(2,1,1);
line('Parent',h_a1,'XData',t_lgn,'YData',gt_lgn(1:nlgn));
line('Parent',h_a1,'XData',t_lgn,'YData',v_lgn,...
    'LineStyle','--','Color','r');
set(h_a1,'XLim',[0 t_lgn(256)],'YLim',[-20 60]);
%%

%static non-linearity
s2 = 10;
x = 0:1:200;
y = 100*x./(s2 + x);

f_vp = max(0,f_v);
f_vn = 100*f_vp./(s2 + f_vp);

h_f2 = figure; 
h_a2 = subplot(2,1,1);
h_a3 = subplot(2,1,2);
nplot = 1000;
line('Parent',h_a2,'XData',(1:nplot)*dt,'YData',w_n(1:nplot));
line('Parent',h_a3,'XData',(1:nplot)*dt,'YData',f_v(1:nplot));
line('Parent',h_a3,'XData',(1:nplot)*dt,'YData',f_vn(1:nplot),...
    'Color','r');
set(h_a3,'YLim',[-200 200]);

%%
[gf2_lgn,f_lgn] = tfestimate(w_n,f_vn,nfft,nfft/2,nfft,fs,'twosided');
 
gt2_lgn = ifft(gf2_lgn);

%proportionality factor
alpha = mean(f_v.*f_vn)/var(f_v);
gt2_lgn = gt2_lgn/alpha;

figure(h_f1);
h_a4 = subplot(2,1,2);
line('Parent',h_a4,'XData',t_lgn,'YData',gt2_lgn(1:nlgn));
line('Parent',h_a4,'XData',t_lgn,'YData',v_lgn,...
    'LineStyle','--','Color','r');
set(h_a4,'XLim',[0 t_lgn(256)],'YLim',[-20 60]);

%%
%recovery of the static non-linearity
f_v2 = conv(gt2_lgn,w_n);
f_v2 = f_v2(1:length(w_n));

my = mean(f_v2);
sy = std(f_v2);
yv = my-5*sy:sy/100:my+5*sy;
fy = normcdf(yv,my,sy);

[n, xout] = hist(f_vn,100);
nn = n/sum(n);
cpdfz = cumsum(nn);

snl = zeros(size(yv));

for i = 1:length(snl)
    if ( fy(i) < cpdfz(1) )
       snl(i) = xout(1);
    elseif ( fy(i) > cpdfz(end) )
        snl(i) = xout(end);
    else
        %inbetween value
        snl(i) = interp1(cpdfz,xout,fy(i));
    end;
end;
 
h_f3 = figure; 
h_a5 = axes;
line('Parent',h_a5,'XData',yv,'YData',snl);
line('Parent',h_a5,'XData',x,'YData',y, ...
    'Color','r','LineStyle','--');
set(h_a5,'XLim',[0 200]);

xlabel(h_a1,'time (ms)');
xlabel(h_a4,'time (ms)');
ylabel(h_a3,'firing rate (spk/s)');
xlabel(h_a2,'time (ms)');
ylabel(h_a1,'firing rate change (spk/s)');
ylabel(h_a4,'firing rate change (spk/s)');
xlabel(h_a3,'time (ms)');

%print(handles.figure1,'-depsc2','lgn_est3.eps');
