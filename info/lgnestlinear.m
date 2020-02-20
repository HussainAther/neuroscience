% Linear estimation of LGN (lateral geniculate nucleus) transfer function from firing rate
fc = 5.5; %Hz
fc1 = 2*pi*fc*1e-3; %converts to kHz and circular frequency

%sampling step in msec
dt = 10;
fs = 1/(dt*1e-3); %sampling frequency
fn = fs/2; %nyquist frequency
nlgn = 32;
t_lgn = (0:nlgn-1)*dt; %in msec
v_lgn1 = exp(-fc1*t_lgn);
v_lgn = t_lgn.*v_lgn1 - (fc1/2) * (t_lgn.^2) .* v_lgn1;

max_fr = 50;
v_lgn = max_fr*(v_lgn./max(v_lgn));


ni = 128;
nfft = 64;
nw = ni*nfft;
sigma_w = 0.25;
w_n = sigma_w*randn(1,nw);


%the matlab definition of the convolution requires the lgn 
%filter to be first and the stimulus to be second. To get
%no phase delay we need to take the full wave form and truncate
%it, as seen by convolving v_lgn with a [1 0 0 ... 0] sequence
f_v = conv(v_lgn,w_n);
f_v = f_v(1:length(w_n));


%[tlgn,flgn] = tfestimate(w_n,f_v,nfft,nfft/2,nfft,fs,'twosided');
%glgn = real(ifft(tlgn));

%generates an ni by nfft matrix, with each line being
%l_1:  e_1 ... e_nfft
%l_ni: e_((ni-1)*nfft + 1) ... e(ni*nfft)
wn2 = reshape(w_n',nfft,ni)';

%compute the covariance of the WN
cwn2 = cov(wn2);

%compute the covariance of the firing rate
fv2 = zeros(ni,nfft);
for i = 1:ni
    m_inter = conv(v_lgn,wn2(i,:));
    fv2(i,:) = m_inter(1:nfft);
end;

c_fv2_wn2 = fv2(1,:)'*wn2(1,:);
for i = 2:ni
    c_fv2_wn2 = c_fv2_wn2 + fv2(i,:)'*wn2(i,:);
end;

c_fv2_wn2 = c_fv2_wn2/(ni-1);

%compute the reconstruction matrix
h_mat = c_fv2_wn2*inv(cwn2);

h_f1 = figure; 
h_a1 = axes;
%compute the mean across diagonals and plot individual 
%diagonal elements
glgn2 = zeros(1,nlgn);
for i = 1:nlgn
    d_v = diag(h_mat,1-i);
    glgn2(i) = mean(d_v);
    
    line('Parent',h_a1,'XData',t_lgn(i)*ones(size(d_v)),'YData',d_v,...
        'LineStyle','none','Marker','.','MarkerSize',4);    
end;

line('Parent',h_a1,'XData',t_lgn,'YData',v_lgn);
line('Parent',h_a1,'XData',t_lgn,'YData',glgn2,'Color','r');
set(h_a1,'XLim',[0 320],'YLim',[-30 60]);
xlabel(h_a1,'time (ms)');
ylabel(h_a1,'firing rate change (spk/s)');

%%
%test to make sure we got everything right
pred = h_mat*wn2(1,:)';

h_f2 = figure; 
h_a2 = subplot(2,1,1);
h_a3 = subplot(2,1,2);
%plot the raw data and reconstruction
line('Parent',h_a2,'XData',(0:nfft-1)*dt,'YData',wn2(1,:));
line('Parent',h_a3,'XData',(0:nfft-1)*dt,'YData',fv2(1,:));
line('Parent',h_a3,'XData',(0:nfft-1)*dt,'YData',pred,'Color','r');
set(h_a2,'XLim',[0 640]);
set(h_a3,'XLim',[0 640],'YLim',[-60 60]);
xlabel(h_a3,'time (ms)');
ylabel(h_a3,'firing rate change (spk/s)');

%%
h_f3 = figure; 
h_a4 = axes;
colormap('gray');
imagesc(h_mat);
colorbar;
%surf(handles.axes1,h_mat');
%colormap(handles.axes1,'gray');
%axes(handles.axes1);
%view(90,90);
%colorbar;
set(h_a4,'XLim',[1 64],'YLim',[1 64],'TickDir','out');
xlabel(h_a4,'column index');
ylabel(h_a4,'row index');

%print(handles.figure1,'-depsc2','lgn_est5.eps');

