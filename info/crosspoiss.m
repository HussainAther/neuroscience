% Generate correlated Poisson spike trains

rho = 0.2; %events/ms, 200 events per second

%probability of keeping a spike
prob_keep = 0.1;

dt = 0.5; %in ms
fs = 1/(dt*1e-3); %sampling frequency in Hz


t_max = 200; %ms
nev = 4*t_max*rho;
isi_v = exprnd(1/rho,1,nev);
spk_v = cumsum(isi_v);
spk_v = spk_v(find(spk_v<t_max));

%initializes the data structure
n_trains = 10;
for i = 1:n_trains
    s(i).tspk = [];
end;

for i = 1:length(spk_v)
    u = find(rand(1,n_trains)<=prob_keep);
    for j = 1:length(u)
        s(u(j)).tspk = [s(u(j)).tspk spk_v(i)];
    end;
end;

hf = figure; 
ha = axes;
for i = 1:length(spk_v)
    line('Parent',ha,'XData',[spk_v(i) spk_v(i)],'YData',[(n_trains+1.05) (n_trains+1.95)],...
        'Color','r');
end;

for i = 1:n_trains
  for j = 1:length(s(i).tspk)
    line('Parent',ha,'XData',[s(i).tspk(j) s(i).tspk(j)],'YData',[i+0.05 (i+0.95)]);
  end;
end;


%uncomment to generate associated EPSP conductance changes and plot them 
% Nt = ceil(t_max/dt);
% 
% %excitatory epsp time constant
% tauE = 1; %ms
%  
% epsp_v = zeros(n_trains,Nt);
%  
% for j = 1:n_trains
%     for i = 1:Nt
%         tc = i*dt;
%         epsp_v(j,i) = 0.5*sum( ( max( (tc-s(j).tspk),0 )/tauE).*exp( 1 - (max( (tc-s(j).tspk),0)/tauE) ) );
%     end;
% 
%     line('Parent',ha,'XData', (1:Nt)*dt, 'YData', epsp_v(j,:) + j,'Color','k','LineStyle','-');
% end;

set(ha,'TickDir','out');
xlabel(ha,'time (ms)');

%uncomment to save figure
%print(hf,'-depsc2','cross_poiss1.eps');

%compute the crosscorrelation between two such spike trains

t_max = 32768; %ms, 2^15
nev = ceil(4*t_max*rho);
isi_v = exprnd(1/rho,1,nev);
spk_v = cumsum(isi_v);
spk_v = spk_v(find(spk_v<t_max));

%initializes the data structure
n_trains = 2;
for i = 1:n_trains
    s(i).tspk = [];
end;

for i = 1:length(spk_v)
    u = find(rand(1,n_trains)<=prob_keep);
    for j = 1:length(u)
        s(u(j)).tspk = [s(u(j)).tspk spk_v(i)];
    end;
end;

%uncomment to plot spike trains
%figure; subplot(n_trains+1,1,1);
%stem(spk_v,ones(size(spk_v)),'Marker','none');

%for i = 1:n_trains
%    subplot(n_trains+1,1,i+1);
%    stem(s(i).tspk,ones(size(s(i).tspk)),'Marker','none');
%end;

%generates spike trains
Nt = ceil(t_max/dt);

spk = zeros(2,Nt);
inds1 = ceil(s(1).tspk/dt);
spk(1,inds1) = 1;
inds2 = ceil(s(2).tspk/dt);
spk(2,inds2) = 1;

%converts to delta function approximation (in spk/s)
spk = spk*fs;

m_spk1 = mean(spk(1,:));
m_spk2 = mean(spk(2,:));
s_spk1 = std(spk(1,:));
s_spk2 = std(spk(2,:));

%zero mean
spk(1,:) = spk(1,:)-m_spk1;
spk(2,:) = spk(2,:)-m_spk2;

window = 2048;
noverlap = 1024;
nfft = 2048;
[Pxy,F] = cpsd(spk(1,:),spk(2,:),window,noverlap,nfft,fs,'twosided');
%figure; plot(F,abs(Pxy));

%compute the inverse Fourier transform, but multiply by nfft, since we
%did not start with the original sequence
Cxy = nfft*ifft(Pxy);
dim = find(size(Cxy) ~= 1);
if ( dim == 1 )
    shift_v = [(nfft/2)-1 0];
else 
    shift_v = [0 (nfft/2)-1];
end;

%correct for the phase shift
Cxy = circshift(Cxy,shift_v);
t_v = [-( (nfft/2)-1 ):(nfft/2)]*dt;
%figure; plot(t_v,Cxy);

Cxyn = Cxy/(std(spk(1,:))*std(spk(2,:)));
hf = figure; 
ha = axes;
plot(t_v,Cxyn,'k');
set(ha,'TickDir','out');
xlabel(ha,'time (ms)');
ylabel(ha,'normalized crosscorrelation');

%uncomment to save figure
%print(hf,'-depsc2','cross_poiss2.eps');

