

t = -2:.005:2; 

% Morlet wavelet
f=4;
csw = cos(2*pi*f*t);
w = 2*( 5/(2*pi*f) )^2;
gaussian = exp( (-t.^2)/w );
MorletWavelet = csw .* gaussian;

subplot(311)
plot(t,MorletWavelet)
set(gca,'ylim',[-1 1.1])


% Haar wavelet
HaarWavelet = zeros(size(t));
HaarWavelet(dsearchn(t',0):dsearchn(t',.5)) = 1;
HaarWavelet(dsearchn(t',.5):dsearchn(t',1)) = -1;

subplot(312)
plot(t,HaarWavelet)
set(gca,'ylim',[-1.1 1.1])



% Mexican hat wavelet
s = .4;
MexicanWavelet = (2/(sqrt(3*s)*pi^.25)) .* (1- (t.^2)/(s^2) ) .* exp( (-t.^2)./(2*s^2) );
subplot(313)
plot(t,MexicanWavelet)
set(gca,'ylim',[-.7 1.5])


% create complex sine wave
csw = exp(1i*2*pi*f*t);

% Morlet wavelet is a sine wave tapered by a Gaussian.
% The variable gaussian was created above.
mwavelet = csw .* gaussian;

clf
plot3(t,real(mwavelet),imag(mwavelet),'k-o','markerfacecolor','m');
xlabel('Time (s)'), ylabel('real part'),
zlabel('imaginary part')
rotate3d % unnecessary in Octave


% create a linear chirp
srate   = 1000;
t       = 0:1/srate:6;
f       = [2 8]; % frequencies in Hz
chirpTS = sin(2*pi.*linspace(f(1),mean(f),length(t)).*t);

% The wavelet has its own time vector, but 
% the sampling rate MUST be the same as for the signal.
wavetime = -2:1/srate:2;

% width of the Gaussian that tapers the sine wave
w = 2*( 4/(2*pi*5) )^2;

% create the complex Morlet wavelet
cmw = exp(1i*2*pi*5.*wavetime) .* exp( (-wavetime.^2)/w );

% half of the length of the wavelet
halfwavsize = floor(length(wavetime)/2);

% Perform time-domain convolution...
% the signal must be zero-padded
chirpTSpad = [zeros(1,halfwavsize) chirpTS zeros(1,halfwavsize)];
% initialize results matrix
convres = zeros(size(chirpTSpad));

% and now the actual convolution
for i=halfwavsize+1:length(chirpTS)+halfwavsize-1
    convres(i) = sum( chirpTSpad(i-halfwavsize:i+halfwavsize) .* cmw );
end
% the final step of convolution is to trim the edges
convres = convres(halfwavsize:end-halfwavsize-1);

% and plot the results
clf
subplot(211)
plot(t,chirpTS)

subplot(212)
plot(t,abs(convres))
xlabel('Time (s)'), ylabel('Amplitude (unscaled)')

% length of convolution
Lconv = length(t)+length(wavetime)-1;

% Convolution as multiplication in the frequency domain
convres2 = ifft( fft(chirpTS,Lconv).*fft(cmw,Lconv) );
convres2 = convres2(halfwavsize:end-halfwavsize-1);

hold on, plot(t,abs(convres2),'r')
Lconv = length(t)+length(wavetime)-1;

% Fourier spectra of the two signals
cmwX = fft(cmw,Lconv);
cmwX = cmwX./max(cmwX);

% and their multiplied inverse
convres4 = ifft( fft(chirpTS,Lconv).*cmwX );
convres4 = convres4(halfwavsize:end-halfwavsize-1);

clf
plot(t,chirpTS), hold on
plot(t,2*abs(convres4),'r')
% number and range of frequencies for analysis
nfrex = 30;
frex  = logspace(log10(2),log10(30),nfrex);

% initialize...
tf = zeros(nfrex,length(chirpTS));

% Fourier spectrum of the signal. Note that because this does not change as
% a function of wavelet frequency, it needs to be computed only once.
chirpTSx = fft(chirpTS,Lconv);

% loop through frequencies
for fi=1:nfrex
    % compute normalized Fourier spectra of wavelet
    w = 2*( 5/(2*pi*frex(fi)) )^2;
    cmwX = fft(exp(1i*2*pi*frex(fi).*wavetime) .* exp( (-wavetime.^2)/w ), Lconv);
    cmwX = cmwX./max(cmwX);
    
    % convolution and extract amplitude
    convres = ifft( chirpTSx .* cmwX );
    tf(fi,:) = 2*abs(convres(halfwavsize:end-halfwavsize-1));
end

clf
contourf(t,frex,tf,40,'linecolor','none')
set(gca,'clim',[0 1]), colorbar
times2plot = dsearchn(t',(1:0.2:4.5)');
freq2plot  = dsearchn(frex',6);

% full data
plot(t,tf(freq2plot,:)), hold on

% downsampled data
plot(t(times2plot),tf(freq2plot,times2plot),'ro-')

clf

% properties of wavelet with 7 cycles
freq2use= 5;
w       = 2*( 7/(2*pi*freq2use) )^2;
mwave7  = exp(1i*2*pi*freq2use.*wavetime) .* exp( (-wavetime.^2)/w );
mwave7x = fft(mwave7,Lconv);
mwave7x = mwave7x./max(mwave7x);
convres = ifft( chirpTSx .* mwave7x );
pow7    = 2*abs(convres(halfwavsize:end-halfwavsize-1));

% same as previous, but using 3 cycles
w       = 2*( 3/(2*pi*freq2use) )^2;
mwave3  = exp(1i*2*pi*freq2use.*wavetime) .* exp( (-wavetime.^2)/w );
mwave3x = fft(mwave3,Lconv);
mwave3x = mwave3x./max(mwave3x);
convres = ifft( chirpTSx .* mwave3x );
pow3    = 2*abs(convres(halfwavsize:end-halfwavsize-1));


subplot(221)
plot(wavetime,real(mwave7)), hold on
plot(wavetime,real(mwave3),'r')
xlabel('Time (s)'), ylabel('Amplitude')
title('Projection of wavelets onto time and real axes')

subplot(222)
hz = linspace(0,srate/2,floor(length(wavetime)/2+1));
x7 = 2*abs(fft(mwave7)); 
x3 = 2*abs(fft(mwave3));
plot(hz,x7(1:length(hz))), hold on
plot(hz,x3(1:length(hz)),'r')
set(gca,'xlim',[0 20])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Non-normalized power spectra of wavelets')
legend({'7 cycles';'3 cycles'})

subplot(212)
plot(t,pow7), hold on
plot(t,pow3,'r')
xlabel('Time (s)'), ylabel('Amplitude')
title([ num2str(freq2use) ' Hz component of chirp' ])
legend({'7 cycles';'3 cycles'})


nfrex = 9;
frex  = logspace(log10(2),log10(30),nfrex);
% n is the parameter that sets the number of cycles
ncyc = logspace(log10(3),log10(12),nfrex);

% create family of wavelets
wavelet_fam = zeros(nfrex,length(wavetime));
for fi=1:nfrex
    % all wavelets in this family are the same except 
    % for the 'n' parameter
    w = 2*( ncyc(fi)/(2*pi*frex(fi)) )^2;
    wavelet_fam(fi,:) = exp(1i*2*pi*frex(fi).*wavetime) ...
        .* exp( (-wavetime.^2)/w );
end

for i=1:9
    subplot(3,3,i)
    plot(wavetime,real(wavelet_fam(i,:)))
    set(gca,'xlim',[-1 1],'ylim',[-1.1 1.1])
end

srate=1000; 
t = 0:1/srate:12;
f = [6 14 25 40 70];
% note the widely varying amplitudes
a = [.001234 1.234 1234 123400 12340000];
% relevant amplitude modulations
m = [-.1 .1 -.2 .2 -.3];

basetidx = dsearchn(t',[2 4]');

signal = zeros(size(t));
for i=1:length(f)
    
    % compute 'base' signal
    signal = signal + a(i).*sin(2*pi*f(i).*t);
    
    % add time-limited modulation to signal
    extrasignal = m(i)*a(i)*sin(2*pi*f(i).*t) .* exp( -(t-7).^2 );
    signal = signal + extrasignal;
end

clf

% plot the time-domain signal
subplot(211)
plot(t,signal)
xlabel('Time (s)'), ylabel('Amplitude')



% setup for wavelet convolution
nfrex       = length(f);
wavetime    = -2:1/srate:2;
Lconv       = length(t)+length(wavetime)-1;
halfwavsize = floor(length(wavetime)/2);
tf          = zeros(nfrex,length(t));
signalX     = fft(signal,Lconv);

% loop over frequencies and perform convolution
for fi=1:nfrex
    w = 2*( 10/(2*pi*f(fi)) )^2;
    cmwX = fft(exp(1i*2*pi*f(fi).*wavetime) .* exp( (-wavetime.^2)/w ), Lconv);
    cmwX = cmwX./max(cmwX);
    convres = ifft( signalX .* cmwX );
    tf(fi,:) = 2*abs(convres(halfwavsize:end-halfwavsize-1));
end

% plot the non-dB-normalized result
subplot(223)
plot(t,tf)
xlabel('Time (s)'), ylabel('Amplitude (raw)')
title('Non-normalized amplitude')

% convert to dB and plot again
subplot(224)
db = 10*log10( bsxfun(@rdivide,tf,mean(tf(:,basetidx(1):basetidx(2)),2)) );
plot(t,db)
set(gca,'ylim',[-2 2])
xlabel('Time (s)'), ylabel('Amplitude (dB)')
title('dB-normalized amplitude')
% The signal will have two overlapping sine waves,
% one with increasing power and one with decreasing power.

% first, the basics:
srate     = 1000;
time      = 0:1/srate:6;
sinefreqs = [5 11]; % in hz
numTrials = 20;
wavetime  = -2:1/srate:2;
nfrex     = 9; 
frex      = logspace(log10(2),log10(30),nfrex);
ncyc      = logspace(log10(3),log10(12),nfrex);
Lconv     = length(time)+length(wavetime)-1;
halfwavsize= floor(length(wavetime)/2);

% second, initialize matrices
signal = zeros(numTrials,length(time));
tf     = zeros(2,nfrex,length(time));

% third, create the signal
ampInc = linspace(0,1,length(time));
ampDec = linspace(1,0,length(time));

for ti=1:numTrials
    sine1 = ampInc.*sin(2*pi*sinefreqs(1)*time + rand*2*pi);
    sine2 = ampDec.*sin(2*pi*sinefreqs(2)*time + rand*2*pi);
    signal(ti,:) =  sine1+sine2;
end

% fourth, perform time-frequency decomposition
dataX1 = fft(signal,Lconv,2);     % FFT along 2nd dimension
dataX2 = fft(mean(signal),Lconv); % note that trials are first averaged

for fi=1:nfrex
    % create wavelet
    w          = 2*( ncyc(fi)/(2*pi*frex(fi)) )^2; % w is width of Gaussian
    cmwX  = fft(exp(1i*2*pi*frex(fi).*wavetime) .* exp( (-wavetime.^2)/w ), Lconv);
    cmwX  = cmwX./max(cmwX);
    
    % convolution of all trials simultaneously using bxsfun
    convres    = ifft( bsxfun(@times,dataX1,cmwX) ,[],2);
    temppower  = 2*abs(convres(:,halfwavsize:end-halfwavsize-1));
    tf(1,fi,:) = mean(temppower,1);
    
    % The 2 lines above are equivalent to the 5 lines below, which uses
    % loops and thus should be avoided when possible.
%     for ti=1:numTrials
%         dataX = fft(signal(ti,:),Lconv);
%         convres = ifft( dataX .* mwave_fft );
%         tf(ti,fi,:) = 2*abs(convres(halfwavsize:end-halfwavsize-1));
%     end

    
    % convolution of trial average
    convres    = ifft( dataX2.*cmwX );
    tf(2,fi,:) = 2*abs(convres(:,halfwavsize:end-halfwavsize-1));
end

% finally (sixth), plot the average and TF power
subplot(221)
plot(time,signal)
set(gca,'ylim',[-1.5 1.5])

subplot(223)
contourf(time,frex,squeeze(tf(1,:,:)),40,'linecolor','none')
set(gca,'clim',[0 .6])
colormap gray


% finally (sixth), plot the average and TF power
subplot(222)
plot(time,mean(signal))
set(gca,'ylim',[-1.5 1.5])

subplot(224)
contourf(time,frex,squeeze(tf(2,:,:)),40,'linecolor','none')
set(gca,'clim',[0 .6])
colormap gray

%% end
