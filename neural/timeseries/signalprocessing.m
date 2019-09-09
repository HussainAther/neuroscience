

%% Chapter 8.1, Figure 8.1

clf

% signal properties
srate = 1000;
t = 0:1/srate:5;

% loop over three different time courses
for i=1:3
    switch i
        case 1
            % signal 1 is linearly increasing
            freqTS = linspace(1,20,length(t));
        case 2
            % signal 2 is a smoothed random sequence
            freqTS = abs(interp1(linspace(t(1),t(end),10),10*rand(1,10),t,'spline'));
        case 3
            % signal 3 is a triangle
            freqTS = abs(mod(t,2)-1)*10;
    end
    
    % These signals are then used as frequency time series.
    % Below is the general formula to create a signal with 
    %  arbitrary changes in instantaneous frequency.
    centfreq = mean(freqTS);
    k = (centfreq/srate)*2*pi/centfreq;
    y = sin(2*pi.*centfreq.*t + k*cumsum(freqTS-centfreq));
    
    % plot sinusoidal signal derived from time-varying instantaneous frequencies
    subplot(3,2,i*2)
    plot(t,y)
    set(gca,'ylim',[-1.1 1.1])
    if i==3
        % show axis labels only on bottom plots
        xlabel('Time (s)'), ylabel('Amplitude')
    end
    
    % plot the instantaneous frequencies
    subplot(3,2,(i-1)*2+1)
    plot(t,freqTS)
    if i==3
        % show axis labels only on bottom plots
        xlabel('Time (s)'), ylabel('Frequency (Hz)')
    end
    
end

%% Chapter 8.2, Figure 8.2

% signal properties
srate = 1000;
t  = 0:1/srate:5;
f  = [4 6];
ff = linspace(f(1),mean(f),length(t));

% create chirp signal
signal = sin(2*pi.*ff.*t);

% the phases are extracted from the Hilbert transform. This function is in
% the Signal processing toolbox in Matlab, the Signal package in Octave, or
% it can be created using instructions in Chapter 7.
phases = angle(hilbert(signal));

% Angular velocity is defined as the first temporal derivative 
% of the unwrapped phase angles.
angVeloc = diff(unwrap(phases));

% convert the results to Hz...
instFreq1 = srate*angVeloc/(2*pi);

subplot(211)
plot(t,signal)
ylabel('Amplitude')

subplot(212)
% This next line will crash if you try to plot the entire time vector.
% Why is this, and what is a solution to this?
plot(t(1:end-1),instFreq1)
xlabel('Time (s)'), ylabel('Frequencies (Hz)')

%% Chapter 8.3, Figure 8.3

% Properties of background activity that will be added to the 
% instantaneous frequency time series created above.
a = [10 2 5 8]/10;
f = [3 1 6 12];

background = zeros(size(t));
for i=1:length(a)
    background = background + a(i)*sin(2*pi*f(i)*t);
end

% New data are created from signal plus background activity.
data = signal + background;

% Compute instantaneous frequencies.
instFreq2 = srate*diff(unwrap(angle(hilbert(data))))/(2*pi);

subplot(211)
plot(t,data)
ylabel('Amplitude')

subplot(212)
plot(t(1:end-1),instFreq2)
xlabel('Time (s)'), ylabel('Frequencies (Hz)')

%% Chapter 8.3, Figure 8.4

clf
plot(t,angle(hilbert(data)))
set(gca,'ytick',-pi:pi/2:pi,'ylim',[-3.3 3.3])
xlabel('Time (s)'), ylabel('Phase angles (rad.)')

%% Chapter 8.3, Figure 8.5

% create complex Morlet wavelet with a peak frequency of 5 Hz
wavetime = -2:1/srate:2;
wavefreq = 5; % Hz
nwavecyc = 5; % number of wavelet cycles
w        = 2*( nwavecyc/(2*pi*wavefreq) )^2;
cmw      = exp(1i*2*pi*wavefreq.*wavetime) .* exp( (-wavetime.^2)/w );

% convolution parameters
halfwavsize = floor(length(wavetime)/2);
Lconv = length(t)+length(wavetime)-1;

% convolution
convres = ifft( fft(data,Lconv).*fft(cmw,Lconv) );
convres = convres(halfwavsize:end-halfwavsize-1);

% compute instantaneous frequencies again
instFreq3 = srate*diff(unwrap(angle(convres)))/(2*pi);

% add additional data point to the end for convenience in plotting
instFreq1(numel(t)) = instFreq1(numel(t)-1);
instFreq2(numel(t)) = instFreq2(numel(t)-1);
instFreq3(numel(t)) = instFreq3(numel(t)-1);

clf
plot(t,instFreq1), hold on,
plot(t,instFreq3,'r')
plot(t,instFreq2,'m')
set(gca,'ylim',[3 7])
xlabel('Time (s)'), ylabel('Frequencies (Hz)')
legend({'Original';'With background, filtered';'With background, unfiltered'})

% Hint: try changing nwavecyc from 5 to 12. Then try changing it to 2. 
% What do the results tell you about the appropriate approach to filtering
% when isolating band-limited instantaneous frequencies in the presence of 
% background activity?

%% Chapter 8.4, Figure 8.6

% add noise to signal
data     = signal + randn(size(signal));
phases   = angle(hilbert(data));
angVeloc = diff(unwrap(phases));

% compute instantaneous frequency
instFreq = srate*angVeloc/(2*pi);
instFreq(end+1) = instFreq(end);

clf
% plot the time domain signal
subplot(411)
plot(t,data)
ylabel('Amplitude')

% plot the phase angles
subplot(412)
plot(t,phases)
set(gca,'ytick',-pi:pi/2:pi,'ylim',[-3.3 3.3])
ylabel('Phases (rad.)')

% plot instantaneous frequencies
subplot(413)
plot(t,instFreq)
ylabel('Frequencies (Hz)')

%% Chapter 8.4, Figure 8.7

% convolution with the same wavelet defined earlier
convres = ifft( fft(data,Lconv).*fft(cmw,Lconv) );
convres = convres(halfwavsize:end-halfwavsize-1);

% recompute instantaneous frequencies
instFreq3 = diff(unwrap(angle(convres)));
instFreq3(length(t)) = instFreq3(end);

% and plot
subplot(414)
plot(t,instFreq1), hold on,
plot(t,srate*instFreq3/(2*pi),'r')
set(gca,'ylim',[3 7])
legend({'original';'With noise, filtered'})
xlabel('Time (s)'), ylabel('Frequency (Hz)')

%% Chapter 8.5, Figure 8.8

srate = 1000;
t = 0:1/srate:5;
signal = zeros(size(t));

% pairs of frequencies to use
fr = [13 16; 24.5 24; 32 39; 50 48];

for i=1:4
    ff = linspace(fr(i,1),fr(i,2)*mean(fr(i,:))/fr(i,2),length(t));
    signal = signal + sin(2*pi.*ff.*t);
end

hz = linspace(0,srate/2,floor(length(t)/2)+1);
x  = fft(signal)/length(t);

clf
subplot(211)
plot(t,signal)
xlabel('Time (s)'), ylabel('Amplitude')

subplot(212)
plot(hz,2*abs(x(1:length(hz))))
set(gca,'xlim',[0 60])
xlabel('Frequency (Hz)'), ylabel('Amplitude')

%% Chapter 8.5, Figure 8.9

% Compute intrinsic mode functions via empirical mode decomposition.
% The function 'emdx' is included in the online Matlab code and can be
% inspected.
imfs = emdx(signal,4); % return first 4 modes

% compute Fourier spectra of all modes
imfsX = fft(imfs,[],2)/length(t);

clf
% plot power
plot(hz,abs(imfsX(:,1:length(hz))).^2)
set(gca,'xlim',[0 60])
xlabel('Frequency (Hz)'), ylabel('Power')

%% Chapter 8.6, exercise 1, Figure 8.11

srate = 1000;
t = 0:1/srate:6;

% create double-chirp by summing two chirps
freqTS   = linspace(10,13,length(t));
centfreq = mean(freqTS);
k        = (centfreq/srate)*2*pi/centfreq;
signal   = sin(2*pi.*centfreq.*t + k*cumsum(freqTS-centfreq));

freqTS   = linspace(35,30,length(t));
centfreq = mean(freqTS);
k        = (centfreq/srate)*2*pi/centfreq;
signal   = signal + sin(2*pi.*centfreq.*t + k*cumsum(freqTS-centfreq));


% instantaneous frequency of broadband signal
if_broadband = srate*diff(unwrap(angle(hilbert(signal))))/(2*pi);
if_broadband(length(t)) = if_broadband(end);

% plot instantaneous frequencies from broadband signal
clf
subplot(221), plot(t,if_broadband)
set(gca,'ylim',[5 40])


% filter at 11.5 and 32.5 Hz and compute instantaneous frequency again
wavetime  = -2:1/srate:2;
Lconv     = length(t)+length(wavetime)-1;
halfwavsize= floor(length(wavetime)/2);
frex = [11.5 32.5];

sigX     = fft(signal,Lconv);
instfrex = zeros(length(frex),length(t)-1);

for fi=1:2
    w       = 2*( 10/(2*pi*frex(fi)) )^2; % w is width of Gaussian
    cmwX    = fft(exp(1i*2*pi*frex(fi).*wavetime) .* exp( (-wavetime.^2)/w ), Lconv);
    cmwX    = cmwX./max(cmwX);
    convres = ifft( sigX.*cmwX );
    convres = convres(halfwavsize:end-halfwavsize-1);
    
    instfrex(fi,:) = srate*diff(unwrap(angle(convres)))/(2*pi);
end
instfrex(:,length(t)) = instfrex(:,end);

subplot(222)
plot(t,instfrex)
set(gca,'ylim',[5 40])



% TF analysis with Morlet wavelet convolution
nfrex = 50;
frex  = linspace(5,40,nfrex);
tf    = zeros(length(frex),length(t));
n     = linspace(8,30,nfrex);

for fi=1:nfrex
    w       = 2*( n(fi)/(2*pi*frex(fi)) )^2; % w is width of Gaussian
    cmwX    = fft(exp(1i*2*pi*frex(fi).*wavetime) .* exp( (-wavetime.^2)/w ), Lconv);
    cmwX    = cmwX./max(cmwX);
    convres = ifft( sigX.*cmwX );
    convres = convres(halfwavsize:end-halfwavsize-1);
    
    tf(fi,:) = abs(convres).^2;
end

subplot(223)
contourf(t,frex,tf,40,'linecolor','none')


subplot(224)
freqrange  = dsearchn(frex',[25 50]');
[vals,idx] = max(tf(freqrange(1):freqrange(2),:));
plot(t,frex(idx+freqrange(1)+1))
hold on

freqrange  = dsearchn(frex',[4 25]');
[vals,idx] = max(tf(freqrange(1):freqrange(2),:));
plot(t,frex(idx+freqrange(1)+1))
set(gca,'ylim',[5 40])

%% end
