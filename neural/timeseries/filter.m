srate = 1000;
t = 0:1/srate:3;
n = length(t);

% create data, comprising a double-sine and noise
signal = 2*sin(2*pi*.8*t) + sin(2*pi*6*t);
noise  = 100*sin(2*pi*50*t);
data   = signal + noise;

% FFT of data (with Hann taper)
% Try this with and without the Hann taper to get a sense of the effect of
% tapering the data. You can easily observe the balance between maximizing
% the signal while introducing edge artifacts, vs. minimizing edge
% artifacts while also attenuating part of the real signal.
hann  = .5*(1-cos(2*pi*(1:n)/(n-1)));
dataX = fft(data.*hann);
hz    = linspace(0,srate,length(t));

% create low-pass filter and apply to data spectrum
filterkernel = (1-1./(1+exp(-hz+40)));
dataX1       = dataX.*filterkernel;

% compute time-domain version of filtered signal
data2 = real(2*ifft(dataX1)); % *2 b/c one-sided filter


figure(1), clf
subplot(311)
plot(t,data)
ylabel('Amplitude')
title('Data (signal + noise)')
set(gca,'ylim',[-120 120])

subplot(312)
plot(t,data2), hold on
plot(t,signal,'r')
ylabel('Amplitude'), xlabel('Time (s)')
title('With Hann taper')
legend({'filtered';'signal without noise'})
set(gca,'ylim',[-5 5])



dataX  = fft(data);
dataX1 = dataX.*filterkernel;
data2  = real(2*ifft(dataX1)); % *2 b/c one-sided filter
subplot(313)
plot(t,data2), hold on
plot(t,signal,'r')
ylabel('Amplitude'), xlabel('Time (s)')
title('Without Hann taper')
legend({'filtered';'signal without noise'})
set(gca,'ylim',[-5 5])



figure(2), clf
subplot(211)
plot(hz,2*abs(dataX/n))
set(gca,'xlim',[0 70])
ylabel('Amplitude')

subplot(212)
plot(hz,2*abs(dataX1/n))
set(gca,'xlim',[0 70])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
% Note the difference in amplitude between the two plots. Scaling the top
% plot to set(gca,'ylim',[0 1]) will reveal the same signal frequency
% components that can be seen in the lower plot.
% This time add random noise to the signal.
% Try changing the line below from random Gaussian noise to random uniform
% noise (function rand). Why does the moving-average filter over-estimate
% the signal, and what does this tell you about when it is appropriate to
% use a moving-average filter?
noise = randn(size(signal));
data  = signal + noise;

% moving-average filter
d=9; % this really means a 19-point mean filter!
dataMean=zeros(size(data));
for i=d+1:length(t)-d-1
    dataMean(i) = mean(data(i-d:i+d));
end

clf
subplot(211)
plot(t,data)
ylabel('Amplitude')
title('Signal and noise')

subplot(212), hold on
plot(t,dataMean,'b')
plot(t,signal,'r')
xlabel('Time (s)'), ylabel('Amplitude')
title('Signal and noise, filtered')
legend({'filtered';'original signal'})

%% Chapter 10.3, Figure 10.3

% Again, data is signal plus noise
noise = randn(size(signal));
data  = signal + noise;
d=9; % filter parameter (19-points)

% Gaussian-weighted filter, using a width of 2
gausfilt = exp(-(-d:d).^2/4);
halfGausSize = floor(length(gausfilt)/2);

% convolution as short-cut to moving-average weighted filter
Lconv    = length(gausfilt)+length(t)-1;
convres  = ifft( fft(data,Lconv).*fft(gausfilt,Lconv) );
dataGaus = convres(halfGausSize:end-halfGausSize-1);

clf
% time domain data
subplot(211)
plot(t,data)
ylabel('Amplitude')

subplot(212), hold on
% note the amplitude scaling according to the height of the Gaussian. 
% This would be avoided by amplitude-scaling the Gaussian so that its
% integral (sum over all points) is zero.
plot(t,dataGaus/sum(gausfilt),'b')
plot(t,signal,'r')
xlabel('Time (s)'), ylabel('Amplitude')
legend({'Filtered';'original'})

%% Chapter 10.4, Figure 10.4

% Now the noise is created differently from previous examples.
% Each index correspond to a prime number is set to 100. This 
% creates large outliers that are always positive.
noise = zeros(size(t));
% necessity of 'find' is version-dependent
noise(isprime(1:length(t))) = 100;
data = signal + noise;

% try mean and median filter
d = 9;
[dataMed,dataMean]=deal(zeros(size(data)));
for i=d+1:length(t)-d-1
  dataMed(i)  = median(data(i-d:i+d));
  dataMean(i) = mean(data(i-d:i+d));
end

clf
% Plot the data. Note that y-axis scaling; the signal itself
% has values between -3 and +3, but the noise goes up to 100.
subplot(211)
plot(t,data)
set(gca,'ylim',[-10 110])
ylabel('Amplitude')

subplot(212)
plot(t,dataMed), hold on
plot(t,dataMean,'r')
plot(t,signal,'k')
set(gca,'ylim',[-5 40])
xlabel('Time (s)'), ylabel('Amplitude')
legend({'median filter';'mean filter';'original'})

%% Chapter 10.4

% Running-median filter using an order>1
n_order = 6;
% d parameter increases with each order
d = round(linspace(5,19,n_order));
dataMedFull = zeros(n_order,length(t));

% loop through order and loop through time points
for oi=1:n_order
    for i=d(oi)+1:length(t)-d(oi)-1
        temp = sort(data(i-d(oi):i+d(oi)));
        dataMedFull(oi,i) = temp(floor(length(temp)/2)+1);
    end
end

% The median was computed 6 times; now average the results 
% from the different orders together.
dataMed = mean(dataMedFull);

clf
plot(t,dataMed,'m')
hold on
plot(t,signal,'r')
xlabel('Time (s)'), ylabel('Amplitude')
legend({'median filter';'original'})

%% Chapter 10.5, Figure 10.5

% create signal
signal = .5*sin(2*pi*60*t) + sin(2*pi*6*t);

% add sharp noise (prime numbers, as used earlier)
noise = zeros(size(t));
% necessity of 'find' is version-dependent
noise(isprime(1:length(t))) = 100;
data = signal + noise;
d=9; % this really means a 19-point median filter!
dataMed=data;
points2filter = find(data>2*std(data)+median(data));

for i=1:length(points2filter)
    centpoint = points2filter(i);
    dataMed(centpoint) = median(data(max(1,centpoint-d):min(length(data),centpoint+d)));
end

clf
subplot(211)
plot(t,data)
set(gca,'ylim',[-10 110])

subplot(212)
plot(t,dataMed), hold on
plot(t,signal,'k')
set(gca,'xlim',[0 .5])

%% 10.5

% The code below shows the performance of the non-threshold-based filter. 
% Although the large peaks are removed, so are the higher-frequency fluctuations,
% which are part of the signal.

dataMed2 = zeros(size(data));
for i=d+1:length(t)-d-1
    dataMed2(i)  = median(data(i-d:i+d));
end
plot(t,dataMed2,'r')

%% Chapter 10.6, Figure 10.7

% x a linear function
x = 1:20;
% and y is a function of x plus noise
y = 2*x+randn(size(x));

% Matlab function to fit a polynomial (1 is for first-order)
p = polyfit(x,y,1);

clf
plot(x,y,'o-'), hold on
plot(x,p(2)+p(1)*x,'r*-')
xlabel('x'), ylabel('y')
legend({'original';'denoised'})

%% Chapter 10.6, Figure 10.8

% create signal and add noise
srate  = 1000;
t      = 0:1/srate:5;
signal = interp1(0:5,randn(6,1),t,'spline');
noise  = 3*randn(size(t));
data   = signal+noise;

% fit polynomial and extract fitted response
polyorder   = 6;
p           = polyfit(t,data,polyorder);
dataPolyFit = polyval(p,t);

clf
subplot(211)
plot(t,data)
ylabel('Amplitude')

subplot(212)
plot(t,dataPolyFit)
hold on
plot(t,signal,'r')
xlabel('Time (s)'), ylabel('Amplitude')
legend({'denoised';'original'})

%% Chapter 10.6

% In the previous cell, the 'noise' had high-frequency characteristics
% and the 'signal' had low-frequency characteristics. If the opposite is
% the case (the low-frequency component is the noise), the data can be
% subtracted from the polynomial fit.

clf
subplot(211)
plot(t,data)

subplot(212)
plot(t,data-dataPolyFit,'k')

%% Chapter 10.7, Figure 10.9

% Matlab and Octave come with several built-in images and datasets.

clf

% soft-code color limit
clim = [3 30];  % for Matlab
% clim = [3 300]; % for Octave

% get signal
image
signal = get(findobj(gcf,'type','image'),'CData');

% data is signal + noise
data = signal + 1000*reshape(isprime(1:numel(signal)),size(signal));

% d (distribution) parameter 
d = 1;

% threshold for filtering is 2 standard deviations above 
% the median of the data.
thresh=2*std(data(:))+median(data(:));

% loop through two dimensions and apply filter 
% if the threshold is exceeded.
dataMed = data;
for i=d+1:size(data,1)-d-1
    for j=d+1:size(data,2)-d-1
        if data(i,j)>thresh
            temp = data(i-d:i+d,j-d:j+d);
            dataMed(i,j) = median(temp(:));
        end
    end
end

% show original image
subplot(131)
imagesc(signal)
set(gca,'clim',clim)
axis image

% show currupted data
subplot(132)
imagesc(data)
set(gca,'clim',clim)
axis image

% and cleaned version. note the residual artifacts around the edges.
subplot(133)
imagesc(dataMed)
set(gca,'clim',clim)
axis image

colormap gray

%% exercise 1

clf

% get image as above
image
signal = get(findobj(gcf,'type','image'),'CData');
data   = signal + 1000*reshape(isprime(1:numel(signal)),size(signal));

% initialize variables. Now d is set to 2. Try setting it to 1 to see what
% happens. What is the difference, and why does it happen?
d       = 2;
thresh  = 2*std(data(:))+median(data(:));
dataMed = data;
dims    = size(data);

% loop through and filter
for i=1:size(data,1)
    for j=1:size(data,2)
        if data(i,j)>thresh
            % this line is the key difference with the previous cell
            temp = data(max(1,i-d):min(dims(1),i+d),max(1,j-d):min(dims(2),j+d));
            dataMed(i,j) = median(temp(:));
        end
    end
end

% plot original
subplot(131)
imagesc(signal)
set(gca,'clim',clim)
axis image

% plot image with noise
subplot(132)
imagesc(data)
set(gca,'clim',clim)
axis image

% filtered image
subplot(133)
imagesc(dataMed)
set(gca,'clim',clim)
axis image

colormap gray

%% end
