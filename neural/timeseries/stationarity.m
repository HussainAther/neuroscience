n = 10000;
nbins = 200;

x = linspace(0,5,n) + randn(1,n);
y = randn(1,n);

% This is a fast way to compute means over non-overlapping windows.
% Note that it will work only if 'nbins' neatly fits into 'n'.
% As an exercise, modify the code below so it will work if n=10001.
timeMeanX = mean(reshape(x,n/nbins,nbins));
timeMeanY = mean(reshape(y,n/nbins,nbins));

subplot(221)
plot(x)
xlabel('Time'), ylabel('Value')
title('Mean non-stationary signal')

subplot(222)
plot(y)
xlabel('Time'), ylabel('Value')
title('Mean stationary signal')

subplot(223)
plot(1:n/nbins:n,timeMeanX)
xlabel('Time'), ylabel('Mean value')

subplot(224)
plot(1:n/nbins:n,timeMeanY)
xlabel('Time'), ylabel('Mean value')

% compute the derivative of mean-non-stationary variable x
x = diff(x);
x(end+1)=x(end);

% recompute means over time
timeMeanX = mean(reshape(x,n/nbins,nbins));

clf
subplot(221)
plot(x)
xlabel('Time'), ylabel('Value')
title('Mean non-stationary signal, made stationary')

subplot(223)
plot(1:n/nbins:n,timeMeanX)
xlabel('Time'), ylabel('Value')

% increasing variance over time
x = linspace(1,5,n) .* randn(1,n);

% variance increases then decreases over time
y = exp(-linspace(-1,1,n).^2) .* randn(1,n);

% Same trick as above to compute variance over non-overlapping windows.
% This will also only work if 'nbins' neatly divides into 'n'.
timevarX = var(reshape(x,n/nbins,nbins),[],1);
timevarY = var(reshape(y,n/nbins,nbins),[],1);

subplot(221)
plot(x)
xlabel('Time'), ylabel('Value')
title('Variance non-stationary signal')

subplot(222)
plot(y)
xlabel('Time'), ylabel('Value')
title('Variance non-stationary signal')

subplot(223)
plot(1:n/nbins:n,timevarX)
xlabel('Time'), ylabel('Variance')

subplot(224)
plot(1:n/nbins:n,timevarY)
xlabel('Time'), ylabel('Value')

srate = 1000;
t = 0:1/srate:5;
n = length(t);

% Chirp of randomly fluctuating frequencies
freqTS = interp1(10*rand(5,1),linspace(1,5,n),'spline');
centfreq = mean(freqTS);
k = (centfreq/srate)*2*pi/centfreq;
y = sin(2*pi.*centfreq.*t + k*cumsum(freqTS-centfreq));

% compute phase velocity (that this value, when
% scaled, is the instantaneous frequency of the signal.
phases1 = diff(unwrap(angle(hilbert(y))));

% now create a signal that is phase stationary
y = sin(2*pi.*10.*t);
phases2 = diff(unwrap(angle(hilbert(y))));

% plot
subplot(211)
plot(t(1:end-1),phases1)
set(gca,'ylim',[0 .15])
ylabel('Phase angle derivative')
title('Phase non-stationary signal')

subplot(212)
plot(t(1:end-1),phases2)
set(gca,'ylim',[0 .15])
xlabel('Time (s)'), ylabel('Phase angle derivative')
title('Phase stationary signal')

% N per section
n = 100000;

% covariance matrices (3x3 because there will be 3 variables)
v{1} = [1 .5  0; .5  1  0;  0  0  1.1];
v{2} = [1.1 .1  0; .1  1.1  0;  0  1.1 .1];
v{3} = [.8  0 .9;  0 1  0; 0 0  1.1];

% initialize matrices
mvsig = zeros(0,length(v));
covar = cell(size(v));

% Loop through covariance segments
for i=1:length(v)
    % create ideal covariance matrix
    c = chol(v{i}*v{i}');
    
    % create random numbers with this covariance structure
    tempd = randn(n,size(v{i},1))*c;
    
    % compute the actual covariance
    covar{i} = (tempd'*tempd)/n;
    mvsig = cat(1,mvsig,tempd);
end
nn = length(mvsig);

% plot the three time series
figure(1), clf
subplot(211)
plot(mvsig)
hold on
plot([n n;n*2 n*2]',get(gca,'ylim'),'k:')

% and plot the covariance matrices within each segment
for i=1:3
    subplot(2,3,i+3)
    imagesc(covar{i})
    set(gca,'clim',[-1.5 1.5])
    title([ 'Covariance of segment ' num2str(i) ])
end
% Convert the signal to 3D as done previously 
% to easily compute mean and variance.
nbins = 200;
d3d   = reshape(mvsig,nn/nbins,nbins,length(v));

% mean stationarity
timeMean = reshape( mean(d3d,1),nbins,length(v) );
subplot(311)
plot(timeMean)
title('Channel means over time')

% variance stationarity
timeVar = reshape( var(d3d,[],1),nbins,length(v) );
subplot(312)
plot(timeVar)
title('Within-channel variance over time')


% initialize time-covariance matrices
timeCovar = zeros(nbins,length(v),length(v));
y = zeros(nbins,2);

% loop over bins
for bini=1:nbins
    
    % compute and store the covariance in this bin
    timeCovar(bini,:,:) = (squeeze(d3d(:,bini,:))'*squeeze(d3d(:,bini,:)))/size(d3d,1);
    
    % Compute 'distance' between the current and previous covariance
    % matrices, using two different distance metrics.
    if bini>1
        % First, Euclidean distance between the covariance matrices.
        y(bini,1) = sum( (reshape(timeCovar(bini,:,:)-timeCovar(bini-1,:,:),1,[]).^2 ));
        
        % Second, ratio of largest to smallest eigenvalues in a joint Eigen
        % decomposition.
        e = eig( squeeze(timeCovar(bini,:,:)) , squeeze(timeCovar(bini-1,:,:)) );
        y(bini,2) = e(end)/e(1);
    end
end

% Plot the distance results. They are scaled by their maximum values to be
% plotted on the same graph. They produce nearly identical results.
subplot(313)
plot(bsxfun(@rdivide,y,max(y)))
title('Covariance stationarity')

% Importantly, only the covariance distance metrics can successfully
% identify the non-stationarities in the covariances.

% To evaluate the statistical significance of stationarity, use the function mutualinformationX.
% For example, to test for mean stationarity, re-run the code and then:
[miX,entropyX,fd_binsX,pvalX] = mutualinformationx(timeMeanX,1:nbins,[],1);
[miY,entropyY,fd_binsY,pvalY] = mutualinformationx(timeMeanY,1:nbins,[],1);
% When the 4th input is a 1, the resulting pval (probability value) will be the probability that a 
% mutual information value as large as the observed value could have been obtained by chance. By
% statistical convention, if that value is smaller than 0.05 (5% probability of occuring by chance),
% it is deemed significant, which is to say, there is significant non-stationarity. Note that because
% of random distribution calculations, the p-value will change slightly each time it is called.


% Next, try re-running the code and then the following lines.
[miX,entropyX,fd_binsX,pvalX] = mutualinformationx(timevarX,1:nbins,[],1);
[miY,entropyY,fd_binsY,pvalY] = mutualinformationx(timevarY,1:nbins,[],1);
% Both mutual information values are statistically significant, indicate non-stationarity. 
% Mutual information is thus able to identify both linear and nonlinear relationships.

%% Chpater 9.6, Figure 9.7

% signal properties
srate = 1000;
t = 0:1/srate:9;
n = length(t);

% create double-chirp
freqTS   = linspace(12,30,n);
centfreq = mean(freqTS);
k        = (centfreq/srate)*2*pi/centfreq;
basesig  = sin(2*pi.*centfreq.*t + k*cumsum(freqTS-centfreq)) + sin(2*pi*10*t);

% plot time domain representation
figure(2),clf
subplot(211)
plot(t,basesig)


% setup parameters for wavelet convolution
wavetime    = -2:1/srate:2;
Lconv       = length(t)+length(wavetime)-1;
halfwavsize = floor(length(wavetime)/2);

nfrex = 50;
frex  = linspace(5,40,nfrex);
ncyc  = linspace(8,30,nfrex);

% loop through 5 simulations; 4 with different noise characteristics, and
% one with no noise for comparison.
for simi=1:5
    
    if simi==1
        % mean and variance stationary
        signal = basesig + randn(1,n)*2;
    elseif simi==2
        % mean non-stationary, variance stationary
        signal = basesig + linspace(1,5,n) + randn(1,n)*2;
    elseif simi==3
        % mean stationary, variance non-stationary
        signal = basesig + linspace(0,5,n) .* randn(1,n)*2;
    elseif simi==4
        % mean and variance non-stationary
        signal = basesig + linspace(0,5,n) + linspace(0,5,n) .* randn(1,n)*2;
    else
        % normal signal, for comparison
        signal = basesig;
    end
    
    % Fourier spectrum of signal
    sigX = fft(signal,Lconv);
    % initialize time-frequency matrix
    tf   = zeros(length(frex),length(t));
    
    % loop through frequencies
    for fi=1:nfrex
        % create complex Morlet wavelet
        w       = 2*( ncyc(fi)/(2*pi*frex(fi)) )^2; % w is width of Gaussian
        cmwX    = fft(exp(1i*2*pi*frex(fi).*wavetime) .* exp( (-wavetime.^2)/w ), Lconv);
        cmwX    = cmwX./max(cmwX); % amplitude scaling
        convres = ifft( sigX.*cmwX );
        convres = convres(halfwavsize:end-halfwavsize-1);
        
        tf(fi,:) = abs(convres).^2;
    end
    
    % Simulations 1-4 go in one plot; 
    % simulation 5 (without noise) goes in a separate plot.
    if simi<5
        figure(1), subplot(2,2,simi)
    else
        figure(2), subplot(223)
    end
    
    % and finally draw the plot
    contourf(t,frex,tf,40,'linecolor','none')
    set(gca,'clim',[0 .3])
end

%% end
