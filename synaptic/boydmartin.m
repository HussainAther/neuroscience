%The following numerical values reproduce the distribution
%of spontaneous e.p.p.s described by Boyd and Martin (1956).
%The first column is the depolarization in mV and the second 
%column is the number of observations. 
spont_data = [0.2 1;0.3 25;0.4 30;0.5 20;0.6 2];

x1 = spont_data(:,1);
dx = x1(2)-x1(1); %recover the sampling interval

y1 = spont_data(:,2);

y1p = y1./sum(y1); %normalize to get probabilities
mu1 = sum(x1.*y1p); %mean depolarization

%reconstruct a vector of observations appropriate for
%normfit
x1_obs = [];
for i = 1:length(x1)
    x1_obs = [ x1_obs x1(i)*ones(1,y1(i))]; 
end;

%estimate mean and standard deviation    
[muhat,sigmahat,muci,sigmaci] = normfit(x1_obs);

%sample the Gaussian at higher resolution
dx_f = 0.01;
x_f = 0.0:dx_f:0.8;
yn = normpdf(x_f,muhat,sigmahat);

%normalize so that the area under both curves is identical
%(area under yn with sampling step dx_f is 1 by construction, 
%i.e., sum(yn*dx_f) = 1)
yn_rs = yn*sum(y1*dx);

%plot everything
h_f1 = figure; 
h_a1 = axes;
hold on;
bar(x1,y1,'k');
plot(x_f,yn_rs,'r');
xlabel('spontaneous e.p.p.s. (mV)');
ylabel('number of observations');

%This is the stimulation data of Boyd and Martin (1956).
%See spont_data file for a detailed description.
stim_data = [0 18;0.3 11;0.4 20;0.5 13;0.6 6;0.7 14;0.8 18;0.9 17;1.0 6;...
    1.1 11;1.2 10;1.3 9;1.4 4;1.5 7;1.6 9;1.7 5;1.8 5;1.9 3;2.0 2;2.1 2;...
    2.2 1;2.3 1;2.4 2;2.5 1;2.6 1;2.7 1;2.8 0;2.9 0;3.0 1];

x2 = stim_data(:,1);
y2 = stim_data(:,2);
n_x2 = length(x2);

y2p = y2./sum(y2);
mu2 = sum(x2.*y2p); %mean depolarization

l_e = mu2/mu1; %quantal size, standard method
l_f = -log(y2(1)/sum(y2)); %method of failures

%relative error
l_err = abs(l_e - l_f)/(0.5*(l_e + l_f));

info_str = sprintf('relative error: %.3f',l_err);
disp(info_str);

%figure(2)
h_f2 = figure; 
h_a2 = axes;
hold on;
bar(x2,y2,'k');

%sample again at high resolution
x3 = 0:dx_f:3.2;
n_x3 = length(x3);

y3 = zeros(1,n_x3);
y3(1) = poisspdf(0,l_e); %discrete probability of failures

%cycle over the peaks and implement the formula in the lecture notes
for i=1:7
    pois_val = poisspdf(i,l_e);
    mu = i*muhat; %mean of the ith peak
    sigma = sqrt(i)*sigmahat; %standard deviation
    yi = normpdf(x3,mu,sigma);
    y3= y3+pois_val*yi; %add everything toghether
end;

%normalize so that the area under everything but the first element 
%is equal to 1 - the probability of the 1st element (failures)
%this is the correct probability density
y3(2:n_x3)= y3(2:n_x3)*(1-y3(1))/sum(y3(2:n_x3)*dx_f);

%correctly normalized probability density
y3n = zeros(1,n_x3);

%the fraction of failures is the total sum of events times the probability
%of failures
y3n(1) = sum(y2)*y3(1);

%normalize so that the area under the continuous part is equal to 
%the area under the sample histogram
y3n(2:n_x3) = y3(2:n_x3)*(sum(y2(2:n_x2)*dx))/(sum(y3(2:n_x3)*dx_f));

hold on;
plot(x3,y3n,'r');
xlabel('evoked e.p.p.s. (mV)');
ylabel('number of observations');
set(h_a2,'XLim',[-0.25 3.25]);

%plot of the mean obtained by by both methods
m_vals = [3.36 3.22; 2.89 3.19; 2.64 2.72; 2.33 2.40; 1.86 1.95; ...
    1.50 1.48; 1.15 1.08; 1.04 1.08; 0.58 0.59; 0.38 0.31];

h_f3 = figure; 
h_a3 = axes;
line('Parent',h_a3,'XData',m_vals(:,1),'YData',m_vals(:,2),...
    'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerSize',4);

line('Parent',h_a3,'XData',[0 4],'YData',[0 4],'LineStyle','--');
xlabel(h_a3,'\lambda (potentials)');
ylabel(h_a3,'\lambda (failures)');

%print(handles.figure1,'-depsc2','boyd_martin2.eps');

