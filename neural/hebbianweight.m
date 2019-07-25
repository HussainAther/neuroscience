% Hebbian synapses in rate model
clear; clf; hold on;
nn=500; %try alos nn=5000
npat=1000;

% random pattern; firing rates are exponentially distributed
afr=40;
r1=-afr.*log(rand(nn,npat));
r1=-afr.*log(rand(1,npat));

w=(r2-afr)*(r1-afr)';
w=w/npat

x=-260:20:260;
[n, x]=hist(w, x);
n=n/sum(n)/20;
h=bar(x, n); set(h, 'facecolor', 'none');

% Fit
a0=[0 40];
a=lsqcurvefit('normal', a0, x, n);
n2=normal(a, x);
plot(x, n2, 'r');
