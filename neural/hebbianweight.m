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
