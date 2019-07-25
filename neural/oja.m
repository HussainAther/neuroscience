% Linear associator with Hebb and weight decay: PCA a la Oja
clear; clf; hold on;

y=0; w=[0.1;0.4]; a=-pi/6;
rot=[cos(a) sin(a);-sin(a) cos(a)] % rotation matrix

for i:1000
    x=.05*randn(2,1); x(1)=4*x(1); % training examples
    x=rot*x; % rotation of training examples
    y=w'x; % network update
    
