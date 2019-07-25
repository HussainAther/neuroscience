% Linear associator with Hebb and weight decay: PCA a la Oja
clear; clf; hold on;

y=0; w=[0.1;0.4]; a=-pi/6;
rot=[cos(a) sin(a);-sin(a) cos(a)] % rotation matrix

for i:1000
    x=.05*randn(2,1); x(1)=4*x(1); % training examples
    x=rot*x; % rotation of training examples
    y=w'x; % network update
    plot(x(1), x(2), '.');
    w=w+.1*y*(x-y*2); % training
    w_traj(:,1)=w; % recording of weight history
end

plot(w_traj(1,:), w_traj(2,:),'r');
plot([0 w(1)], [0 w(2)], 'k', 'linewidth', 2);
axis([-1 1 -1 1);
plot([-1 1], [0 0], 'k'); 
plot([0 0], [-1 1], 'k');
