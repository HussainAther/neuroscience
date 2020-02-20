% Compute eigenvalues and eigenststates of a matrix
% and find the area of the ellipse and parallelogram from them.
N = 100;
figure(1)
ang = linspace(0,2*pi,N);
x = [cos(ang); sin(ang)];	% a circle
plot(x(1,:),x(2,:),'k.')
hold on
B = [2 -1;-1 2];
y = B*x;			% a distorted circle
plot(y(1,:),y(2,:),'r.')
[V,D] = eig(B);
plot([0 D(1,1)*V(1,1)],[0 D(1,1)*V(2,1)],'k')
plot([0 D(2,2)*V(1,2)],[0 D(2,2)*V(2,2)],'k')
hold off
axis equal
set(gca,'tickdir','out')

figure(2)
t = linspace(0,1,25);
x = [t ones(1,25) fliplr(t) zeros(1,25)
     ones(1,25) fliplr(t) zeros(1,25) t];   % a square
plot(x(1,:),x(2,:),'k.')
hold on
B = [2 -1;-1 2];
y = B*x;				% a distorted square
plot(y(1,:),y(2,:),'r.') 
hold off
axis equal
set(gca,'tickdir','out')


