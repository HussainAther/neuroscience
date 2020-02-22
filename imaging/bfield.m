%  Plot the magnetic field near a deoxygenated blood vessel
figure(1)
[x,y] = meshgrid(-1.5:.21:1.5,-1.5:.11:1.5);
Bx = 4/3*(x.^2+y.^2<=1) + ((x.^2-y.^2)./(x.^2+y.^2).^2/3 + 1).*(x.^2+y.^2>1);
By = ((2*x.*y)./(x.^2+y.^2).^2/3).*(x.^2+y.^2>1);
quiver(x,y,Bx,By,'k')
hold on
t = 0:pi/256:2*pi;
plot(cos(t),sin(t),'r')
hold off
axis image
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)

