% Convex convolution
t = 0:.01:10;
t1 = 1;
t2 = 2;
f = exp(-t/t1);
g = exp(-t/t2);
fcg = (t1*t2)*(exp(-t/t1)-exp(-t/t2))/(t1-t2);
plot(t,f,'k')
hold on
plot(t,g,'k--')
plot(t,fcg,'r')
hold off
box off
legend('exp(-t)','exp(-t/2)','their convolution')
xlabel('t','fontsize',14)

