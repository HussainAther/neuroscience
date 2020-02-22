%  Plot two approximations involving the the zeroth order
%  Bessel function

figure(1)
x = 0:1e-5:2;
y1 = besselj(0,x);
y2 = 1 - x.^2/4;
plot(x,y1,'k','linewidth',1.5)
hold on
plot(x,y2,'r','linewidth',1.5)
legend('J_0(x)','1-x^2/4','location','SW')
text(1.8,0.9,'(A)','fontsize',20)
xlabel('x','fontsize',14)
%ylabel('J_0','fontsize',14)
hold off
box off

figure(2)
u = 1e-5:1e-5:1;
for m=1:40,

    phi(m) = (9+m)/10;
    fc(m) = trapz(u,(1-besselj(0,phi(m)*u))./(u.^2));
    est(m) = phi(m) - 1;
   
end
plot(phi,fc,'k','linewidth',1.5)
hold on
plot(phi,est,'r','linewidth',1.5)
legend('F(\phi)','\phi-1','location','SE')
text(1.1,3.5,'(B)','fontsize',20)
xlabel('\phi','fontsize',14)
hold off
box off
