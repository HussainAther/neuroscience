% Illustrate the MRI process of slice selection

close all

figure(1)    % draw cube and slice
L = 1;
x = [L -L -L L L];
y = [L L -L -L L];
z = [1 1 1 1 1]*(-L);
plot3(x,y,z,'k')
hold on
plot3(x,y,-z,'k')
plot3([L L],[L L],[-L L],'k')
plot3(-[L L],-[L L],[-L L],'k')
plot3([L L],-[L L],[-L L],'k')
plot3(-[L L],[L L],[-L L],'k')
text(L,L,L+.1,'(L,L,L)')
text(-L,-L,-L-.1,'(-L,-L,-L)')
text(L,-L,-L-.1,'(L,-L,-L)')
text(0,-L,-L-.1,'x')
text(-L,0,-L-.1,'y')
text(L+.1,-L,0.125,'z_0')
zlo = 0*[1 1 1 1 1];
zhi = 0.25*[1 1 1 1 1];
plot3(x,y,zlo,'k')
plot3(x,y,zhi,'k')
patch(x,y,zhi,[1 1 1]/2)
patch(x,y,zlo,[1 1 1]/2)
patch([-L -L -L -L],[-L -L L L],[0 .25 .25 0],[1 1 1]/2)
patch([-L L L -L],[-L -L -L -L],[0 0 .25 .25],[1 1 1]/2)
text(L-.2,-L,L+.5,'(A)','fontsize',20)
axis equal
hold off
axis off
view(-38,16)

figure(2)    % draw sinc function
t = -1:.0001:1;
q = sin(6*pi*t)./(6*pi*t);
N = length(t);
q((N-1)/2+1) = 1;   % plug the NaN
plot(t,q,'k','linewidth',1.5)
xlabel('t  (ms)','fontsize',14)
ylabel('q','fontsize',14)
text(0.75,0.8,'(B)','fontsize',20)
box off

figure(3)     % draw Fourier Transform of sinc function
m = [-(N-1)/2:(N-1)/2]/2;
plot(m,abs(fftshift(fft(q))),'k','linewidth',1.5)
xlabel('\omega  (kHz)','fontsize',14)
ylabel('$|\hat q|$','interpreter','latex','fontsize',14)
xlim([-10 10])
text(7,1700,'(C)','fontsize',20)
box off

