% Solve the static and (2) dynamic problems of K. Zhang, 1996, 
%     J Neurosci, 16(6), 2112-2126
% for a head-direction (head direction) cell ensemble
%  hdnet(dt,Tfin,pinc,pr)
%  "Representation of spatial orientation by the intrinsic dynamics of the head-direction cell ensemble: a theory"
%  e.g.,   hdnet(1,600,25)

function hdnet(dt,Tfin,pinc)

dtheta = pi/360;
th = -pi:dtheta:pi-dtheta;
th = th';

K = 8;
A = 1;
B = 39*exp(-K);
f = A + B*exp(K*cos(th));
N = length(f);

p = struct('beta', 0.8,'b',10,'c',0.5,'a',6.34);  % sigmoid parameters

u = isig(p,f);  % siginv(f)

fhat = fft(f)/N;
uhat = fft(u)/N;
fh2max = max(abs(fhat).^2);
lam = fh2max*1e-3;

What = uhat.*fhat./(lam + abs(fhat).^2);

W = real(N*ifft(What));

figure(1)
plot(th,fftshift(W),'k')
box off
set(gca,'tickdir','out')
xlabel('\theta','fontsize',14)
ylabel('W','fontsize',14)

figure(2)
plot(th,f,'k')
hold on
Wcf = N*ifft(What.*fhat);
frec = sig(p,Wcf);    % sig(w star f)
plot(th,frec,'r')
box off
hold off
set(gca,'tickdir','out')
legend('desired','actual')
xlim([-pi pi])
xlabel('\theta','fontsize',14)
ylabel('f  (Hz)','fontsize',14)

% and now some dynamics

tau = 10;
f0 = f/3 + circshift(f,180)/2 + 2*rand(N,1);
f  = f0;
u = isig(p,f);
figure(3)
plot3(th,zeros(N,1),real(f),'k');
hold on

Nt = ceil(Tfin/dt);

for j=2:Nt,
    
    t = (j-1)*dt;
    
    u = ((tau/dt)*u + ifft(What.*fft(f)) ) / (1 + tau/dt);
    
    f = sig(p,u);
    
    if mod(j,pinc)==0
       plot3(th,t*ones(N,1),real(f),'k');
    end
    
end

hold off
xlabel('\theta','fontsize',14)
ylabel('t  (ms)','fontsize',14)
zlabel('f  (Hz)','fontsize',14)
xlim([-pi pi])

% now demonstrate dynamic shift

dW = diff([W; W(1)])/dtheta;
gamma = 0.063;
w = W + gamma*dW;
figure(1)
hold on
plot(th,fftshift(w),'r')
hold off
legend('stationary','dynamic')
axis tight

f = f0;
u = isig(p,f);
figure(4)
plot3(th,zeros(N,1),real(f),'k');
hold on
what = fft(w)/N;

for j=2:Nt,
    
    t = (j-1)*dt;
    
    u = ((tau/dt)*u + ifft(what.*fft(f)) ) / (1 + tau/dt);
    
    f = sig(p,u);
    
    if mod(j,pinc)==0
       plot3(th,t*ones(N,1),real(f),'k');
    end
    
end

hold off
xlabel('\theta','fontsize',14)
ylabel('t  (ms)','fontsize',14)
zlabel('f  (Hz)','fontsize',14)
xlim([-pi pi])

figure(5)
u = -1:.01:1.5;
f = sig(p,u);
plot(u,f,'k')
xlabel('average synaptic input, u','fontsize',14)
ylabel('firing rate, f  (Hz)','fontsize',14)
box off 
set(gca,'tickdir','out')

function val = sig(p,u)
val = p.a*(log(1+exp(p.b*(u+p.c)))).^p.beta;

function val = isig(p,f)
val = (1/p.b)*log(exp((f/p.a).^(1/p.beta))-1) - p.c;
