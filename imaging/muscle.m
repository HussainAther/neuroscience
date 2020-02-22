%  smooth muscle cell from 
%  K Dormanns, EMJ van Disseldorp, RG Brown and T David
%  Neurovascular coupling and the influence of luminal agonists
%  via the endothelium, J Theoretical Biology 364 (2015) 49-70
%
%  usage:    muscle(dt,tfin,stim)
%
%  where:    dt = time step in seconds, e.g. dt = 0.005
%            tfin = final time in seconds, e.g., tfin = 200
%            stim = stimulus of potassium in perivascular space and or
%                    calcium and/or IP3 through gap junctions with 
%                    endothelial cell,
%                    e.g. stim = struct('kp',10,'caec',2,'ip3ec',4)
%
function muscle(dt,tfin,stim)

K = [0 0.5 0.4 0.1 0.5 0 0.1];    % Hai-Murphy Crossbridge coeeficients
jfin = ceil(tfin/dt);
t = zeros(jfin,1);      % time

cam = zeros(jfin,1);     % [Ca^++]_muscle
casr = zeros(jfin,1);     % [Ca^++]_sarcoplasmic reticulum
ip3m = zeros(jfin,1);     % [IP_3]_muscle
Vm = zeros(jfin,1);     % muscle, membrane potential
rad = zeros(jfin,1);     % radius of blood vessel
Fr = zeros(jfin,1);     % fraction of bridges

% establish resting values

cam(1) = 0.22; CAEC = cam(1);
casr(1) = 1.25;
ip3m(1) = 0.05;  IP3EC = ip3m(1);
Vm(1) = -25;    VEC = Vm(1);
rad(1) = 19;
thick = 1.5;
kp = 3;        % perivascular K^+ concentration, mM

ym =  cam(1)^2/(cam(1)^2 + 0.13*exp(-(Vm(1)+27)/12));  % KCa gate

K(1) = 17*cam(1)^3;
K(6) = K(1);
HMA = [-(K(1)+K(2)+K(3)) K(4)-K(1) -K(1); 
         K(3) -(K(4)+K(5)) K(6);
         0  K(5) -(K(6)+K(7))];
HMb = [-K(1); 0; 0];

HMrest = HMA \ HMb;
Mp = HMrest(1);
AMp = HMrest(2);
AM = HMrest(3);
Fr(1) = AMp + AM;

for j=1:jfin-1,
    
    t(j+1) = dt*j;
    
    % compute muscle fluxes
    
    i2 = ip3m(j)^2;
    JIP3m = 0.23*i2/(1+i2);          % (A.42)   micro molar / s
    
    c2 = cam(j)^2;
    JSRup = 2.025*c2/(1+c2);        % (A.43)  micro molar / s
    
    cs2 = casr(j)^2;
    JCICRsr = 55*(cs2/(4+cs2))*(c2^2/(0.9^4+c2^2));   % (A.44) micro molar / s
    
    Jcampump = 0.24*cam(j)*( 1 + (Vm(j)+100)/250 );   % (A.45)
    
    JSRleak = 0.025*casr(j);    % (A.46)
    
    JCAm = (1.29e-3)*(Vm(j)-100)/(1+exp(-((Vm(j)+24)/8.5)));   % (A.47)
    
    JNaCam = (3.16e-3)*cam(j)*(Vm(j)+30)/(cam(j)+0.5);     % (A.48)
    
    Jstretchm = (6.1e-3)*(Vm(j)+18)/(1+exp(-7.4e-3*(30*rad(j)/thick-500)));  %(A.49)
    
    JNaKm = 4.32e-2;    % (A.50)
    
    JClm = (1.34e-3)*(Vm(j)+25);   % (A.51)
    
    JKCam = (4.46e-3)*ym*(Vm(j)+94);    % (A.52)
    
    Kactm =  c2/(c2 + 0.13*exp(-(Vm(j)+27)/12));    % (A.70) + table below it
    
    tkp = kp + stim.kp*(t(j)>60)*(t(j)<120);
    EKIR = 4.5*tkp - 112;                    % Table 1, p 53 and (2.7)
    gKIR = exp((0.42)*tkp - (7.4e-2)*Vm(j) - 12.6);    % Table 1 and (2.8)
    JKIRm = (750/1970)*gKIR*(Vm(j)-EKIR);     % (A.53)
    
    JIP3deg = 0.1*ip3m(j);
    
    % advance muscle processes
    
    ym = ym + dt*45*(Kactm-ym);      % (A.35) gating of muscle KCa channel
    
    CAECt = CAEC + stim.caec*(t(j)>60)*(t(j)<120);
    cam(j+1) = cam(j) + dt*( JIP3m - JSRup + JCICRsr - JCAm + ...
        JNaCam - Jcampump + 0.1*Jstretchm + JSRleak - 0.05*(cam(j)-CAECt)); 
      
    casr(j+1) = casr(j) + dt*( JSRup - JCICRsr - JSRleak );
    
    IP3ECt = IP3EC + stim.ip3ec*(t(j)>60)*(t(j)<120);
    ip3m(j+1) = ip3m(j) + dt*( -0.05*(ip3m(j)-IP3ECt) - JIP3deg );
    
    Vm(j+1) = Vm(j) + dt*( 1970*(-JNaKm - JClm - 2*JCAm - ...
           JNaCam - JKCam - JKIRm - Jstretchm)- 0.5*(Vm(j)-VEC));
    
    K(1) = 17*cam(j)^3;
    K(6) = K(1);
    Mp = Mp + dt*((K(4)-K(1))*AMp + K(1) - K(1)*AM - (K(1)+K(2)+K(3))*Mp);
    AMp = AMp + dt*(K(3)*Mp + K(6)*AM - (K(4)+K(5))*AMp);
    AM = AM + dt*(K(5)*AMp - (K(6)+K(7))*AM);
    Fr(j+1) = AMp + AM;
    
    thick = -rad(j) + sqrt(rad(j)^2 + 2*(20)*(3) + (3)^2);  % (A.78)
    EFr = 66 + Fr(j)*(233-66);
    R0 = 20 + Fr(j)*(0.6-1)*20;
    rad(j+1) = rad(j) + dt*(20/10)*(rad(j)*4/thick - EFr*(rad(j)-R0)/R0);
    
end

figure
plot(t,cam,'k','linewidth',1.5)
hold on
plot(t,casr,'k--','linewidth',1.5)
plot(t,ip3m,'r','linewidth',1.5)
plot(t,Fr,'r--','linewidth',1.5)
hold off
legend('[Ca]_m','[Ca]_{sr}','[IP_3]_m','AM^*','location','best')
xlabel('t   (s)','fontsize',14)
ylabel('(\muM)','fontsize',14)

figure
[ax,h1,h2] = plotyy(t,Vm,t,rad);
%ylim(ax(1),[-60 -10])
set(ax(1),'ycolor',[1 0 0])
set(h1,'color','red','linewidth',1.5)
set(ax(2),'ycolor',[0 0 0])
set(h2,'color','black','linewidth',1.5)
xlabel('t   (s)','fontsize',14)
ylabel(ax(1),'V_m   (mV)','fontsize',14,'color','red')
ylabel(ax(2),'\rho   (\mum)','fontsize',14)

return
%  Plot.
close all
stim = struct('kp',10,'caec',0,'ip3ec',0);
muscle(.005,200,stim)
figure(1)
text(20,1,'(A)','fontsize',20)
figure(2)
text(20,-10,'(B)','fontsize',20)

stim = struct('kp',0,'caec',2,'ip3ec',0);
muscle(.005,200,stim)
figure(3)
text(20,1,'(C)','fontsize',20)
figure(4)
text(20,-15,'(D)','fontsize',20)
