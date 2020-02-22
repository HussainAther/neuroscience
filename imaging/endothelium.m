%  Endothelial cell from 
%  K Dormanns, EMJ van Disseldorp, RG Brown and T David
%  Neurovascular coupling and the influence of luminal agonists
%  via the endothelium, J Theoretical Biology 364 (2015) 49-70
%
%  usage:    endothelium(dt,tfin,Jip3lum)
%
%  where:    dt = time step in seconds, e.g. dt = 0.005
%            tfin = final time in seconds, e.g., tfin = 100
%            Jip3lum = flux of IP3 triggered by lumen
%                    e.g. Jip3lum = 0.18 or 0.4
%
function endothelium(dt,tfin,Jip3lum)

close all
jfin = ceil(tfin/dt);
t = zeros(jfin,1);      % time

cae = zeros(jfin,1);     % [Ca^++]_endothelium
caer = zeros(jfin,1);     % [Ca^++]_endooplasmic reticulum
ip3e = zeros(jfin,1);     % [IP_3]_endothelium
Ve = zeros(jfin,1);     % endothelium, membrane potential

% establish resting values

cae(1) = 0.22; % 0.22; % try 0.078
CAM = cae(1);
caer(1) = 1.25; % 1.25;  % try 0.121
ip3e(1) = 0.05; %0.05;  % try 0.0175 
IP3M = ip3e(1);
Ve(1) = -75;    
VM = Ve(1);

rad = 20;
thick = 2;

for j=1:jfin-1,
    
    t(j+1) = dt*j;
    
    % compute endothelial fluxes
    
    i2 = ip3e(j)^2;
    JIP3e = 0.23*i2/(1+i2);          % (A.55)   micro molar / s
    
    c2 = cae(j)^2;
    JERup = 0.5*c2/(1+c2);        % (A.56)  
    
    ce2 = caer(j)^2;
    JCICRer = 5*(ce2/(4+ce2))*(c2^2/(0.9^4+c2^2));   % (A.57) 
    
    Jcaepump = 0.24*cae(j);   % (A.58)
    
    Jstretche = (6.1e-3)*(Ve(j)+18)/(1+exp(-7.4e-3*(30*rad/thick-500)));  %(A.59)
    
    JERleak = 0.025*caer(j);    % (A.60)
    
    logcae = log10(cae(j));
    Jcate = 6.6e-4*(50-Ve(j))*(1+tanh((logcae+0.18)/0.37))/2;  % (A.61)
    
    top = (logcae+0.4)*(Ve(j)+80.8) - 53.3;
    bot = 1.32e-3*(Ve(j)+53.3*(logcae+0.4)+80.8)^2 + 0.3;
    JBKe = 0.2*(1+tanh(top/bot));      % (A.63)
    
    JSKe = 0.3*(1+tanh((logcae+0.28)/0.389));    % (A.64)
    
    JKe = 6927*(Ve(j)+80)*(JBKe+JSKe);       % (A.62)
    
    JRe = 955*(Ve(j)+31.1);                 % (A.65)
    
    JIP3dege = 0.1*ip3e(j);                 % (A.66)
    
    J0e = 0.029;  % from nvu code
    
    % advance endothelial processes
    
    cae(j+1) = cae(j) + dt*( JIP3e - JERup + JCICRer - Jcaepump + ...
         JERleak + Jcate + J0e + Jstretche  - 0.05*(cae(j)-CAM)); 
      
    caer(j+1) = caer(j) + dt*( JERup - JCICRer - JERleak );
    
    ip3e(j+1) = ip3e(j) + dt*( Jip3lum - JIP3dege - 0.05*(ip3e(j)-IP3M) );
    
    Ve(j+1) = Ve(j) + dt*( -(JKe + JRe)/25.8 - 0.5*(Ve(j)-VM)); 
    
end

figure
plot(t,cae,'k','linewidth',1.5)
hold on
plot(t,caer,'k--','linewidth',1.5)
plot(t,ip3e,'r','linewidth',1.5)
%plot(t,Ve,'r--','linewidth',1.5)
hold off
legend('[Ca]_e','[Ca]_{er}','[IP_3]_e','location','best')
xlabel('t   (s)','fontsize',14)
ylabel('(\muM)','fontsize',14)
text(5,2.6,'(B)','fontsize',20)
box off

%figure
%plot(t,Ve,'k');
%xlabel('t   (s)','fontsize',14)
%ylabel('V_e   (mV)','fontsize',14)

return

