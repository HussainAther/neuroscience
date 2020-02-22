% Solve for the extracellular potential in the neighborhood
% of an active cable, and examine the dependence on gK and gNa

function stEcabextra
close all
g = struct('K', 36, 'Na', 120, 'Cl', 0.3);
[t,phibase] = stEcabextrag(g);
g.K = 1.1*g.K;
[t,phiHiK] = stEcabextrag(g);
g.K = 36;
g.Na = 1.1*g.Na;
[t,phiHiNa] = stEcabextrag(g);
figure
plot(t,squeeze(phiHiNa(1,5,:)-phibase(1,5,:)),'k')
hold on
plot(t,squeeze(phiHiK(1,5,:)-phibase(1,5,:)),'r')
legend('Hi g_{Na}','Hi g_K')
hold off
xlabel('t')
ylabel('difference in \phi')

function [t,phi] = stEcabextrag(g)

cab = struct('rad',1e-4,'ell',1e-1,'dx',1e-4,'dt',0.01);
stim = struct('t1',1,'t2',2,'amp',2e-4,'loc',0.05,'Tfin',8);

pinc = 20;

Cm = 1;		% micro F / cm^2
R2 = 0.3; %0.034;		% k Ohm cm
dx = cab.dx;
dt = cab.dt;
Nx = cab.ell/dx;		% patch length
A = 2*pi*cab.rad*dx;		% patch surface area
x = dx/2:dx:cab.ell-dx/2;	% vector of patch midpoints

Nt = ceil(stim.Tfin/dt);

E = struct('K', -77, 'Na', 56, 'Cl', -68);

eloc = round(Nx*stim.loc/cab.ell);
Iapp = stim.amp/A;

e = ones(Nx,1);
B = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
B(1,1) = 1/dx/dx;
B(Nx,Nx) = 1/dx/dx;
B = (cab.rad/2/R2)*B;
dB = diag(B);

options = optimset('Jacobian','on');
V = fsolve(@(V) Iss(V,E,g,B),-70*e,options);    % initial conditions

n = an(V)./(an(V)+bn(V)); 
m = am(V)./(am(V)+bm(V)); 
h = ah(V)./(ah(V)+bh(V)); 

t = 0;
if pinc > 0
   figure(1)
   plot3(x,t*ones(Nx,1),V,'k')
   hold on
end

I = zeros(Nx,1);
phi = zeros(4,9,Nt);
sig = 0.003; % S/cm % page 14 Magneto ref
ycenters = ((1:Nx)'-1/2)*cab.ell/Nx;

for j=2:Nt,

      t = (j-1)*dt;

      et = max(t-stim.t1,0);
      I(eloc) = Iapp.*et.*exp(1-et);
      %I(eloc) = Iapp.*(t-dt/2>stim.t1).*(t-dt/2<stim.t2);

      a = an(V);  c = (a+bn(V))/2;
      n = ( (1/dt-c).*n + a) ./ (1/dt + c); n4 = n.^4;

      a = am(V);  c = (a+bm(V))/2;
      m = ( (1/dt-c).*m + a) ./ (1/dt + c);

      a = ah(V);  c = (a+bh(V))/2;
      h = ( (1/dt-c).*h + a) ./ (1/dt + c); m3h = m.^3.*h;

      d = g.Na.*m3h + g.K.*n4 + g.Cl;

      f = g.Na.*m3h*E.Na + g.K.*n4*E.K + g.Cl.*E.Cl + I;

      B(1:Nx+1:end) = dB + d + 2*Cm/dt;         % update the diagonal

      Vmid = B\(2*Cm*V/dt + f);
      
      V = 2*Vmid - V;

      if mod(j,pinc) == 0
         plot3(x,t*ones(Nx,1),V,'k')
      end
  
      % collect the membrane currents

    INa = g.Na*m3h.*(V-E.Na);
    IK = g.K*n4.*(V-E.K);
    ICl = g.Cl*(V-E.Cl);
    ICap = 2*Cm*(V-Vmid)/dt;
    Imem = A*(INa + IK + ICl - I + ICap);
    for ix = 1:4
        xrec = 0.005*ix;
        for iy = 1:9
            yrec = 0.01*iy;
            dscale = 4*pi*sig*sqrt(xrec^2+(yrec-ycenters).^2);
            phi(ix,iy,j) = sum(Imem./dscale);
        end
    end

end

if pinc
   xlabel('x  (cm)','fontsize',16)
   ylabel('t  (ms)','fontsize',16)
   zlabel('V  (mV)','fontsize',16)
   hold off
end

t = linspace(0,stim.Tfin,Nt)';
tim = (1:Nt)*dt;
figure(2)

for ix=1:4, 
    for iy = 1:9, 
        plot(tim+(5/4)*(ix-1)*Nt*dt,squeeze(phi(ix,iy,:))+5*(iy-1),'k','linewidth',2)
        hold on
    end
end

X = [-2 -1 -1 -2];
Y = [-5 -5 45 45];
fill(X,Y,'k')   % draw cable
for j=1:9
    text(-8,(j-1)*5,[num2str(j*100) ' \mum']); % vertical labels
end
for j=1:4,
    text(1+(j-1)*10,-3,[num2str(50*j) ' \mum']) % horizontal labels
end

plot([8 11],[39 39],'k')
text(8,38,'3 ms')

plot([9 9],[41 44],'k')
text(9.5,43,'3 \muV')

hold off
axis off
box off
%axis tight

return

function [val, jac] = Iss(V,E,g,B)
Nx = length(V);
a = am(V); b = bm(V);
m = a./(a+b);
dm = (dam(V).*b - a.*dbm(V))./(a+b).^2;

a = ah(V); b = bh(V);
h = a./(a+b);
dh = (dah(V).*b - a.*dbh(V))./(a+b).^2;

a = an(V); b = bn(V);
n = a./(a+b);
dn = (dan(V).*b - a.*dbn(V))./(a+b).^2;

m3h = m.^3.*h;
n4 = n.^4;

val = B*V + g.Na.*m3h.*(V-E.Na) + g.K.*n4.*(V-E.K) + g.Cl.*(V-E.Cl);

dj = g.Na.*((3*dm.*m.^2.*h + m.^3.*dh).*(V-E.Na) + m3h) + ...
     g.K.*(4*dn.*n.^3.*(V-E.K) + n4) + g.Cl;
jac = B + spdiags(dj,0,Nx,Nx);

function val = an(v)
val = .01*(10-(v+71))./(exp(1-(v+71)/10)-1);

function val = dan(v)
tmp = exp(-(61+v)/10);
val = -( tmp.*(71+v) - 10 )./(tmp-1).^2/1000;

function val = bn(v)
val = .125*exp(-(v+71)/80);

function val = dbn(v)
val = -exp(-(v+71)/80)/640;

function val = am(v)
val = .1*(25-(v+71))./(exp(2.5-(v+71)/10)-1);

function val = dam(v)
tmp = exp(-(46+v)/10);
val = -( tmp.*(56+v) - 10 )./(tmp-1).^2/100;

function val = bm(v)
val = 4*exp(-(v+71)/18);

function val = dbm(v)
val = -(2/9)*exp(-(v+71)/18);

function val = ah(v)
val = 0.07*exp(-(v+71)/20);

function val = dah(v)
val = -(7/2000)*exp(-(v+71)/20);

function val = bh(v)
val = 1./(exp(3-(v+71)/10)+1);

function val = dbh(v)
tmp = exp(-(v+41)/10);
val = tmp./(tmp+1).^2/10;

