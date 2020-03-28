%  Voltage response and phase plane for 1 Morris Lecar cell
%
%  usage   ml1pp(dt, Tfin, IsA)      IsA  = nA/cm^2
%
%  example    ml1pp(.001,8,550)
%

function ml1pp(dt, Tfin, IsA)

Cm = 1;
g = struct('K',20,'Ca',15,'Cl',5,'syn',0*30);
E = struct('K',-80,'Ca',100,'Cl',-50,'syn',-80);

Nt = ceil(Tfin/dt);

V = zeros(Nt,1); w = V;

V(1) = -50;
w(1) = minf(V(1));

for j=2:Nt,

    tmpt = tauw(V(j-1));
    tmpm = minf(V(j-1));
    tmps = sinf(fliplr(V(j-1)));
  
    w(j) = (tmpt.*w(j-1) + tmpm*dt)./(dt + tmpt);

    top = Cm*V(j-1) + dt*(g.Ca*tmpm*E.Ca + g.K*w(j)*E.K + g.Cl*E.Cl + ...
                          g.syn*tmps*E.syn + IsA);
    bot = Cm + dt*(g.Ca*tmpm + g.K*w(j) + g.Cl + g.syn*tmps);

    V(j) = top./bot;

end

t = linspace(0,Tfin,Nt);

figure(1)
plot(t',V,'k')
box off
set(gca,'tickdir','out')
xlabel('t  (ms)','fontsize',14)
ylabel('V  (mV)','fontsize',14)


figure(2)
plot(V,w,'r')
hold on
v = -80:.1:80;
plot(v,minf(v),'k--')
w = (IsA - g.Ca*minf(v).*(v-E.Ca) - g.Cl*(v-E.Cl))./(g.K*(v-E.K));
plot(v,w,'k:')
ylim([-.1 1.1])
set(gca,'tickdir','out')
xlabel('V  (mV)','fontsize',14)
ylabel('n','fontsize',14)
box off
hold off

return

function val = minf(V)
val = (1 + tanh(V/15))/2;
return

function val = tauw(V)
val = 1./cosh(V/30);
return

function val = sinf(V)
val = (1 + tanh(V/15))/2;
return

