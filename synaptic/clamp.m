% Simulate a voltage clamp experiment
%
% usage:  clamp(dt,Tfin)
%
% e.g. clamp(.01,15)
%

function clamp(dt,Tfin)

vK = -6;	% mV
GK =  36;	% mS/(cm)^2

vNa = 127;      % mV
GNa = 120;      % mS/(cm)^2

for vc = 8:10:90,

    j = 2;
    t(1) = 0;
    v(1) = 0;
    n(1) = an(0)/(an(0)+bn(0));  % 0.3177;
    m(1) = am(0)/(am(0)+bm(0));  % 0.0529;
    h(1) = ah(0)/(ah(0)+bh(0));  % 0.5961;

    gK(1) = GK*n(1)^4;
    gNa(1) = GNa*m(1)^3*h(1);

    while j*dt < Tfin,

        t(j) = j*dt;

        v(j) = vc*(t(j)>2)*(t(j)<Tfin);

        n(j) = ( n(j-1) + dt*an(v(j)) )/(1 + dt*(an(v(j))+bn(v(j))) ); 
        m(j) = ( m(j-1) + dt*am(v(j)) )/(1 + dt*(am(v(j))+bm(v(j))) );
        h(j) = ( h(j-1) + dt*ah(v(j)) )/(1 + dt*(ah(v(j))+bh(v(j))) );

        gK(j) = GK*n(j)^4;
        gNa(j) = GNa*m(j)^3*h(j);

        j = j + 1;

    end

    subplot(3,1,1); plot(t,v); hold on
    subplot(3,1,2); plot(t,gK); hold on
    subplot(3,1,3); plot(t,gNa); hold on
    %pause

end

subplot(3,1,1)
ylabel('v','fontsize',16)
hold off

subplot(3,1,2)
ylabel('g_K','fontsize',16)
hold off

subplot(3,1,3)
xlabel('t  (ms)','fontsize',16)
ylabel('g_{Na}','fontsize',16)
hold off

% rate functions from page 519 of HH
% (with polarity switched to agree with modern usage)

function val = an(v)
val = .01*(10-v)./(exp(1-v/10)-1);

function val = bn(v)
val = .125*exp(-v/80);

function val = am(v)
val = .1*(25-v)./(exp(2.5-v/10)-1);

function val = bm(v)
val = 4*exp(-v/18);

function val = ah(v)
val = 0.07*exp(-v/20);

function val = bh(v)
val = 1./(exp(3-v/10)+1);

