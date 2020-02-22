%  Compute and Plot the eigenvalues of the Quasi-Active Cable
%

function quasicabspec

syms V m n h

gK = 36;             % mS / (cm)^2
gNa = 120;           % mS / (cm)^2
gl = 1/15; %0.3;            % mS / (cm)^2
VK = -77;               % mV
VNa = 56;              % mV
VL = -68;            % mV

tau = 1/gl;
rad = 1e-4;
R2 = 0.3; %0.034;
la2 = rad/2/R2/gl;

an = .01*(61+V)/(1-exp(-(V+61)/10));
bn = .125*exp(-(V+71)/80);
ninf = an/(an+bn);

am = .1*(46+V)/(1-exp(-(46+V)/10));
bm = 4*exp(-(V+71)/18);
minf = am/(am+bm);

ah = 0.07*exp(-(V+71)/20);
bh = 1/(exp(-(41+V)/10)+1);
hinf = ah/(ah+bh);

I = gK*(ninf^4)*(V-VK) + gNa*(minf^3)*hinf*(V-VNa) + gl*(V-VL);

i = inline(char(I));
Vr = fsolve(@(V) i(V), -71);

nr = subs(ninf,V,Vr);
mr = subs(minf,V,Vr);
hr = subs(hinf,V,Vr);

F(1) = -gK*(n^4)*(V-VK) - gNa*(m^3)*h*(V-VNa) - gl*(V-VL);

F(4) = an*(1-n) - bn*n;

F(2) = am*(1-m) - bm*m;

F(3) = ah*(1-h) - bh*h;

J = jacobian(F,[V m h n]);

B = subs(J,{V,m,h,n},{Vr,mr,hr,nr});

[W,Z] = eig(B)
c = W\[1 0 0 0]'
taum = -1/B(2,2); tauh = -1/B(3,3); taun = -1/B(4,4);
vNa = Vr - VNa; vK = Vr-VK;
dm = B(2,1)*taum; dh = B(3,1)*tauh; dn = B(4,1)*taun;

gamma = mr^3*hr*gNa + nr^4*gK + gl;

syms z

lab = {'(A)' '(B)'};

for k=1:2
    figure(k)
    ell = k*1e-1;
    for n=0:4
        P = (z+gamma+(n*pi/ell)^2*la2/tau).*(1+z*taum).*(1+z*tauh).*(1+z*taun) + ...
            3*mr^2*hr*vNa*dm*gNa.*(1+z*tauh).*(1+z*taun) + ...
            mr^3*vNa*dh*gNa.*(1+z*taum).*(1+z*taun) + ...
            4*nr^3*vK*dn*gK.*(1+z*taum).*(1+z*tauh);
        cP = sym2poly(P);
        zn = roots(cP)
        plot(real(zn),imag(zn),'ko','markersize',5*(n+1))
        hold on
    end
    hold off
    axis([-.35 -.05 -.25 .25])
    text(-.33,.2,lab{k},'fontsize',20)
    xlabel('real','fontsize',14)
    ylabel('imaginary','fontsize',14)
    %figname = ['quasicabspec' num2str(k)];
    %print('-depsc',figname) 
end

