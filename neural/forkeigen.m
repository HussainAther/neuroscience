% Depict fork eigenfunctions and strength for full
% and reduced representations
%

close all

a = 1e-4*[1 1 1];
ell = 1e-4*[250 250 250];
dx = 1e-4;
N = ell/dx;
A3 = 2*pi*a(3)*dx;
As = 4*pi*1e-6;
rho = A3/As;
R2 = 0.3; %0.034;
gL = 1/15; %0.3;
lam = a/(2*R2*gL)/dx^2;   % lambda^2
r = a/a(3);
Hd = [2*lam(1)*ones(1,N(1)) 2*lam(2)*ones(1,N(2)) 2*lam(3)*ones(1,N(3)+1)];
Hd(1) = lam(1);
Hd(N(1)+1) = lam(2);
Hd(N(1)+N(2)+1) =  lam*r';
Hd(end) = rho*lam(3);

Hu = [-lam(1)*ones(1,N(1)-1) 0 -lam(2)*ones(1,N(2)) -lam(3)*ones(1,N(3))];
Hl = [-lam(1)*ones(1,N(1)-1) 0 -lam(2)*ones(1,N(2)-1) -r(2)*lam(2) -lam(3)*ones(1,N(3))];
Hl(end) = rho*Hl(end);

H = diag(Hd) + diag(Hu,1) + diag(Hl,-1);

H(N(1)+N(2)+1,N(1)) = -r(1)*lam(1);
H(N(1),N(1)+N(2)+1) = -lam(1);

Dd = [a(1)*ones(1,N(1)) a(2)*ones(1,N(2)) a(3)*ones(1,N(3)) a(3)/rho];
D = diag(Dd);
sD = diag((Dd).^(1/2));
sDi = diag((Dd).^(-1/2));
S = sD*H*sDi;

[Q,Z] = eig(S);

z = diag(Z);

[sz,si] = sort(z);

sz(1:10)

Qs = Q(:,si);

x3 = 1e4*(0:dx:ell(3));
x1 = 1e4*(ell(3):dx:ell(3)+ell(1));
x2 = 1e4*(ell(3):dx:ell(3)+ell(2));

for k=1:9,
    q = sDi*Qs(:,k+1);
    q1 = fliplr(q(1:N(1))');
    q2 = fliplr(q(N(1)+1:N(1)+N(2))');
    q3 = fliplr(q(N(1)+N(2)+1:end)');
    subplot(3,3,k)
    plot(x3,q3,'r')
    hold on 
    plot(x2,[q3(end) q2],'k')
    plot(x1,[q3(end) q1],'r--')
    hold off
    
end
    
W = sDi*Qs;

NN = size(W,1);

figure(2)

Vstr = zeros(NN);
Vstr01 = zeros(NN);

for m1=1:NN,

for m2=1:NN,
 
   for n=1:NN

       Vstr(m1,m2) = Vstr(m1,m2) + W(NN,n)*(W(m1,n)+W(m2,n))/(1+z(n));
  
   end

   for n=1:2
       Vstr01(m1,m2) = Vstr01(m1,m2) + W(NN,n)*(W(m1,n)+W(m2,n))/(1+z(n));
   end;
   
end

end     

Vstr = (1e-4)*Vstr/(2*pi*dx);

%for printing purposes, we subsample every 10 points and use the 
%painters render
ss_v = 1:10:NN;
mesh(ss_v,ss_v,Vstr(ss_v,ss_v));

colormap(1-gray)
xlabel('compartment 1','fontsize',14)
ylabel('compartment 2','fontsize',14)
zlabel('Strength  (mVms)','fontsize',14)
view(-45,45)
set(gca,'ZLim',[2 5]);

figure(3)

Vstr01 = (1e-4)*Vstr01/(2*pi*dx);

%full scale plot
%mesh(Vstr01);

%for printing, see above
mesh(ss_v,ss_v,Vstr01(ss_v,ss_v));

colormap(1-gray)
xlabel('compartment 1','fontsize',14)
ylabel('compartment 2','fontsize',14)
zlabel('Strength  (mVms)','fontsize',14)
view(-45,45)
set(gca,'ZLim',[2 5]);
