%{Generate a Poisson process with unit rate by generating a gamma order 2 spike train with 
parameter rho = 40. Then compute the transformed spike times and estimate numerically 
the resulting interspike interval distribution. Superpose on your plot the expected 
exponential distribution with unit rate.
%}

%generate nsamp samples of a gamma distribution
%of order 2 with 50 ms mean isi
a = 2; %equivalent of n
b = 25e-3; %equivalent of 1/varrho

nsamp = 100000;
vect = gamrnd(a,b,1,nsamp);
%mean(vect); %debug

[n,xout] = hist(vect,100);

%normalize to obtain a probability 
%density
dx = xout(2)-xout(1);
C = sum(n)*dx;
nr = n/C;
bar(xout,nr);

%compare with theoretical pdf
y = gampdf(xout,a,b);
axes(gca);
hold on;
plot(xout,y,'r');

rho = 1/b;
vect1 = cumsum(vect);
vtransf = zeros(size(vect));

%compute transformed times
vtransf(1) = rho*vect1(1) - log(1+rho*vect1(1));
for i = 2:length(vect)
    vtransf(i) = vtransf(i-1) + rho*(vect1(i) - vect1(i-1)) - log(1 + rho*(vect1(i)-vect1(i-1)));
end;

%compute transformed isis
dvtransf = diff(vtransf);

[nt,xt] = hist(dvtransf,100);
dxt = xt(2) - xt(1);
C1 = sum(nt)*dxt;
ntr = nt/C1;

figure;
bar(xt,ntr);

yt = exppdf(xt,1);
hold on;
plot(xt,yt,'r');
