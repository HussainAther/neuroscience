%{Power spectrum of a gamma renewal process of order 2, with a mean rate of 80 spk/s. 
Power spectrum of a gamma process of order 10.
%}

%plot the power spectrum of a gamma order
%two and 10 process.

close all; clear all;

%paramters of gamma process
rho = 160; %in spk/s
n_g = 2;

%mean rate in spk/s
chi = rho/n_g; 
f_v = [0:1:1000]; %in Hz = 1/s
cf_v = 2*pi*f_v;
cf_v2 = cf_v.^2;

%compute chi_c
k_v = 1;
n_k = length(k_v);
zk = exp(2*pi*i*k_v/n_g);
sk = rho*(zk-1);
skzk = sk.*zk;
sum_v = zeros(size(f_v));

for j = 1:n_k
    sum_v = sum_v + skzk(j)./(cf_v2 + sk(j)^2);
end;
ps_v = chi*(1-2*chi*sum_v);

%get rid of potential imaginary part due 
%to numerical errors
ps_v = real(ps_v);
figure(1);
plot(f_v,ps_v,'k');

set(gca,'TickDir','out');
ylabel('power spectral density ((spk/s)2/Hz)');
xlabel('frequency (Hz)');

%uncomment to print figure
%print(1,'-depsc2','gamma_ps_p1.eps');


%clear all;
rho = 800; %in spk/s
n_g = 10;

%mean rate in spk/sec
chi = rho/n_g;

k_v = [1:1:9];
n_k = length(k_v);
zk = exp(2*pi*i*k_v/n_g);
sk = rho*(zk-1);
skzk = sk.*zk;
sum_v = zeros(size(f_v));

for j = 1:n_k
    sum_v = sum_v + skzk(j)./(cf_v2 + sk(j)^2);
end;
ps_v = chi*(1-2*chi*sum_v);

%get rid of potential imaginary part due 
%to numerical errors
ps_v = real(ps_v);
figure(2);
plot(f_v,ps_v,'k');
set(gca,'TickDir','out');
xlabel('frequency (Hz)');
ylabel('power spectral density ((spk/s)2/Hz)');
set(gca,'TickDir','out');

%uncomment to print figure
%print(2,'-depsc2','gamma_ps_p2.eps');
