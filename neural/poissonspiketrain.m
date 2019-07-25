% Poisson psike train with refractoriness
clear; clf; hold on;
fr_mean=15/1000;

% Generate presynpatic Poisson spike trains.
lambda=1/fr_mean;
ns=1000;
isi1=-lambda.*log(rand(ns,1));

% Delete spikes in refractory period
is=0; for i=1:ns;
    if rand>exp(-isi1(i)^2/32);
        is=is+1;
        isi(is)=isi1(i);
    end
end

hist(isi, 50);
cv=std(isi)/mean(isi);
