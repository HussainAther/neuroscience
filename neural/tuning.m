% Tuning curve of a neuron

%peak firing rate 
pfr = 30; %in spk/s
vfr = 5; %variability in spk/s

%plot a tuning curve
dx = 10; %in deg
x = -180+dx:dx:180; %in deg
xr = x * (pi/180);

fr = pfr * cos(xr);

mfr = zeros(size(fr));
sfr = zeros(size(fr));

for i = 1:length(x)
    %draw a 100 random numbers
    fr_vals = randn(1,100)*vfr + fr(i);
    fr_vals = max(fr_vals,0);
    mfr(i) = mean(fr_vals);
    sfr(i) = std(fr_vals);
end;

h_f1 = figure;
h_a1 = subplot(2,2,1);
line('Parent',h_a1,'XData',x,'YData',mfr, ...
    'Marker','o','MarkerFaceColor','k','MarkerSize',4);
for i = 1:length(x);
    line('Parent',h_a1,'XData',[x(i) x(i)],'YData',[mfr(i)-sfr(i) mfr(i)+sfr(i)]);
end;

set(h_a1,'XLim',[-110 110],'YLim',[0 40]);
xlabel(h_a1,'stimulus angle (deg)');
ylabel(h_a1,'firing rate (spk/s)');


%plot stimulus and reconstruction
n_rand = 10; %number of random samples

%orthonormal basis corresponding to phi1 and phi2
e1 = [1;0];
e2 = [0;1];

%use a stimulus at 45 deg
x = 45;
xr = x*(pi/180);

%stimulus vector
fr1 = pfr*cos(xr);
fr2 = pfr*sin(xr);
v = [fr1;fr2];

%generate n_rand random responses
rfr1 = randn(1,n_rand)*vfr + fr1;
rfr2 = randn(1,n_rand)*vfr + fr2;

%compute the corresponding vector estimates
fr_mat = e1*rfr1 + e2*rfr2;

h_a2 = subplot(2,2,2);
line('Parent',h_a2,'XData',[0 v(1)],'YData',[0 v(2)],'Color','r');
for i = 1:n_rand
    line('Parent',h_a2,'XData',[0 rfr1(i)],'YData',[0 rfr2(i)]);
end;

set(h_a2,'XLim',[0 30],'YLim',[0 30]);
xlabel(h_a2,'firing rate (spk/s)');
ylabel(h_a2,'firing rate (spk/s)');

%compute the average mean square error between stimulus and reconstruction
n_rand = 100; %number of random samples

%use a range between 0 and 90 deg, since by symmetry the rest will be
%identical
x = -180+dx:dx:180; %in deg
xr = x * (pi/180);

rms_ang = zeros(size(x));
rms_errn = zeros(size(x));
se_errn = zeros(size(x));

for i=1:length(x)
    %stimulus vector
    fr1 = pfr*cos(xr(i));
    fr2 = pfr*sin(xr(i));
    v = [fr1;fr2];

    %generate n_rand random responses
    rfr1 = randn(1,n_rand)*vfr + fr1;
    rfr2 = randn(1,n_rand)*vfr + fr2;

    %compute the corresponding vector estimates
    fr_mat = e1*rfr1 + e2*rfr2;
    err_coord = fr_mat - repmat(v,[1 n_rand]);
    err_tot = sqrt(sum(err_coord.^2,1)); %distance between random sample and stimulus

    %root mean square error
    rms_err = mean(err_tot);
    std_err = std(err_tot);
    
    %normalize by vector length
    rms_errn(i) = rms_err/pfr;
    se_errn(i) = std_err/(pfr*sqrt(n_rand));
    
    %corresponding angle
    rms_ang(i) = atan(rms_err/pfr) * (180/pi);
end;

h_a3 = subplot(2,2,3);
line('Parent',h_a3,'XData',x,'YData',rms_errn, ...
    'Marker','o','MarkerFaceColor','k','MarkerSize',4);
for i = 1:length(x)
    line('Parent',h_a3,'XData',[x(i) x(i)],'YData',[rms_errn(i)-se_errn(i) rms_errn(i)+se_errn(i)]);
end;
set(h_a3,'YLim',[0.15 0.25]);
xlabel(h_a3,'stimulus angle (deg)');
ylabel(h_a3,'relative RMS');

h_a4 = subplot(2,2,4);
line('Parent',h_a4,'XData',x,'YData',rms_ang, ...
    'Marker','o','MarkerFaceColor','k','MarkerSize',4);
xlabel(h_a4,'stimulus angle (deg)');
ylabel(h_a4,'angular error (deg)');

%print(handles.figure1,'-depsc2','tuning_rec.eps');
