%{
A series of experiments first reported in 1942 investigated the threshold of human subjects for 
detecting brief, weak light flashes. The experimental conditions were carefully optimized to maximize 
the sensitivity of human subjects. Prior to the task, the subjects were kept in the dark for at least 
30 mins to ensure full dark-adaptation of their visual system. The flashes were delivered at a horizontal 
distance 20 degrees away from the fovea in a region where the density of rod photoreceptors is high. 

The area covered by the stimulus (10 mins of arc) was also optimized to yield the highest sensitivity. 
Stimuli were presented for 1 ms and the wavelength of the light stimulus was 510 μm, a value at which the 
eye is known to be most sensitive for dim vision (Figure 19.2A). In the experiments, the energy of the 
light flash or equivalently the mean number of photons delivered at the cornea was varied and the 
frequency at which the observers detected the flashes was recorded.

This function plots the results.
"hsp" = "human single photon"
%}

function hsp

%first data set subject sh

%first column is mean no. of quanta
%second column is frequency of seeing in percent
%table V of HSP 1942 paper. Data set 1 for SH
sh_dat1 = ...
    [46.9 0.0;...
    73.1 9.4;...
    113.8 33.3;...
    177.4 73.5;...
    276.1 100.0;...
    421.7 100.0];

%log base 10 of the mean number of photons
x_sh1 = log10(sh_dat1(:,1));

%convert to percents
y_sh1 = sh_dat1(:,2)/100;

%plot the data
h_f1 = figure; 
h_a1 = subplot(2,2,1);
line('Parent',h_a1,'XData',x_sh1,'YData',y_sh1,'Marker','x','LineStyle','none');

%typical log10 range for the poisson cpf
x_cp1 = -0.25:0.01:1.25;

%compute cpf at 6 and no absorption factor, no dark photons
y_cp1 = p_c(x_cp1,6,1,0);

%index for the approximate mid-value of the cumulative function
[v_m1, ind_m1] = min(abs(y_cp1 - 0.5));
x_m1 = x_cp1(ind_m1);

%find the mid value of the data by linear interpolation. Since the last two
%values of y_hsp1 are 1.0, use only the first one for a vector of length 5.
x_mid1 = interp1(y_sh1(1:5),x_sh1(1:5),0.5);

%Use this value to translate the poisson cpf to get a good match
x_d1 = x_mid1 - x_m1;

line('Parent',h_a1,'XData',x_cp1+x_d1,'YData',y_cp1,'Color','r');
set(h_a1,'XLim',[1 3]);

%first data set subject sh
%load ss_dat1.txt
%first column is mean no. of quanta
%second column is frequency of seeing in percent
%table V of HSP 1942 paper. Data set 1 for SS
ss_dat1 = ...
    [24.1 0.0;...
    37.6 4.0;...
    58.6 18.0;...
    91.0 54.0;...
    141.9 94.0;...
    221.3 100.0];

%log base 10 of the mean number of photons
x_ss1 = log10(ss_dat1(:,1));

%convert to percents
y_ss1 = ss_dat1(:,2)/100;

%plot the data
h_a2 = subplot(2,2,2);
line('Parent',h_a2,'XData',x_ss1,'YData',y_ss1,'Marker','x','LineStyle','none');

%compute cpf at 6 and no absorption factor, no dark photons
y_cp2 = p_c(x_cp1,7,1,0);

%index for the approximate mid-value of the cumulative function
[v_m2,ind_m2] = min(abs(y_cp2 -0.5));
x_m2 = x_cp1(ind_m2);

%find the mid value of the data by linear interpolation. Since the last two
%values of y_hsp1 are 1.0, use only the first one for a vector of length 5.
x_mid2 = interp1(y_ss1,x_ss1,0.5);

%Use this value to translate the poisson cpf to get a good match
x_d2 = x_mid2 - x_m2;

line('Parent',h_a2,'XData',x_cp1+x_d2,'YData',y_cp2,'Color','r');
set(h_a2,'XLim',[1 3]);

%data set subject mhp
%first column is mean no. of quanta
%second column is frequency of seeing in percent
%table V of HSP 1942 paper. Data set 1 for MHP
mhp_dat1 = ...
    [37.6 6.0;...
    58.6 6.0;...
    91.0 24.0;...
    141.9 66.0;...
    221.3 88.0;...
    342.8 100.0];

%log base 10 of the mean number of photons
x_mhp1 = log10(mhp_dat1(:,1));

%convert to percents
y_mhp1 = mhp_dat1(:,2)/100;

%plot the data
h_a3 = subplot(2,2,3);
line('Parent',h_a3,'XData',x_mhp1,'YData',y_mhp1,'Marker','x','LineStyle','none');

%compute cpf at 5 and no absorption factor, no dark photons
y_cp3 = p_c(x_cp1,5,1,0);

%index for the approximate mid-value of the cumulative function
[v_m3, ind_m3] = min(abs(y_cp3-0.5));
x_m3 = x_cp1(ind_m3);

%find the mid value of the data by linear interpolation. Since the first two
%values of y_mhp1 are 0.6, use only the second one for a vector of length 5.
x_mid3 = interp1(y_mhp1(2:6),x_mhp1(2:6),0.5);
%x_mid3 = x_mhp1(3);

%Use this value to translate the poisson cpf to get a good match
x_d3 = x_mid3 - x_m3;

line('Parent',h_a3,'XData',x_cp1+x_d3,'YData',y_cp3,'Color','r');
set(h_a3,'XLim',[1 3]);

h_a4 = subplot(2,2,4);
x_cp = -1:0.01:2;
for i = 1:8
    y = p_c(x_cp,i,1,0);
    line('Parent',h_a4,'XData',x_cp,'YData',y);
end;

set(h_a4,'XLim',[-1 2]);
xlabel(h_a4,'log10 (mean nphotons/flash)')
ylabel(h_a3,'detection probability');

%print(handles.figure1,'-depsc2','hsp.eps');
end

function y = p_c(x,c,alpha,d_p)
%
% function y = p_c(x, c, alpha, d_p)
%  
% cumulative poisson function
% for photon counts
%  
% x is the log10 photon count at the retina
% alpha is the absorption yield factor
% c is the threshold for seeing 
% d_p is the dark photon number
  
%covert back x to linear space
xl = 10.^x;

%corresponding poisson parameter

p_p = alpha*xl + d_p;

y = 1 - poisscdf(c,p_p);
end
