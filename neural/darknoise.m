%{The dark light hypothesis. The experiment discussed in "hsp.m" suggests that most of the variability in the observers’ responses 
is due to noise in the physical stimulus rather than biological noise. 

As pointed out a decade later, the experimental design and its interpretation have, however, several shortcomings:
1. If rods are indeed sensitive to single photons, why would observers not be as well, given that biological noise is assumed to be nonexistent?
2. The experiment described above is by itself somewhat ambiguous: an observer could always lower its threshold and thus give the appearance of “seeing” better.

A solution to these two problems is obtained by interpreting the results differently and by proposing a modified model of photon absorption. 
Although rods may be sensitive to single photons, it could be that several rods must be activated simultaneously when a weak flash is detected 
to overcome biological noise. One plausible source of noise is the random spontaneous decay of the rod photopigments (rhodopsin) in the absence 
of light. This decay would give the illusion of photon arrival and thus the registration of a single photon would in turn be unreliable to signal 
the presence of weak light flashes. Other sources of noise might result from central nervous system processing and can be lumped together with 
spontaneous rhodopsin decay for modeling purposes.

Let us assume that in the absence of light the mean number of absorbed photons (dark light) is x and follows a Poisson distribution. 
When presented with “blank” trials where no flash occurs, an observer is expected to report a light flash (even if none occurred) in a 
fraction of the trials because of this noise. If we call PFA the probability of such “false-alarms,” it is given by
P_FA = summ_k>=k_0 of ((x^k/k!)e^-x)
for amount of noise (x) and detection threshold (k_0).

In the presence of a light flash,  the mean number of absorbed photons will be due both to absorption related to the light flash, αn, 
and to the noise, x. If both processes follow independent Poisson distributions, their sum is also Poisson with mean a = αn + x. Thus:

P_D(a) = summ_k>=k_0 of (((αn + x)^k/k!)e^(-αn+x))

This function shows the example fit of the data taken from the same experiments with the dark noise model 
(the parameters are as follows: α = 0.13, x = 8.9, k0 = 21) and a  Fit of the dark noise model for one subject 
asked to be very conservative in detecting the flashes (crosses and red curve, probability of false-alarms equal 
to zero) and less conservative (pluses and black curve, probability of false-alarms equal to 0.1). 
The parameters are as follows: α = 0.13, x = 9.8, k0 = 19, 17.
%}

%this m-file reproduces the Barlow Fig. 2

function darknoise

%second data set
ss_dat2 = [ ...
    23.5 0.0;...
    37.1 0.0;...
    58.5 12.0;...
    92.9 44.0;...
    148.6 94.0;...
    239.3 100.0];

%log base 10 of the mean number of photons
x_hsp2s = log10(ss_dat2(:,1));

%convert to percents
y_hsp2s = ss_dat2(:,2)/100;

h_f1 = figure; 
h_a1 = subplot(1,2,1);
line('Parent',h_a1,'XData',x_hsp2s,'YData',y_hsp2s,...
    'Marker','x', 'LineStyle','none');

%compare with the Barlow model
x_hb = 1.3:0.01:2.7;
y_hb = p_c(x_hb,21,0.13,8.9);
 
line('Parent',h_a1,'XData',x_hb,'YData',y_hb,'Color','r');


%Barlow model and fig. 1 
x_hb1 = 1.3:0.01:2.3;
y_hb1 = p_c(x_hb1,19,0.13,9.8);
 
%note: we use 9.8 instead of 8.9 for the dark noise since it obviously 
%works much better with the data. 

h_a2 = subplot(1,2,2);
line('Parent',h_a2,'XData',x_hb1,'YData',y_hb1,'Color','r');

%Barlow model and fig. 1 
x_hb2 = 1.3:0.01:2.3;
y_hb2 = p_c(x_hb2,17,0.13,9.8);
 
line('Parent',h_a2,'XData',x_hb2,'YData',y_hb2,'Color','k');

%coordinates measured in points, relative to (1.5 0.2)
%0.2 in x direction is 28.88 pts
%0.2 in y direction is 28.58 pts
cf_x = 0.2/28.88;
cf_y = 0.2/28.58;

ptsx = [ 18.16; 38.1; 61.62; 82.16; 105.38 ];
ptsy = [ -7.74; 9.82; 54.48; 97.64; 107.17 ];

ptsx_c = 1.5 + ptsx*cf_x;
ptsy_c = 0.2 + ptsy*cf_y;

line('Parent',h_a2,'XData',ptsx_c,'YData',ptsy_c,...
    'LineStyle','none','Marker','x');
set(h_a2,'XLim',[1.3 2.3]);

ptsy2 = [ 12.21; 44.06; 83.65; 103.89; 111.04];
ptsy2_c = 0.2 + ptsy2*cf_y;
line('Parent',h_a2,'XData',ptsx_c,'YData',ptsy2_c,...
    'LineStyle','none','Marker','+');

xlabel(h_a1,'log10 nphotons');
ylabel(h_a1,'detection probability');

%print(handles.figure1,'-depsc2','darknoise.eps');

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
  
%convert back x to linear space
xl = 10.^x;

%corresponding poisson parameter

p_p = alpha*xl + d_p;

y = 1 - poisscdf(c,p_p);
