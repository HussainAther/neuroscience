
%{
Model lateral geniculate nucleus (LGN) neuron that receives a target bar
of light preceded by two masking bars of light. If positioned in space
and time appropriately, the target bar should be less visible.
}%

clear all;
close all;

%Initialize parameters.
sigma_c = 1.4; %width of center portion of spatial r.f. [degrees]
sigma_s = 2.1; %determines width of surround portion of spatial r.f. [deg]
A_c = 1; %strength of center portion of spatial r.f.
A_s = 0.9; %strength of surround portion of spatial r.f.
dx = 0.05; %resolution of spatial grid
x_min = -4*sigma_s; %x-value roughly giving -x border of r.f.
x_max = 4*sigma_s; %x-value roughly giving +x border of r.f.
x_vect = x_min:dx:x_max; %values of x over which to compute D_x
alpha = 1/10; %determines length of temporal filter [ms^-1]
dt = 1; %ms
tau_vect_max = 20/alpha; %tau value giving border of temporal r.f.
tau_vect = 0:dt:tau_vect_max; %value of tau over which to compute D_t

%Define and then plot spatial kernel D_x.
D_x_center = A_c*exp(-(x_vect.^2)/(2*(sigma_c^2)))/sqrt(2*pi*sigma_c^2);
D_x_surround = A_s*exp(-(x_vect.^2)/(2*(sigma_s^2)))/sqrt(2*pi*sigma_s^2);
D_x_vect = D_x_center - D_x_surround;
figure(1)
subplot(3,1,1)
%plot(x_vect, D_x_center, "r--") %plot center with red dashed lines
hold on
%plot(x_vect, D_x_surround, "k--") %plot surround with black dashed lines
plot(x_vect, D_x_vect) %plot full receptive field with solid blue line
xlabel("x (deg)")
ylabel("D_x")

%Define and then plot temporal kernel D_t.
D_t_vect = alpha*exp(-alpha*tau_vect).*((alpha*tau_vect).^5/(5*4*3*2) - ...
 (alpha*tau_vect).^7/(7*6*5*4*3*2));
subplot(3,1,2)
plot(tau_vect,D_t_vect)
set(gca,"XDir","reverse")
xlabel("tau (ms)")
ylabel("D_t")

%Define and then plot the full spatio-temporal kernel D(x,tau).
D_xt_mat = D_x_vect"*D_t_vect; %full 2-D r.f., with x as 1st dimension & t as 2nd dimension
subplot(3,1,3)
contour(tau_vect,x_vect,D_xt_mat,12); %makes a contour plot of the data
colorbar

set(gca,"Xdir","reverse")
xlabel("tau (ms)")
ylabel("x (deg)")

%spatial parameters & plot of spatial locations of stimuli
Target_LeftEnd = -0.5; %position of start of bar [deg]
Target_RightEnd = 0.5; %position of end of bar [deg]
Target_x_vect = [zeros(1,(Target_LeftEnd - x_min)/dx)... %nothing flashed here
                 ones(1,((Target_RightEnd - Target_LeftEnd)/dx)+1)... %Target position
                 zeros(1,(x_max - Target_RightEnd)/dx)]; %nothing flashed here
Mask1_LeftEnd = -1.5;
Mask1_RightEnd = -0.5;
Mask2_LeftEnd = 0.5;
Mask2_RightEnd = 1.5;
Mask_x_vect = [zeros(1,(Mask1_LeftEnd - x_min)/dx) ... %nothing flashed here
               ones(1,((Mask1_RightEnd - Mask1_LeftEnd)/dx)+1) ... %Target position
               zeros(1,((Mask2_LeftEnd - Mask1_RightEnd)/dx)-1) ... %nothing flashed here
               ones(1,((Mask2_RightEnd - Mask2_LeftEnd)/dx)+1) ... %Target position
               zeros(1,(x_max - Mask2_RightEnd)/dx)]; %nothing flashed here

figure(2)
subplot(3,1,1)
plot(x_vect,Target_x_vect,"o")
xlabel("x (deg)")
ylabel("Target(1=Lt,0=Dk)")

%temporal parameters & plot of temporal locations of stimuli
tmax = 800; %ms
Target_on = 400; %ms
Target_off = 600; %ms
t_vect = 0:dt:tmax;
Target_t_vect = [zeros(1,Target_on/dt) ... %bar initially off from t=0 to t=Target_on-dt
                 ones(1,(Target_off-Target_on)/dt) ... %then Target on
                 zeros(1,((tmax-Target_off)/dt)+1)]; %then Target off again
Mask_on = 200; %ms
Mask_off = 400; %ms
Mask_t_vect = [zeros(1,Mask_on/dt) ... %bar initially off from t=0 to t=Mask_on-dt
               ones(1,(Mask_off-Mask_on)/dt) ... %then Mask on
               zeros(1,((tmax-Mask_off)/dt)+1)]; %then Mask off again
subplot(3,1,2)
plot(t_vect,Target_t_vect)
xlabel("t (ms)")
ylabel("Target")
%define and plot stimulus in both space & time
Target_xt_mat = Target_x_vect"*Target_t_vect;
Mask_xt_mat = Mask_x_vect"*Mask_t_vect;
Both_xt_mat = Target_xt_mat + Mask_xt_mat; %mask and target both presented
subplot(3,1,3)
contourf(t_vect,x_vect,Target_xt_mat); %makes a contour plot of the data
colorbar
xlabel("t (ms)")
ylabel("x (deg)")

%Run model.
r0 = 10e-3; %background rate (ms^-1)
Target_L_x = dx*D_x_vect*Target_x_vect"; %spatial linear filter integral
Mask_L_x = dx*D_x_vect*Mask_x_vect"; %spatial linear filter integral
%prepare to do temporal filter integral
Target_L_t_vect = zeros(1,length(t_vect)); %set up vector to hold temporal filter values
Mask_L_t_vect = zeros(1,length(t_vect)); %set up vector to hold temporal filter values
Stim_t_NegTimes = zeros(1,tau_vect_max+1); %need to add stimulus values at t<0 so
 %can get what influenced cell at times t<tau_vect_max
Target_t_vect_long = [Stim_t_NegTimes Target_t_vect]; %includes Stimulus at times t<0
Mask_t_vect_long = [Stim_t_NegTimes Mask_t_vect]; %includes Stimulus at times t<0
i = 0;
for t=0:dt:tmax
   i = i+1;
   Target_L_t_vect(i)=dt*D_t_vect*Target_t_vect_long(i+tau_vect_max:-1:i)";
   Mask_L_t_vect(i)=dt*D_t_vect*Mask_t_vect_long(i+tau_vect_max:-1:i)";
end
Target_r_vect_NoThresh = 1000*(r0 + Target_L_x*Target_L_t_vect); %1000 converts to Hz
Target_r_vect = max(Target_r_vect_NoThresh,0); %thresholding
Mask_r_vect_NoThresh = 1000*(r0 + Mask_L_x*Mask_L_t_vect); %1000 converts to Hz
Mask_r_vect = max(Mask_r_vect_NoThresh,0); %thresholding
Both_r_vect_NoThresh = 1000*(r0 + Target_L_x*Target_L_t_vect + Mask_L_x*Mask_L_t_vect); %linear
filter portions add linearly
Both_r_vect = max(Both_r_vect_NoThresh,0); %thresholding
%plot model results
figure(3)
subplot(3,1,1)
contourf(tau_vect,x_vect,D_xt_mat,12); %makes a contour plot of the data
colorbar
set(gca,"Xdir","reverse")
xlabel("tau (ms)")
ylabel("x (deg)")
subplot(3,1,2)
%contourf(t_vect,x_vect,Target_xt_mat); %makes a contour plot of the data
%contourf(t_vect,x_vect,Mask_xt_mat); %makes a contour plot of the data
contourf(t_vect,x_vect,Both_xt_mat); %makes a contour plot of the data
colorbar
xlabel("t (ms)")
ylabel("x (deg)")
subplot(3,1,3)
plot(t_vect,Target_r_vect_NoThresh,"k:")
hold on;
plot(t_vect,Mask_r_vect_NoThresh,"r:")
plot(t_vect,Both_r_vect_NoThresh,":")
p1 = plot(t_vect,Target_r_vect,"k");
p2 = plot(t_vect,Mask_r_vect,"r");
p3 = plot(t_vect,Both_r_vect);
xlabel("time (ms)");
ylabel("rate (Hz)");
legend([p1 p2 p3], "Target only", "Mask only", "Target+Mask")
grid on
hold off;
