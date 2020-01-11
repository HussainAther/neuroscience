
%{
Model lateral geniculate nucleus (LGN) neuron that receives a target bar
of light preceded by two masking bars of light. If positioned in space
and time appropriately, the target bar should be less visible.
}%

clear all;
close all;

% Initialize parameters.
sigma_c = 1.4; % width of center portion of spatial r.f. [degrees]
sigma_s = 2.1; % determines width of surround portion of spatial r.f. [deg]
A_c = 1; % strength of center portion of spatial r.f.
A_s = 0.9; % strength of surround portion of spatial r.f.
dx = 0.05; % resolution of spatial grid
x_min = -4*sigma_s; % x-value roughly giving -x border of r.f.
x_max = 4*sigma_s; % x-value roughly giving +x border of r.f.
x_vect = x_min:dx:x_max; % values of x over which to compute D_x
alpha = 1/10; % determines length of temporal filter [ms^-1]
dt = 1; % ms
tau_vect_max = 20/alpha; % tau value giving border of temporal r.f.
tau_vect = 0:dt:tau_vect_max; % value of tau over which to compute D_t

% Define and then plot spatial kernel D_x.
D_x_center = A_c*exp(-(x_vect.^2)/(2*(sigma_c^2)))/sqrt(2*pi*sigma_c^2);
D_x_surround = A_s*exp(-(x_vect.^2)/(2*(sigma_s^2)))/sqrt(2*pi*sigma_s^2);
D_x_vect = D_x_center - D_x_surround;
figure(1)
subplot(3,1,1)
%plot(x_vect, D_x_center, "r--") % plot center with red dashed lines
hold on
%plot(x_vect, D_x_surround, "k--") % plot surround with black dashed lines
plot(x_vect, D_x_vect) % plot full receptive field with solid blue line
xlabel("x (deg)")
ylabel("D_x")
