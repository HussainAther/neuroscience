
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
