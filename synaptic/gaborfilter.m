%{
the receptive fields of simple cells in V1 reflect the orientation and spatial frequency
preference of the neurons. One way to model this is to use the Gabor function, which is
basically a two-dimensional Gaussian modulated by a sinusoid.
%}
% gabor_filter.m
function f = gabor_filter(OR, SF)
% Creates a Gabor filter for orientation and spatial frequency
% selectivity of orientation OR (in radians) and spatial frequency SF.
%
% set parameters
sigma_x=7;% standard deviation of 2D Gaussian along x-dir
sigma_y=17;% standard deviation of 2D Gaussian along y-dir
%
% create filter
[x,y]=meshgrid(-20:20);
X=x*cos(OR)+y*sin(OR); %rotate axes
Y=-x*sin(OR)+y*cos(OR);
f=(1/(2*pi*sigma_x*sigma_y)).*exp(-(1/2)*(((X/sigma_x).^2)+ ...
    ((Y/sigma_y).^2))).*sin(2*pi*SF*X);
%gabor_conv.m
clear all; close all
I=imread('rose.jpg');
OR=0; SF=.01;
G=gabor_filter(OR,SF);
figure
subplot(1,3,1); imagesc(G); axis square; colorbar; title ('Gabor function')
subplot(1,3,2); imagesc(I); title('original image')
subplot(1,3,3); imagesc(conv2_mirrored(double(I),G));
colormap(bone); title(['Convolved image OR=',num2str(OR),' SF=', num2str(SF)])
