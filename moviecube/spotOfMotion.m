function [movieFrames] = spotOfMotion(maxY, maxX, nFramesToMake, backgroundGray, posY, posX, diameter,  angle, spatialPeriod_pixPerCycle, speed_PixPerFrame, contrast );
%
% [movieFrames] = spotOfMotion(maxY, maxX, nFramesToMake, backgroundGray, posY, posX, diameter,  angle, spatialPeriod_pixPerCycle, speed_PixPerFrame, contrast );
%  
% Will make a circular spot of 'motion' (drifting sinusoid) centered at posX, posY.   
% The spot will stay on (and in motion) for the number of indicated nFramesToMake
%
% parameters:
% maxY = number of pixels along row dimension (vertical dimension of each frame)
% maxX = number of pixels along column dimension (horizontal dimension of each frame)
% nFramesToMake  = number of frames of this movie to make. 
% backgroundGray = gray level of the background, range: 0 (black) - 255(white)
% posY = the vertical location of the center of the spot (0 is top)
% posX = the horizontal location of the center of the spot (0 is left)
% diameter = diameter of the spot in pixels
% angle = drift angle of the sinusoid (i.e. the direction of motion)
% spatialPeriod_pixPerCycle = number of pixels covered by each cycle of the
%   sinusoid  (this is basically how 'thin' or 'fat' the lines are)
% speed_PixPerFrame = number of pixels the sinusoid drifts on each frame
% contrast = the depth of modulation of the sinusoid.  (1 = full modulation
%   from black to white, 0 = no modulation = gray screen)
%
% drift angle conventions (deg):
% 0 deg is drift upward 
% 90 deg is drift rightward
% 180 is drift downward
% 270 is drift leftward


angle = angle+90;

[X,Y] = meshgrid([1:maxX],[1:maxY]);

r = diameter/2;
cx = posX;
cy = posY;
Xcentered = (X-cx);
Ycentered = (Y-cy);
D = (sqrt((Xcentered.^2) + (Ycentered.^2)));
z = find(D(:)<=r);
imageApeture = zeros(maxY,maxX);
imageApeture(z) = 1;
backZ = find(imageApeture==0);


clear movieFrames;
%% added March 2007
movieFrames = zeros(maxY,maxX,nFramesToMake,'uint8');

P = ([X*cosd(angle)]+[Y*sind(angle)]);
for frame = 1:nFramesToMake
    phasePix = speed_PixPerFrame*(frame-1);
    Z = sind( ( (P + phasePix)*(360/spatialPeriod_pixPerCycle)) );
    %size(Z)
    %size(imageApeture)
    temp = Z*contrast;     %% range here should be -1 to +1 (*contrast)
    temp = ((temp+1)/2);                 %% [0 - 1]
    temp = temp*255;
    
    %% apply apeture
    temp = imageApeture.*temp;
    temp(backZ) = backgroundGray;
    
    rtemp = round(temp);
    z = find(rtemp(:)>255); rtemp(z) = 255;
    z = find(rtemp(:)<0); rtemp(z) = 0;
    movieFrames(:,:,frame) = rtemp(:,:);
    
    %figure(1); clf; imagesc(temp); pause;
    
end


return;;
