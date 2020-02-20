function [movieFrames] = spotOfLight(maxY, maxX, nFramesToMake, backgroundGray, posY, posX, diameter, spotGray )
%
% [movieFrames] = spotOfLight(maxY, maxX, nFramesToMake, backgroundGray,posX, posY, diameter, spotGray )
%  
% Will make a circular spot of 'light'  (can be black, white, or any gray level)
% The spot will stay on for the number of indicated nFramesToMake
%
% parameters:
% maxY = number of pixels along row dimension (vertical dimension of each frame)
% maxX = number of pixels along column dimension (horizontal dimension of each frame)
% nFramesToMake  = number of frames of this movie to make. 
% backgroundGray = gray level of the background, range: 0 (black) - 255(white)
% posY = the vertical location of the center of the spot (0 is top)
% posX = the horizontal location of the center of the spot (0 is left)
% diameter = diameter of the spot in pixels
% spotGray = gray level of the spot (0 (black) = 255 (white))



[X,Y] = meshgrid([1:maxX],[1:maxY]);

r = diameter/2;
cx = posX;
cy = posY;
Xcentered = (X-cx);
Ycentered = (Y-cy);
D = (sqrt((Xcentered.^2) + (Ycentered.^2)));
z = find(D(:)<=r);
imageArray = ones(maxY,maxX)*backgroundGray;
imageArray(z) = spotGray;

clear movieFrames;
%% added March 2007
movieFrames = zeros(maxY,maxX,nFramesToMake,'uint8');

for frame = 1:nFramesToMake 
   
    rtemp = round(imageArray);
    z = find(rtemp(:)>255); rtemp(z) = 255;
    z = find(rtemp(:)<0); rtemp(z) = 0;
    movieFrames(:,:,frame) = rtemp(:,:);
  
end
