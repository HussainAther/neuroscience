function [movieFrames] = driftingBar(maxY, maxX, nFramesToMake, backgroundGray, angle, barWidthPix, barLengthPix, barGray, speed_PixPerFrame)
%
% [movieFrames] = driftingBar(maxY, maxX, nFramesToMake, backgroundGray, angle, barWidthPix, barLengthPix, barGray, speed_PixPerFrame)
%  
% Will make a single bar that drifts across the screen in a direction indicated by 'angle'.
% The bar length is perpendicular to the direction of drift.  The bar is
% windowed by a circular apeture of diameter min(maxX,maxY).  The output values in each frame
% range from 0 (black) to 255 (white)
%
% parameters:
% maxY = number of pixels along row dimension (vertical dimension of each frame)
% maxX = number of pixels along column dimension (horizontal dimension of each frame)
% nFramesToMake  = number of frames of this movie to make. 
% backgroundGray = gray level of the background, range: 0 (black) - 255(white)
% angle = drift angle of the bar (i.e. the direction of motion)
% barWidthPix = the width of the bar in pixels (min = 1)
% barLengthPix = teh length of the bar in pixels (min = 1)
% speed_PixPerFrame = number of pixels the bar drifts on each frame
%
% drift angle conventions (deg):
% 0 deg is drift upward 
% 90 deg is drift rightward
% 180 is drift downward
% 270 is drift leftward

angle = -(angle+180);

if ((barGray<0) | (barGray > 255)) errore('Bar gray value out of range (0-255)');end;
if ((backgroundGray<0) | (backgroundGray > 255)) errore('Background gray value out of range (0-255)');end;

[X,Y] = meshgrid([1:maxX],[1:maxY]);

r = min([maxX maxY])/2;
cx = round(maxX/2);
cy = round(maxY/2);
Xcentered = (X-cx);
Ycentered = (Y-cy);
D = (sqrt((Xcentered.^2) + (Ycentered.^2)));
z = find(D(:)<=(round(r)-1));
imageApeture = zeros(maxY,maxX);
imageApeture(z) = 1;
backZ = find(imageApeture==0);


imageArrayNull = ones(maxY,maxX)*backgroundGray;
imageSize = size(imageArrayNull);

velocityX_pixPerFrame = cos(angle)*speed_PixPerFrame;
velocityY_pixPerFrame = sin(angle)*speed_PixPerFrame;

clear movieFrames;
%% added March 2007
movieFrames = zeros(maxY,maxX,nFramesToMake,'uint8');


%% startlocation of the center point
angleRads = (angle/180)*pi;
startAngleRads = (angleRads+pi);        %% flip 180 degs
barCenterLastX = cy + (cos(startAngleRads)*r);
barCenterLastY = cx + (sin(startAngleRads)*r);
barCenterLast = [barCenterLastX barCenterLastY];

car = cos(angleRads);
sar = sin(angleRads);
car2 = cos(angleRads-(pi/2));
sar2 = sin(angleRads-(pi/2));
car3 = cos(angleRads+(pi/2));
sar3 = sin(angleRads+(pi/2));

if (abs(car)<0.001) car = 0; end;
if (abs(sar)<0.001) sar = 0; end;
if (abs(car2)<0.001) car2 = 0; end;
if (abs(sar2)<0.001) sar2 = 0; end;
if (abs(car3)<0.001) car3 = 0; end;
if (abs(sar3)<0.001) sar3 = 0; end;

offsetX = cosd(angle)*speed_PixPerFrame;
offsetY = sind(angle)*speed_PixPerFrame;



for frame = 1:nFramesToMake

    %% find index values not that are part of the bar
    imageArray = imageArrayNull; 
    
    %% move the bar center (new frame)
    barCenter = barCenterLast + [offsetX offsetY];
    
    %d = sqrt( (barCenter-[cx cy]).^2);
    %if (d>(r+1)) 
        %% do nothing, bar is outside apeture
    %   break; end;
    
    for w = -(round(barWidthPix/2)-1):0.5:(round(barWidthPix/2)-1)
    
        barCenterCurrent = barCenter + [ car*w sar*w  ];  
    
        ind = round(barCenterCurrent);
        z = find(ind<=0);
        zz = find(ind>imageSize);
        if ((length(z)==0)&(length(zz)==0))
            imageArray(ind(1),ind(2)) = 1;
        end
        
        for l = 0:0.5:min([r (barLengthPix/2)])
            dX = car2*l;
            dY = sar2*l;
            pix = barCenterCurrent + [dX dY];
            ind = round(pix);
            z = find(ind<=0);
            zz = find(ind>imageSize);
            if ((length(z)==0)&(length(zz)==0))
                imageArray(ind(1),ind(2)) = 1;
            end
        end
        for l = 0:0.5:min([r ((barLengthPix/2)-1)])
            dX = car3*l;
            dY = sar3*l;
            pix = barCenterCurrent + [dX dY];
            ind = round(pix);
            z = find(ind<=0);
            zz = find(ind>imageSize);
            if ((length(z)==0)&(length(zz)==0))
                imageArray(ind(1),ind(2)) = 1;
            end
        end
        
    end
    
    temp = imageArray.*imageApeture;
    z = find(temp(:) == 1);
    temp(z) = barGray;
    temp(backZ) = backgroundGray;       %% outside of apeture 
    
    movieFrames(:,:,frame) = temp;

    barCenterLast = barCenter;
    
end




