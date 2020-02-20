function [movieCube, startTimeMS, endTimeMS] = makeBasicMovieCube

%
%  You should design your own experiment -- this is just an example
%  experiment that you can also do with your fly neurons.
%
%  The goal of this example experiment is to test the sensitivity of each
%  fly neuron to motion direction.  That is, do neurons in the fly visual
%  system care about the direction of visual motion?  If they do care, what
%  direction(s) of motion do they like?  Do all neurons like the same
%  direction of motion?
%
%  To try to answer this questions, we need to carefully examine the
%  response of each neuron we record in response to different directions of
%  motion.
%
%  To do this, we need a visual motion stimulus.  In this example, we
%  construct a bar of light (high luminance relative to the background
%  luminance) and we then move that bar across the screen.  When we put the 
%  screen in front of the fly, this means that we will move the bar across
%  retina of the fly.  Note, that this is just one type of motion and you 
%  can think of  all kinds of stimuli that contain what visual scientists 
%  call 'motion energy' (e.g. drifting dots, running dogs, waterfalls, etc.).
%
%  Now that we have decided on a type of motion stimulus, we want to hold
%  that constant and systematically test the effect of motion direction.
%  Each of these tests is refered to as a 'condition'.  In most science
%  experiments we have at least two conditions -- e.g. a treatment
%  condition (new drug) and a control condition (placebo).  In this
%  experiment, we will have 9 conditions:  8 directions of motion (of the
%  bar) and one control condition (no motion, i.e. blank screen).
%
%  Here we build a movie to show a bar moving in 8 directions and to
%  contain some part(s) of the movie where no stimulus is on the screen.
%  Once we show this stimulus to a fly, we can then go back and figure out
%  the neurons response in each of these conditions to learn something
%  about our original questions listed above.
%
%  here we show each bar for 1000 ms.  The total length of the movie is ~10 sec

framesPerSecond = 30;   %% rate at which the movie player will change the image ('frame') on the screen
frameDurationMsec = (1/framesPerSecond)*1000;       %% total time each frame will be shown (msec)


msecBetweenConditions = 200;        %% how long between conditions (msec)
nFramesBetweenConditions = round(msecBetweenConditions/frameDurationMsec);  

msecForEachCondition = 1000;         %% how long each condition will last (msec)
nFramesForEachCondition = round(msecForEachCondition/frameDurationMsec); 

%% stimulus info
BACKGROUND_GRAY = 0;   %% black background
directionsToTest = [ 0 45 90 135 180 225 270 315];   %% the 8 directions to test (in degrees)
numberOfConditions = length(directionsToTest) + 1;   %% all conditions plus a blank condition

%% figure out how many frames we need to show all our conditions
%% (note that we assume some blank frames before and after each condition,
%% so we need the +1 below)
totalFramesNeeded = round((numberOfConditions*nFramesForEachCondition) + ((numberOfConditions+1)*nFramesBetweenConditions)); 

%% size of the movie
maxX = 320;     %% x pixels
maxY = 240;     %% y pixels

%% make a completely 'blank' movie cube
movieCube = ones(maxY,maxX,totalFramesNeeded)*BACKGROUND_GRAY;


lastFrame = 1;      %% just to get going

for condition = 1:numberOfConditions
   
    %% step forward in the movie to leave a blank between conditions
    nextFrame = lastFrame+nFramesBetweenConditions;
   
    firstFrame = nextFrame;                                 %% start at last position in movie
    lastFrame = firstFrame+nFramesForEachCondition-1;       %% the position in the movie where this condition will end
     
    if (condition < numberOfConditions)
        %% drifting bar
        direction = directionsToTest(condition);        
        r = min([maxX maxY]);
        bar_move_pix_per_frame = round(r/nFramesForEachCondition);  %% speed of bar movement  
        %% note, this will make the bar cross the min distance across the screen in the total time
        %%    alotted per condition
        %%  Thus, if the screen subtends ~60 deg of visual angle, and the time per condition is 1000 msec,
        %   then then speed of the bar on teh retina is approximately 60 deg per second.
        barWid = r/10;
        barLength = r/2;    
        barGray = 255;
        movieCube(:,:,firstFrame:lastFrame) = driftingBar(maxY,maxX, nFramesForEachCondition, BACKGROUND_GRAY, direction, barWid, barLength, barGray, bar_move_pix_per_frame);
    else
        %% only the last condition will put us here
        %% do not make any images (this will effectively leave these frames
        %% blank
    end
    
    % let's keep track of our condition start and end times for each condition we create so that we can
    % analyze the data later
    startTimeMS(condition) = (firstFrame-1)*frameDurationMsec;  %% the time of the start of frame 1 is 0 msec, so we need a -1.
    endTimeMS(condition) = lastFrame*frameDurationMsec;
    
end



