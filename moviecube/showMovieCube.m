function showMovieCube(movieCube, fps)
% SHOWMOVIECUBE Play a movie interactively
%
% This is a wrapper around mplay which lets us play a movie
% cube. Our movieCube needs to be of the form
% 
%     M(y, x, frame)
%


mplay(movieCube, fps)
