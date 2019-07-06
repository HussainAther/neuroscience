function schema
%SCHEMA Schema for PLAYER class
  
%   Copyright 2003 The MathWorks, Inc.
%   $Revision: $  $Date: $

package = findpackage('mmovie');
thisclass = schema.class(package,'player');

p = schema.prop(thisclass,'hfig','MATLAB array');
p = schema.prop(thisclass,'fcns','MATLAB array');
