%bdSetValues  Write all values in a pardef/vardef/lagdef cell array.
%Usage:
%   yyydef = bdSetValues(xxxdef,vec)
%where
%   xxxdef is the incoming pardef, vardef or lagdef cell array.
%   vec is an array of new values.
%   yyydef is the returned cell array.
%
%EXAMPLE
%  pardef = [ struct("name","a", "value", 1);
%             struct("name","b", "value",[2,3,4]);
%             struct("name","c", "value",[5 6; 7 8]) ];
%  bdGetValues(pardef)
%
%  ans =
%    1 2 3 4 5 6 7 8
%
%  pardef = bdSetValues(pardef,[8 7 6 5 4 3 2 1]);
%  bdGetValues(pardef)
%
%  ans =
%     8 7 6 5 4 3 2 1
function yyydef = bdSetValues(xxxdef,vec)

    % Verify that vec is the correct size for xxxdef
    if numel(vec) ~= numel(bdGetValues(xxxdef))
        yyydef = [];
        throwAsCaller(MException("bdtoolkit:bdSetValues","Number of new values must match the number of values in xxxdef"));
    end

    % Piecewise copy of vec entries into xxxdef.value entries
    yyydef = xxxdef;
    nvar = numel(yyydef);
    offset = 0;
    for indx=1:nvar
        len = numel(yyydef(indx).value);
        yyydef(indx).value(:) = vec([1:len]+offset);
        offset = offset+len;
    end
end
