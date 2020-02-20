%bdSetValue  Write a value in a pardef/vardef/lagdef cell array.
%Usage:
%   yyydef = bdSetValue(xxxdef,'name',val)
%where
%   xxxdef is the incoming pardef, vardef or lagdef cell array.
%   name is the string name of the element to be updated.
%   val is the new value to be applied.
%   yyydef is the returned cell array.
%
%EXAMPLE
%  pardef = [ struct('name','a', 'value', 1);
%             struct('name','b', 'value',[2,3,4]);
%             struct('name','c', 'value',[5 6; 7 8]) ];
%  pardef = bdSetValue(pardef,'b',[3 6 9]);
%  bdGetValue(pardef,'b')
%
%  ans =
%     3     6     9
function yyydef = bdSetValue(xxxdef,name,val)
    yyydef = xxxdef;
    nvar = numel(yyydef);
    for indx=1:nvar
        if strcmp(yyydef(indx).name,name)==1
            yyydef(indx).value = val;
            return
        end
    end
    warning([name, ' not found in xxxdef']);
end
