%bdGetValue  Read a named value from a pardef/vardef/lagdef cell array.
%Usage:
%   [val,idx] = bdGetValue(xxxdef,'name')
%where
%   xxxdef is the incoming pardef, vardef, lagdef or auxdef cell array.
%   name is the string name of the element to be updated.
%   val is the corresponding value found in the cell array.
%   idx is the corresponding index of the relevant cell array.
%   Both val and indx are returned empty if no matching name was found.
%
%EXAMPLE
%  pardef = [ struct('name','a', 'value', 1);
%             struct('name','b', 'value',[2,3,4]);
%             struct('name','c', 'value',[5 6; 7 8]) ];
%  val = bdGetValue(pardef,'b')
%
%  val =
%     2     3     4
function [val,idx] = bdGetValue(xxxdef,name)
    val = [];
    idx = [];
    nvar = numel(xxxdef);
    for indx=1:nvar
        if strcmp(xxxdef(indx).name,name)==1
            val = xxxdef(indx).value;
            idx = indx;
            return
        end
    end
    %warning('bdGetValue() failed to find a matching name');
end
