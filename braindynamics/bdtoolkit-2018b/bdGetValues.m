%bdGetValues  Read all values from a pardef/vardef/lagdef cell array. 
%   Returns the contents of all xxxdef.value fields as single monolithic
%   column vector
%Usage:
%   Y = bdGetValues(xxxdef)
%where
%   xxxdef is any array of structs that contains a field called 'value', 
%       such as sys.vardef, sys.pardef and sys.lagdef.
%
%EXAMPLE
%   vardef = [ struct('name','a', 'value',1);
%              struct('name','b', 'value',[2 3 4]);
%              struct('name','c', 'value',[5 7 9 11; 6 8 10 12]);
%              struct('name','d', 'value',13); ];
%   Y0 = bdGetValues(vardef)
%
%   Returns Y0 as the column vector [1:13]'
function vec = bdGetValues(xxxdef)
    % extract the value fields of vardef as a cell array
    vec = {xxxdef.value}';

    % convert each cell entry to a column vector
    for indx=1:numel(vec)
        vec{indx} = reshape(vec{indx}, [], 1);
    end

    % concatenate the column vectors to a simple vector
    vec = cell2mat(vec);
end
