function [index, sideTypeNames] = sideTypeIndex(sideType)
%SIDETYPEINDEX Return canonical index for a named kernel side type.
%
%   index = sideTypeIndex(sideType)
%   [index, sideTypeNames] = sideTypeIndex(sideType)
%   [~, sideTypeNames] = sideTypeIndex()
%
% This is the authoritative list used for indexing kernel side types.

sideTypeNames = {'Diff', 'Change', 'NoChange', 'Left', 'Right', 'RF', 'Opp', 'Chosen', 'NotChosen'};

if nargin < 1 || isempty(sideType)
    index = [];
    return
end
index = find(strcmp(sideTypeNames, sideType), 1);
assert(~isempty(index), 'Unknown sideType "%s". Valid side types are: %s', sideType, strjoin(sideTypeNames, ', '));

end