function cleanupObj = initProjectPath()
% initProjectPath  Add this folder and all subfolders to the MATLAB path.
%
% Usage:
%   cleanupObj = initProjectPath();
%
% Only folders added by this call are removed when cleanupObj is cleared.
% This makes nested calls safe.

projectRoot = fileparts(mfilename('fullpath'));

% All folders requested for this project
projectFolders = strsplit(genpath(projectRoot), pathsep);
projectFolders = projectFolders(~cellfun('isempty', projectFolders));

% Determine which folders are not already on the MATLAB path
currentFolders = strsplit(path, pathsep);
foldersAdded = projectFolders(~ismember(projectFolders, currentFolders));

% Add only previously absent folders
if ~isempty(foldersAdded)
    addpath(foldersAdded{:});
end

% Remove only folders added by this call
if nargout > 0
    cleanupObj = onCleanup(@() removeAddedFolders(foldersAdded));
end
end


function removeAddedFolders(foldersAdded)
% Remove folders still present, avoiding warnings if already removed.

if isempty(foldersAdded)
    return
end

currentFolders = strsplit(path, pathsep);
foldersToRemove = foldersAdded(ismember(foldersAdded, currentFolders));

if ~isempty(foldersToRemove)
    rmpath(foldersToRemove{:});
end
end