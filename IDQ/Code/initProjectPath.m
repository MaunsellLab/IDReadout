function cleanupObj = initProjectPath()
% initProjectPath  Add this folder and all subfolders to the MATLAB path.
%
% Usage:
%   cleanupObj = initProjectPath();   % temporary, removed when cleanupObj clears

projectRoot = fileparts(mfilename('fullpath'));

% Add root and all nested subfolders
projectPath = genpath(projectRoot);
addpath(projectPath);

% return cleanup object so the path change is temporary
if nargout > 0
  cleanupObj = onCleanup(@() rmpath(projectPath));
end
end