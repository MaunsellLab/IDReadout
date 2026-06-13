function domainRoot = domainFolder(callerFullPath)
% domainFolder  Return the top-level domain containing the caller.
%
% Example:
%   domainRoot = domainFolder(mfilename('fullpath'));
%
% A caller anywhere under IDReadout/IDR returns IDReadout/IDR.
% A caller anywhere under IDReadout/IDQ returns IDReadout/IDQ.

arguments
  callerFullPath {mustBeTextScalar}
end

project = matlab.project.rootProject;

if isempty(project)
  error("IDReadout:ProjectNotOpen", ...
    "The IDReadout MATLAB Project is not open.");
end

projectRoot = string(project.RootFolder);
callerFullPath = string(callerFullPath);
callerFolder = string(fileparts(callerFullPath));

relativeFolder = extractAfter( ...
  callerFolder, projectRoot + filesep);

parts = split(relativeFolder, filesep);

if isempty(parts) || parts(1) == ""
  error("IDReadout:CallerOutsideProject", ...
    "Caller is not inside the IDReadout project: %s", ...
    callerFullPath);
end

domainRoot = fullfile(projectRoot, parts(1));

% Optional protection against calls from Utilities or another root folder.
validDomains = ["IDR", "IDQ"];

if ~ismember(parts(1), validDomains)
  error("IDReadout:InvalidDomain", ...
    "'%s' is not a recognized analysis domain.", parts(1));
end
end