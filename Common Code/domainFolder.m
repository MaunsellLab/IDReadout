function domainRoot = domainFolder(callerFullPath, requestedDomain)
% domainFolder  Return an inferred or explicitly requested project domain.
%
% Examples:
%   domainFolder(mfilename('fullpath'))
%   domainFolder(mfilename('fullpath'), 'IDQ')
%   domainFolder(mfilename('fullpath'), 'IDR')

arguments
  callerFullPath {mustBeTextScalar}
  requestedDomain {mustBeTextScalar} = ""
end

project = matlab.project.rootProject;

if isempty(project)
  error("IDReadout:ProjectNotOpen", ...
    "The IDReadout MATLAB Project is not open.");
end

projectRoot = string(project.RootFolder);
callerFullPath = string(callerFullPath);
callerFolder = string(fileparts(callerFullPath));

if callerFolder ~= projectRoot && ...
    ~startsWith(callerFolder, projectRoot + filesep)
  error("IDReadout:CallerOutsideProject", ...
    "Caller is not inside the IDReadout project: %s", ...
    callerFullPath);
end

validDomains = ["IDR", "IDQ", "Common Code"];

if strlength(string(requestedDomain)) > 0
  domainName = string(requestedDomain);
else
  relativeFolder = extractAfter( ...
    callerFolder, projectRoot + filesep);
  parts = split(relativeFolder, filesep);

  if isempty(parts) || parts(1) == ""
    error("IDReadout:CallerOutsideDomain", ...
      "Caller is not inside a project domain: %s", callerFullPath);
  end

  domainName = parts(1);
end

match = find(strcmpi(domainName, validDomains), 1);
if isempty(match)
  error("IDReadout:InvalidDomain", ...
    "'%s' is not a recognized project domain.", domainName);
end

% Use canonical capitalization.
domainName = validDomains(match);
domainRoot = fullfile(projectRoot, domainName);
end