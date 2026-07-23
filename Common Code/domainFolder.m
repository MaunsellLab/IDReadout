function domainRoot = domainFolder(callerFullPath,requestedDomain,context)
% domainFolder  Resolve an experimental domain or isolated simulation root.
%
% Experimental, backward-compatible forms:
%   root = domainFolder(mfilename('fullpath'));
%   root = domainFolder(mfilename('fullpath'),'IDR');
%   root = domainFolder(mfilename('fullpath'),'IDQ');
%   root = domainFolder(mfilename('fullpath'),'Common Code');
%
% Context-aware forms:
%   root = domainFolder(mfilename('fullpath'),'IDR',context);
%   runRoot = domainFolder(mfilename('fullpath'),'Tests',context);
%
% With an experimental context, IDR, IDQ, and Common Code resolve beneath
% the project root. With a synthetic context they resolve beneath the named
% simulation run. Tests can only be requested with a validated synthetic
% context and resolves to that run's root, never to the unrestricted Tests
% folder.

if nargin < 2
  requestedDomain = "";
end
if nargin < 3
  context = [];
end

mustBeTextScalar(callerFullPath);
mustBeTextScalar(requestedDomain);
requestedDomain = string(requestedDomain);

project = matlab.project.rootProject;
if isempty(project)
  error("IDReadout:ProjectNotOpen", ...
    "The IDReadout MATLAB Project is not open.");
end
projectRoot = canonicalPath(project.RootFolder);

if isempty(context)
  mode = "experimental";
else
  context = validateAnalysisContext(context);
  if canonicalPath(context.ProjectRoot) ~= projectRoot
    error("IDReadout:ContextProjectMismatch", ...
      "The context belongs to a different MATLAB Project.");
  end
  mode = string(context.Mode);
end

if requestedDomain == ""
  if mode == "synthetic"
    error("IDReadout:SyntheticDomainRequired", ...
      "Synthetic path resolution requires an explicit domain.");
  end
  requestedDomain = inferCallerDomain(callerFullPath,projectRoot);
end

validDomains = ["IDR","IDQ","Common Code"];
if requestedDomain == "Tests"
  if isempty(context) || mode ~= "synthetic"
    error("IDReadout:TestsContextRequired", ...
      "Tests can only be resolved with a validated synthetic context.");
  end
  domainRoot = char(canonicalPath(context.RunRoot));
  return
end
if ~ismember(requestedDomain,validDomains)
  error("IDReadout:InvalidDomain", ...
    "'%s' is not a recognized analysis domain.",requestedDomain);
end

if mode == "experimental"
  domainRoot = char(canonicalPath(fullfile(projectRoot,requestedDomain)));
else
  domainRoot = char(canonicalPath(fullfile(context.RunRoot,requestedDomain)));
end
end

function domain = inferCallerDomain(callerFullPath,projectRoot)
callerFolder = canonicalPath(fileparts(string(callerFullPath)));
prefix = projectRoot + filesep;
if callerFolder == projectRoot || ~startsWith(callerFolder,prefix)
  error("IDReadout:CallerOutsideProject", ...
    "Caller is not inside the IDReadout project: %s",callerFullPath);
end
relativeFolder = extractAfter(callerFolder,prefix);
parts = split(relativeFolder,filesep);
domain = parts(1);
if ~ismember(domain,["IDR","IDQ","Common Code"])
  error("IDReadout:InvalidDomain", ...
    "'%s' is not a recognized analysis domain.",domain);
end
end

function path = canonicalPath(path)
path = string(char(java.io.File(char(path)).getCanonicalPath()));
end
