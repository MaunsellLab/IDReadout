function path = analysisPath(context,domain,varargin)
% analysisPath  Construct and validate a context-bound analysis path.
%
%   path = analysisPath(context,'IDR','Data','AcrossSessionSummaries')
%
% Path components must be relative names. Absolute components and parent
% traversal are rejected. This function does not create folders.

context = validateAnalysisContext(context);
mustBeTextScalar(domain);
root = string(domainFolder(mfilename('fullpath'),domain,context));
path = root;
for k = 1:numel(varargin)
  mustBeTextScalar(varargin{k});
  component = string(varargin{k});
  if strlength(component)==0 || isAbsolutePath(component) || ...
      any(split(component,filesep)=="..")
    error("IDReadout:UnsafePathComponent", ...
      "Analysis path components must be nonempty relative paths without '..'.");
  end
  path = fullfile(path,component);
end
path = validateAnalysisPath(context,path);
path = char(path);
end

function tf = isAbsolutePath(path)
file = java.io.File(char(path));
tf = file.isAbsolute();
end
