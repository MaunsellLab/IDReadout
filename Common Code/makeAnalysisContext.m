function context = makeAnalysisContext(mode,varargin)
% makeAnalysisContext  Define an explicit experimental or synthetic run.
%
%   context = makeAnalysisContext("experimental")
%
%   context = makeAnalysisContext("synthetic", ...
%     'RunID',"MTMax_20260721_001",'Seed',1729)
%
% Synthetic contexts require a RunID and random seed. Creating a context
% has no file-system side effects; initializeSyntheticRun creates the run.

mustBeTextScalar(mode);
mode = lower(string(mode));
if ~ismember(mode,["experimental","synthetic"])
  error("IDReadout:InvalidAnalysisMode", ...
    "Mode must be 'experimental' or 'synthetic'.");
end

p = inputParser;
addParameter(p,'RunID',"",@(x) ischar(x)||isstring(x));
addParameter(p,'Seed',[],@(x) isempty(x) || ...
  (isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=2^32-1&&fix(x)==x));
addParameter(p,'MaxRunBytes',2e9,@(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'SaveLatents',false,@(x) islogical(x)&&isscalar(x));
parse(p,varargin{:});
opts = p.Results;

project = matlab.project.rootProject;
if isempty(project)
  error("IDReadout:ProjectNotOpen", ...
    "The IDReadout MATLAB Project is not open.");
end
projectRoot = canonicalPath(project.RootFolder);
testsRoot = canonicalPath(fullfile(projectRoot,'Tests'));

runID = string(opts.RunID);
if mode == "synthetic"
  if strlength(runID)==0
    error("IDReadout:SyntheticRunIDRequired", ...
      "Synthetic contexts require a nonempty RunID.");
  end
  if isempty(regexp(char(runID),'^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$','once'))
    error("IDReadout:InvalidSyntheticRunID", ...
      "RunID must contain only letters, digits, period, underscore, or hyphen.");
  end
  if isempty(opts.Seed)
    error("IDReadout:SyntheticSeedRequired", ...
      "Synthetic contexts require an explicit random seed.");
  end
  runRoot = canonicalPath(fullfile(testsRoot,'Simulation Runs',runID));
else
  if strlength(runID)>0 || ~isempty(opts.Seed)
    error("IDReadout:ExperimentalRunOptions", ...
      "RunID and Seed apply only to synthetic contexts.");
  end
  runRoot = "";
end

context = struct();
context.SchemaVersion = 1;
context.Mode = mode;
context.RunID = runID;
context.Seed = opts.Seed;
context.MaxRunBytes = double(opts.MaxRunBytes);
context.SaveLatents = opts.SaveLatents;
context.ProjectRoot = projectRoot;
context.TestsRoot = testsRoot;
context.RunRoot = runRoot;
context = validateAnalysisContext(context);
end

function path = canonicalPath(path)
path = string(char(java.io.File(char(path)).getCanonicalPath()));
end
