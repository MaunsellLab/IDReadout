function context = validateAnalysisContext(context)
% validateAnalysisContext  Validate context schema and containment rules.

if ~isstruct(context) || ~isscalar(context)
  error("IDReadout:InvalidAnalysisContext", ...
    "Analysis context must be a scalar structure.");
end
required = {'SchemaVersion','Mode','RunID','Seed','MaxRunBytes', ...
  'SaveLatents','ProjectRoot','TestsRoot','RunRoot'};
missing = required(~isfield(context,required));
if ~isempty(missing)
  error("IDReadout:InvalidAnalysisContext", ...
    "Analysis context is missing: %s",strjoin(missing,', '));
end
if context.SchemaVersion ~= 1
  error("IDReadout:UnsupportedContextSchema", ...
    "Unsupported analysis-context schema version: %g",context.SchemaVersion);
end

context.Mode = lower(string(context.Mode));
context.RunID = string(context.RunID);
context.ProjectRoot = canonicalPath(context.ProjectRoot);
context.TestsRoot = canonicalPath(context.TestsRoot);
context.RunRoot = string(context.RunRoot);

if ~ismember(context.Mode,["experimental","synthetic"])
  error("IDReadout:InvalidAnalysisMode", ...
    "Context mode must be experimental or synthetic.");
end
expectedTestsRoot = canonicalPath(fullfile(context.ProjectRoot,'Tests'));
if context.TestsRoot ~= expectedTestsRoot
  error("IDReadout:InvalidTestsRoot", ...
    "Context TestsRoot does not match the project Tests folder.");
end
if ~isnumeric(context.MaxRunBytes) || ~isscalar(context.MaxRunBytes) || ...
    ~isfinite(context.MaxRunBytes) || context.MaxRunBytes<=0
  error("IDReadout:InvalidRunQuota","MaxRunBytes must be positive and finite.");
end
if ~islogical(context.SaveLatents) || ~isscalar(context.SaveLatents)
  error("IDReadout:InvalidSaveLatents","SaveLatents must be a logical scalar.");
end

if context.Mode == "experimental"
  if strlength(context.RunID)>0 || strlength(context.RunRoot)>0 || ...
      ~isempty(context.Seed)
    error("IDReadout:InvalidExperimentalContext", ...
      "Experimental contexts cannot contain RunID, RunRoot, or Seed.");
  end
  context.RunRoot = "";
  return
end

if isempty(regexp(char(context.RunID), ...
    '^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$','once'))
  error("IDReadout:InvalidSyntheticRunID","Invalid synthetic RunID.");
end
if ~isnumeric(context.Seed) || ~isscalar(context.Seed) || ...
    ~isfinite(context.Seed) || context.Seed<0 || context.Seed>2^32-1 || ...
    fix(context.Seed)~=context.Seed
  error("IDReadout:InvalidSyntheticSeed","Invalid synthetic random seed.");
end
expectedRunRoot = canonicalPath(fullfile(context.TestsRoot, ...
  'Simulation Runs',context.RunID));
context.RunRoot = canonicalPath(context.RunRoot);
if context.RunRoot ~= expectedRunRoot
  error("IDReadout:InvalidSyntheticRunRoot", ...
    "Synthetic RunRoot is not the canonical folder for this RunID.");
end
end

function path = canonicalPath(path)
path = string(char(java.io.File(char(path)).getCanonicalPath()));
end
