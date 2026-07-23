function path = validateAnalysisPath(context,path)
% validateAnalysisPath  Enforce experimental/synthetic path containment.

context = validateAnalysisContext(context);
mustBeTextScalar(path);
path = canonicalPath(path);

if context.Mode == "synthetic"
  runRoot = canonicalPath(context.RunRoot);
  if ~isWithin(path,runRoot)
    error("IDReadout:SyntheticPathEscape", ...
      "Synthetic analysis attempted access outside its run: %s",path);
  end
else
  projectRoot = canonicalPath(context.ProjectRoot);
  testsRoot = canonicalPath(context.TestsRoot);
  if ~isWithin(path,projectRoot)
    error("IDReadout:ExperimentalPathEscape", ...
      "Experimental analysis attempted access outside the project: %s",path);
  end
  if isWithin(path,testsRoot)
    error("IDReadout:ExperimentalTestsAccess", ...
      "Experimental analysis cannot access synthetic Tests data: %s",path);
  end
end
end

function tf = isWithin(path,root)
tf = path==root || startsWith(path,root+filesep);
end

function path = canonicalPath(path)
path = string(char(java.io.File(char(path)).getCanonicalPath()));
end
