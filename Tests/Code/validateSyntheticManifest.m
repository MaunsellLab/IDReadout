function manifest = validateSyntheticManifest(context)
% validateSyntheticManifest  Verify that a run has matching provenance.

context = validateAnalysisContext(context);
if context.Mode ~= "synthetic"
  error("IDReadout:SyntheticContextRequired", ...
    "Manifest validation requires a synthetic context.");
end
path = analysisPath(context,'Tests','manifest.mat');
if ~isfile(path)
  error("IDReadout:MissingSyntheticManifest", ...
    "Synthetic run has no manifest: %s",context.RunRoot);
end
S = load(path,'manifest');
if ~isfield(S,'manifest') || ~isstruct(S.manifest)
  error("IDReadout:InvalidSyntheticManifest", ...
    "Synthetic manifest is missing or invalid: %s",path);
end
manifest = S.manifest;
required = {'ManifestSchemaVersion','DataOrigin','RunID','RandomSeed', ...
  'MaxRunBytes','RunRootAtCreation'};
if ~all(isfield(manifest,required)) || ...
    manifest.ManifestSchemaVersion~=1 || ...
    string(manifest.DataOrigin)~="synthetic" || ...
    string(manifest.RunID)~=context.RunID || ...
    manifest.RandomSeed~=context.Seed || ...
    canonicalPath(manifest.RunRootAtCreation)~=canonicalPath(context.RunRoot)
  error("IDReadout:SyntheticManifestMismatch", ...
    "Synthetic manifest does not match the supplied context.");
end
end

function path = canonicalPath(path)
path = string(char(java.io.File(char(path)).getCanonicalPath()));
end
