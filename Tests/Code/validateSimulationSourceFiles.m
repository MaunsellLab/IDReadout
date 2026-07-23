function records = validateSimulationSourceFiles(context,sourceFiles)
% validateSimulationSourceFiles  Require exact, unchanged manifest sources.
%
% Synthetic analyses may read experimental inputs only when every input was
% declared in initializeSyntheticRun. Both canonical path and SHA-256 must
% match the manifest, preventing an analysis from silently mixing sources.

manifest = validateSyntheticManifest(context);
sourceFiles = string(sourceFiles(:));
if isempty(sourceFiles)
  records = struct([]);
  return
end
if ~isfield(manifest,'SourceFiles') || isempty(manifest.SourceFiles)
  error('IDReadout:SimulationSourceNotDeclared', ...
    'The simulation manifest declares no experimental source files.');
end

declared = manifest.SourceFiles;
declaredPaths = strings(numel(declared),1);
for k = 1:numel(declared)
  declaredPaths(k) = canonicalPath(declared(k).Path);
end

records = repmat(struct('Path',"",'SHA256',""),numel(sourceFiles),1);
for k = 1:numel(sourceFiles)
  if ~isfile(sourceFiles(k))
    error('IDReadout:MissingSimulationSource', ...
      'Simulation source file does not exist: %s',sourceFiles(k));
  end
  path = canonicalPath(sourceFiles(k));
  row = find(declaredPaths == path);
  if numel(row) ~= 1
    error('IDReadout:SimulationSourceNotDeclared', ...
      'Source was not uniquely declared in the simulation manifest: %s',path);
  end
  digest = fileSHA256(path);
  if ~isfield(declared(row),'SHA256') || ...
      lower(string(declared(row).SHA256)) ~= digest
    error('IDReadout:SimulationSourceChanged', ...
      'Source SHA-256 differs from the simulation manifest: %s',path);
  end
  records(k).Path = path;
  records(k).SHA256 = digest;
end
end

function digest = fileSHA256(path)
md = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(path,'rb');
if fid < 0
  error('IDReadout:SourceReadFailed','Could not read source file: %s',path);
end
cleanup = onCleanup(@() fclose(fid));
while ~feof(fid)
  bytes = fread(fid,1024*1024,'*uint8');
  if ~isempty(bytes), md.update(typecast(bytes,'int8')); end
end
raw = typecast(md.digest(),'uint8');
digest = lower(string(reshape(dec2hex(raw,2).',1,[])));
clear cleanup
end

function path = canonicalPath(path)
path = string(char(java.io.File(char(path)).getCanonicalPath()));
end
