function cache = buildCache(filePaths, readerFcn, outputPath, varargin)
%BUILDCACHE Combine project-specific session records into cache v2.
%
% readerFcn(path) must return a scalar cache fragment satisfying the same
% field contract as a full cache. The fragment may omit schemaVersion and
% provenance. This function concatenates fragments without interpreting raw
% task fields.

p = inputParser;
p.addParameter('Replace', false, @(x)islogical(x) && isscalar(x));
p.addParameter('Provenance', struct(), @isstruct);
p.parse(varargin{:});
opts = p.Results;

if isfile(outputPath) && ~opts.Replace
  error('maxpoolfit:buildCache:Exists', ...
    'Output exists. Set Replace=true to overwrite: %s', outputPath);
end
if isstring(filePaths), filePaths = cellstr(filePaths); end
assert(iscell(filePaths) && ~isempty(filePaths), ...
  'filePaths must be a nonempty cell array or string array.');
assert(isa(readerFcn, 'function_handle'), 'readerFcn must be a function handle.');

parts = cell(numel(filePaths),1);
for k = 1:numel(filePaths)
  parts{k} = readerFcn(filePaths{k});
  maxpoolfit.validateCache(partialToCache(parts{k}), 'AllowEmpty', false);
end

cache = parts{1};
audits=cell(numel(parts),1);
if isfield(parts{1},'extractionAudit'), audits{1}=parts{1}.extractionAudit; end
for k = 2:numel(parts)
  cache.trials = [cache.trials; parts{k}.trials]; %#ok<AGROW>
  cache.streamDirDeg = [cache.streamDirDeg; parts{k}.streamDirDeg]; %#ok<AGROW>
  cache.stepCohPC = [cache.stepCohPC; parts{k}.stepCohPC]; %#ok<AGROW>
  cache.preCohPC = [cache.preCohPC; parts{k}.preCohPC]; %#ok<AGROW>
  cache.driftDirDeg = [cache.driftDirDeg; parts{k}.driftDirDeg]; %#ok<AGROW>
  cache.preDriftCohPC = [cache.preDriftCohPC; parts{k}.preDriftCohPC]; %#ok<AGROW>
  if isfield(cache,'raw') || isfield(parts{k},'raw')
    assert(isfield(cache,'raw') && isfield(parts{k},'raw'), ...
      'All cache fragments must either contain raw fields or omit them.');
    cache.raw = appendRaw(cache.raw, parts{k}.raw);
  end
  if isfield(parts{k},'extractionAudit'), audits{k}=parts{k}.extractionAudit; end
end

cache.schemaVersion = "MaxPoolingTrialCache-v2";
cache.createdAt = datetime('now');
cache.createdBy = mfilename;
cache.sourceFiles = string(filePaths(:));
cache.provenance = opts.Provenance;
cache.extractionAudits=audits(~cellfun(@isempty,audits));
if isfield(cache,'extractionAudit'), cache=rmfield(cache,'extractionAudit'); end
maxpoolfit.validateCache(cache);
save(outputPath, 'cache', '-v7.3');
end

function cache = partialToCache(part)
cache = part;
cache.schemaVersion = "MaxPoolingTrialCache-v2";
end

function out = appendRaw(a,b)
assert(isequal(sort(fieldnames(a)),sort(fieldnames(b))), ...
  'Raw cache fragments have different fields.');
out = a;
names = fieldnames(a);
for k = 1:numel(names)
  f = names{k};
  if iscell(a.(f)) || isstring(a.(f)) || istable(a.(f)) || ...
      (isnumeric(a.(f)) && size(a.(f),1) > 1)
    out.(f) = [a.(f); b.(f)];
  else
    assert(isequaln(a.(f),b.(f)), 'Raw metadata field %s differs.', f);
  end
end
end
