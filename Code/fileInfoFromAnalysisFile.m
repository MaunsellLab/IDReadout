function fileInfoRow = fileInfoFromAnalysisFile(filePath)
% fileInfoFromAnalysisFile  Return one-row metadata table for an analysis MAT file.
%
% Probe-session derived files are identified by the presence of the
% sessionProbeHeader variable, not by filename conventions.
%
% Parent acquisition context is represented by:
%   parentNProbeDirections
%   parentProbeDirectionsDeg
%
% Do not use ambiguous fields such as isProbeSpecific.  A file either has a
% sessionProbeHeader or it does not; parent probe-direction context is a
% separate property.

filePath = char(filePath);
[folder, fileName, ext] = fileparts(filePath);
fileNameExt = [fileName ext];

vars = whos('-file', filePath);
varNames = {vars.name};
hasSessionHeader = any(strcmp(varNames, 'sessionHeader'));
hasSessionProbeHeader = any(strcmp(varNames, 'sessionProbeHeader'));

S = struct();
loadVars = {};
if hasSessionHeader
  loadVars{end+1} = 'sessionHeader'; %#ok<AGROW>
end
if hasSessionProbeHeader
  loadVars{end+1} = 'sessionProbeHeader'; %#ok<AGROW>
end
if ~isempty(loadVars)
  S = load(filePath, loadVars{:});
end

sessionHeader = struct();
sessionProbeHeader = struct();
if isfield(S, 'sessionHeader')
  sessionHeader = S.sessionHeader;
end
if isfield(S, 'sessionProbeHeader')
  sessionProbeHeader = S.sessionProbeHeader;
end

probeDirDeg = NaN;
probeTag = '';
parentNProbeDirections = NaN;
parentProbeDirectionsDeg = [];
parentIsSingleProbe = false;
parentIsInterleavedProbe = false;
prefCohNoisePC = NaN;
probeCohNoisePC = NaN;
nTrials = NaN;
nNoiseTrials = NaN;
parentFileName = '';

if hasSessionProbeHeader
  probeDirDeg = localNumericField(sessionProbeHeader, 'probeDirDeg');
  probeTag = localCharField(sessionProbeHeader, 'probeTag');
  parentNProbeDirections = localNumericField(sessionProbeHeader, 'parentNProbeDirections');
  parentProbeDirectionsDeg = localNumericVectorField(sessionProbeHeader, 'parentProbeDirectionsDeg');
  probeCohNoisePC = localNumericField(sessionProbeHeader, 'probeCohNoisePC');
  nTrials = localNumericField(sessionProbeHeader, 'nTrials');
  nNoiseTrials = localNumericField(sessionProbeHeader, 'nNoiseTrials');
  parentFileName = localCharField(sessionProbeHeader, 'parentFileName');
  prefCohNoisePC = localNumericField(sessionProbeHeader, 'prefCohNoisePC');
end

% If this is a parent/session-level file rather than a probe-session file,
% use sessionHeader for parent acquisition context when available.
if isnan(parentNProbeDirections) && hasSessionHeader
  parentNProbeDirections = localNumericField(sessionHeader, 'nProbeDirections');
end
if isempty(parentProbeDirectionsDeg) && hasSessionHeader
  parentProbeDirectionsDeg = localNumericVectorField(sessionHeader, 'probeDirectionsDeg');
end

if isfinite(parentNProbeDirections)
  parentIsSingleProbe = parentNProbeDirections == 1;
  parentIsInterleavedProbe = parentNProbeDirections > 1;
end

if isnan(prefCohNoisePC) && hasSessionHeader
  prefCohNoisePC = localNumericField(sessionHeader, 'prefCohNoisePC');
end
if isnan(nTrials) && hasSessionHeader
  nTrials = localNumericField(sessionHeader, 'numberOfTrials');
end
if isnan(nNoiseTrials) && hasSessionHeader
  nNoiseTrials = localNumericField(sessionHeader, 'nNoiseTrials');
end
if isempty(parentFileName) && hasSessionHeader
  parentFileName = localCharField(sessionHeader, 'fileName');
end

excludeReasons = {};
fileInfoRow = table( ...
  {filePath}, {folder}, {fileNameExt}, ...
  hasSessionHeader, hasSessionProbeHeader, ...
  probeDirDeg, {probeTag}, parentNProbeDirections, {parentProbeDirectionsDeg}, ...
  parentIsSingleProbe, parentIsInterleavedProbe, ...
  prefCohNoisePC, probeCohNoisePC, nTrials, nNoiseTrials, {parentFileName}, ...
  false, {excludeReasons}, ...
  'VariableNames', {'filePath','folder','fileName', ...
                    'hasSessionHeader','hasSessionProbeHeader', ...
                    'probeDirDeg','probeTag','parentNProbeDirections','parentProbeDirectionsDeg', ...
                    'parentIsSingleProbe','parentIsInterleavedProbe', ...
                    'prefCohNoisePC','probeCohNoisePC','nTrials','nNoiseTrials','parentFileName', ...
                    'isExcluded','excludeReasons'});
end

% -------------------------------------------------------------------------
function v = localNumericField(S, fieldName)
v = NaN;
if ~isstruct(S) || ~isfield(S, fieldName)
  return;
end
x = localExtractValue(S.(fieldName));
if isnumeric(x) || islogical(x)
  if ~isempty(x)
    v = double(x(1));
  end
end
end

% -------------------------------------------------------------------------
function v = localNumericVectorField(S, fieldName)
v = [];
if ~isstruct(S) || ~isfield(S, fieldName)
  return;
end
x = localExtractValue(S.(fieldName));
if isnumeric(x) || islogical(x)
  v = double(x(:)');
end
end

% -------------------------------------------------------------------------
function s = localCharField(S, fieldName)
s = '';
if ~isstruct(S) || ~isfield(S, fieldName)
  return;
end
x = localExtractValue(S.(fieldName));
if isstring(x)
  if ~isempty(x)
    s = char(x(1));
  end
elseif ischar(x)
  s = x;
elseif isnumeric(x)
  s = num2str(x(1));
end
end

% -------------------------------------------------------------------------
function v = localExtractValue(x)
if isstruct(x) && isfield(x, 'data')
  v = x.data;
else
  v = x;
end
end
