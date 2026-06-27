function sessionHeader = makeSessionHeader(header, trialMeta)
% makeSessionHeader  Build stable session-level metadata for MT readout analyses.
%
% sessionHeader is the analysis-facing metadata inherited from the converted
% Lablib header. It intentionally contains only session-level quantities. Do
% not add trial-level fields such as cohNoise.

if nargin < 2 || isempty(trialMeta)
  trialMeta = struct();
end

sessionHeader = struct();

% Parent-session probe-direction provenance. These describe the recording
% session as acquired, before trial correction, exclusion, or probe splitting.
sessionHeader.nProbeDirections = getMetaField(trialMeta, 'nProbeDirections', NaN);
sessionHeader.probeDirectionsDeg = getMetaField(trialMeta, 'probeDirectionsDeg', []);
sessionHeader.probeTags = getMetaField(trialMeta, 'probeTags', {});
sessionHeader.nNoiseTrials = getMetaField(trialMeta, 'nNoiseTrials', NaN);

% Stable fields used by downstream analyses or provenance displays.
copyFields = { ...
  'date', ...
  'subject', ...
  'taskName', ...
  'numberOfTrials', ...
  'frameRateHz', ...
  'preStepMS', ...
  'stepMS', ...
  'cohNoiseFrameMS', ...
  'prefDirDeg', ...
  'prefCohNoisePC' ...
  };

for k = 1:numel(copyFields)
  f = copyFields{k};
  if isfield(header, f)
    sessionHeader.(f) = localDataValue(header.(f));
  end
end

[~, sessionHeader.fileName] = fileparts(header.fileName);

% Modern files often store the coherence-noise amplitudes in blockStatus.
% Bring the preferred-stream amplitude into sessionHeader using the standard
% field name, because it is a session-level constant in this analysis.
if ~isfield(sessionHeader, 'prefCohNoisePC')
  [v, ok] = localGetBlockStatusValue(header, 'prefCohNoisePC');
  if ok
    sessionHeader.prefCohNoisePC = v;
  end
end
if ~isfield(sessionHeader, 'cohNoiseFrameMS')
  [v, ok] = localGetBlockStatusValue(header, 'cohNoiseFrameMS');
  if ok
    sessionHeader.cohNoiseFrameMS = v;
  end
end
if ~isfield(sessionHeader, 'prefCohNoisePC') && isfield(header, 'prefNoiseCohPC')
  sessionHeader.prefCohNoisePC = localDataValue(header.prefNoiseCohPC);
end

sessionHeader.createdBy = mfilename;
sessionHeader.createdDate = datetime('now');
sessionHeader.metadataContract = 'session-level analysis metadata; no trial-level fields';
end

% -------------------------------------------------------------------------
function [v, ok] = localGetBlockStatusValue(header, fieldName)
ok = false;
v = [];
if isfield(header, 'blockStatus') && isfield(header.blockStatus, 'data') && ...
    isfield(header.blockStatus.data, fieldName)
  x = header.blockStatus.data.(fieldName);
  v = localExtractValue(x);
  ok = true;
end
end

% -------------------------------------------------------------------------
function v = localDataValue(x)

v = x;
while isstruct(v) && isfield(v, 'data')
  v = v.data;
end
if isnumeric(v) && ~isscalar(v)
  v = v(1);
end
end

% -------------------------------------------------------------------------
function v = localExtractValue(x)

v = x;

while isstruct(v) && isfield(v, 'data')
  v = v.data;
end

if isnumeric(v) && ~isscalar(v)
  v = v(1);
end
end

% -------------------------------------------------------------------------

function v = getMetaField(S, fieldName, defaultValue)

if isstruct(S) && isfield(S, fieldName)
  v = S.(fieldName);
else
  v = defaultValue;
end

if isnumeric(v)
  v = v(:)';
end
end