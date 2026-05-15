function [exclude, reason] = excludeFile(header)
% excludeFile  Return true if this session should be excluded from summary analyses.
%
% Exclusion rules:
%   1. Explicitly excluded source files
%   2. prefCohNoisePC must equal 10
%   3. probeCohNoisePC must be present and positive
%
% This function is used during preprocessing/averaging. It should not exclude
% sessions merely because a probe offset or probe amplitude is not eligible
% for a later paired-probe readout fit. Those rules belong in
% updateAcrossOffsetSummaries.
%
% Input:
%   header   session header struct
%
% Output:
%   exclude  logical
%   reason   char string explaining exclusion, empty if included

excludedFiles = { ...
  'IDReadout_Meetz_20260114.dat', ...
  'IDReadout_Meetz_20260114_2.dat', ...
  'IDReadout_Meetz_20260114_3.dat'};

targetPrefNoisePC  = 10;
tol = 1e-6;
exclude = false;
reason = '';

% Rule 1: explicit file exclusion
if isfield(header, 'fileName')
  thisFile = localExtractHeaderValue(header.fileName);
  if isstring(thisFile)
    thisFile = char(thisFile);
  end
  if ismember(thisFile, excludedFiles)
    exclude = true;
    reason = 'explicitly excluded file';
    return
  end
end

% Rule 2: require pref noise amplitude = 10
[prefNoise, ok] = localGetNoisePC(header, 'pref');
if ~ok || ~isfinite(prefNoise)
  exclude = true;
  reason = 'missing pref noise amplitude';
  return
end
if abs(prefNoise - targetPrefNoisePC) > tol
  exclude = true;
  reason = sprintf('prefCohNoisePC %.6g != %.6g', prefNoise, targetPrefNoisePC);
  return
end

% Rule 3: probe noise amplitude must be present and positive.
% Do not enforce the paired-probe 10/sqrt(2) convention here: 180 deg
% single-stream sessions legitimately use probeCohNoisePC = 10 and should
% still be processed through kernel generation and averaging.  Paired-probe
% amplitude rules belong in updateAcrossOffsetSummaries.
[probeNoise, ok] = localGetNoisePC(header, 'probe');
if ~ok || ~isfinite(probeNoise)
  exclude = true;
  reason = 'missing probe noise amplitude';
  return
end
if probeNoise <= 0
  exclude = true;
  reason = sprintf('probeCohNoisePC %.6g is not positive', probeNoise);
  return
end

end


function [noisePC, ok] = localGetNoisePC(header, streamName)
% localGetNoisePC  Return pref/probe coherence-noise amplitude.
%
% Prefer the blockStatus location used by computeSessionKernels, but retain
% backward compatibility with older top-level header fields.

ok = false;
noisePC = NaN;

switch lower(streamName)
  case 'pref'
    blockField = 'prefCohNoisePC';
    legacyField = 'prefNoiseCohPC';
  case 'probe'
    blockField = 'probeCohNoisePC';
    legacyField = 'probeNoiseCohPC';
  otherwise
    error('excludeFile:UnknownStream', 'Unknown streamName: %s', streamName);
end

if isfield(header, 'blockStatus') && ...
    isfield(header.blockStatus, 'data') && ...
    isfield(header.blockStatus.data, blockField)
  noisePC = header.blockStatus.data.(blockField);
  ok = true;
  return
end

if isfield(header, blockField)
  noisePC = localExtractHeaderValue(header.(blockField));
  ok = true;
  return
end

if isfield(header, legacyField)
  noisePC = localExtractHeaderValue(header.(legacyField));
  ok = true;
end
end


function val = localExtractHeaderValue(x)
% localExtractHeaderValue  Pull scalar/string value from a header field.

if isstruct(x) && isfield(x, 'data')
  val = x.data;
else
  val = x;
end

if isnumeric(val) && ~isscalar(val)
  val = val(1);
end
end