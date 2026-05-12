function [exclude, reason] = excludeFile(header)
% excludeFile  Return true if this session should be excluded from summary analyses.
%
% Exclusion rules:
%   1. Explicitly excluded source files
%   2. prefCohNoisePC must equal 10
%   3. probeCohNoisePC must equal 10/sqrt(2)
%   4. probeDirDeg must be a paired-probe offset: 0 < probeDirDeg < 180
%
% Rule 4 intentionally excludes the legacy single-stream 180 deg condition.
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
targetProbeNoisePCs = [7, 10 / sqrt(2)];  % historical stored value and exact intended value
tol = 1e-4;tol = 1e-6;

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

%   3. probeCohNoisePC must equal the paired-probe convention:
%      either historical stored value 7 or exact intended value 10/sqrt(2)
[probeNoise, ok] = localGetNoisePC(header, 'probe');
if ~ok || ~isfinite(probeNoise)
  exclude = true;
  reason = 'missing probe noise amplitude';
  return
end
if all(abs(probeNoise - targetProbeNoisePCs) > tol)
  exclude = true;
  reason = sprintf('probeCohNoisePC %.6g not in accepted values [%s]', ...
    probeNoise, sprintf('%.6g ', targetProbeNoisePCs));
  return
end

% Rule 4: require paired-probe offsets only.
% Current paired convention is 0 < probeDirDeg < 180. The legacy 180 deg
% condition used a single stream and is excluded from the cleaned analysis.
if ~isfield(header, 'probeDirDeg')
  exclude = true;
  reason = 'missing probeDirDeg';
  return
end

probeDirDeg = abs(double(localExtractHeaderValue(header.probeDirDeg)));
if ~(isfinite(probeDirDeg) && probeDirDeg > 0 && probeDirDeg < 180)
  exclude = true;
  reason = sprintf('probeDirDeg %.6g is not a paired-probe offset', probeDirDeg);
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