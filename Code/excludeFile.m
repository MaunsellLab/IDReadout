function exclude = excludeFile(header)
% excludeFile  Return true if this session should be excluded from summary analyses.
%
% Exclusion rules:
%   1. Explicitly excluded source files
%   2. prefNoiseCohPC must equal 10
%
% Input:
%   header   session header struct, expected to contain:
%              header.fileName
%              header.prefNoiseCohPC.data

excludedFiles = { ...
  'IDReadout_Meetz_20260114.dat', ...
  'IDReadout_Meetz_20260114_2.dat', ...
  'IDReadout_Meetz_20260114_3.dat'};

exclude = false;

% Rule 1: explicit file exclusion
if isfield(header, 'fileName')
  thisFile = localExtractHeaderValue(header.fileName);
  if isstring(thisFile)
    thisFile = char(thisFile);
  end
  if ismember(thisFile, excludedFiles)
    exclude = true;
    return
  end
end

% Rule 2: require pref noise amplitude = 10
if isfield(header, 'prefNoiseCohPC')
  prefNoise = localExtractHeaderValue(header.prefNoiseCohPC);
  exclude = ~(isfinite(prefNoise) && abs(prefNoise - 10) < 1e-9);
else
  % If missing, safest behavior is to exclude
  exclude = true;
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