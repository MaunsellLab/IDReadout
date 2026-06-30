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

[folder, fileName, ext] = fileparts(filePath);
fileNameExt = [fileName ext];
S = load(filePath, 'sessionHeader');
vars = who('-file', filePath);
if ismember('sessionProbeHeader', vars)
  T = load(filePath, 'sessionProbeHeader');
  S.sessionProbeHeader = T.sessionProbeHeader;
end
nTrials = S.sessionHeader.numberOfTrials';
nNoiseTrials = S.sessionHeader.nNoiseTrials';
animal = S.sessionHeader.animal;
parentFileName = S.sessionHeader.fileName';
parentProbeDirectionsDeg = S.sessionHeader.probeDirectionsDeg';
parentNProbeDirections = S.sessionHeader.nProbeDirections;
parentIsSingleProbe = parentNProbeDirections == 1;
parentIsInterleavedProbe = parentNProbeDirections > 1;
prefCohNoisePC = S.sessionHeader.prefCohNoisePC;
if isfield(S, 'sessionProbeHeader')
  probeDirDeg = S.sessionProbeHeader.probeDirDeg;
  probeTag = S.sessionProbeHeader.probeTag;
  probeCohNoisePC = S.sessionProbeHeader.probeCohNoisePC';
else
  probeCohNoisePC = nan;
  probeDirDeg = nan;
  probeTag = "";
end

excludeReasons = {};

fileInfoRow = table( ...
{filePath}, {folder}, {fileNameExt}, ...
animal, probeDirDeg, {probeTag}, parentNProbeDirections, {parentProbeDirectionsDeg}, ...
parentIsSingleProbe, parentIsInterleavedProbe, ...
prefCohNoisePC, probeCohNoisePC, nTrials, nNoiseTrials, {parentFileName}, ...
false, {excludeReasons}, ...
'VariableNames', {'filePath','folder','fileName', 'animal', 'probeDirDeg','probeTag','parentNProbeDirections', ...
        'parentProbeDirectionsDeg', 'parentIsSingleProbe','parentIsInterleavedProbe', 'prefCohNoisePC', ...
        'probeCohNoisePC','nTrials','nNoiseTrials','parentFileName', 'isExcluded','excludeReasons'});
end

% -------------------------------------------------------------------------
% function v = localNumericField(S, fieldName)
% v = NaN;
% if ~isstruct(S) || ~isfield(S, fieldName)
%   return;
% end
% x = localExtractValue(S.(fieldName));
% if isnumeric(x) || islogical(x)
%   if ~isempty(x)
%     v = double(x(1));
%   end
% end
% end
% 
% % -------------------------------------------------------------------------
% function v = localNumericVectorField(S, fieldName)
% v = [];
% if ~isstruct(S) || ~isfield(S, fieldName)
%   return;
% end
% x = localExtractValue(S.(fieldName));
% if isnumeric(x) || islogical(x)
%   v = double(x(:)');
% end
% end
% 
% % -------------------------------------------------------------------------
% function s = localCharField(S, fieldName)
% s = '';
% if ~isstruct(S) || ~isfield(S, fieldName)
%   return;
% end
% x = localExtractValue(S.(fieldName));
% if isstring(x)
%   if ~isempty(x)
%     s = char(x(1));
%   end
% elseif ischar(x)
%   s = x;
% elseif isnumeric(x)
%   s = num2str(x(1));
% end
% end

% -------------------------------------------------------------------------
function v = localExtractValue(x)
if isstruct(x) && isfield(x, 'data')
  v = x.data;
else
  v = x;
end
end