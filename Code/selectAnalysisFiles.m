function [files, fileInfo] = selectAnalysisFiles(dataFolder, varargin)
% selectAnalysisFiles  Select MT readout analysis files by explicit metadata.
%
%   [files, fileInfo] = selectAnalysisFiles(dataFolder, ...)
%
% File structure is represented by hasSessionHeader / hasSessionProbeHeader.
% Parent acquisition context is represented by parentNProbeDirections and
% parentProbeDirectionsDeg. Filename conventions are used only for hard
% exclusions before loading metadata.
%
% Supported selection arguments:
%   FilePattern                       default '*.mat'
%   FileSelectionArgs                 cell array of additional args
%   RequireSessionHeader              require sessionHeader variable
%   RequireSessionProbeHeader         require sessionProbeHeader variable
%   ProbeDirDeg                       require this derived probe direction
%   ParentNProbeDirections            exact parent n-probe-direction match
%   MinParentNProbeDirections         minimum parent n-probe-directions
%   MaxParentNProbeDirections         maximum parent n-probe-directions
%   ApplyExperimentalValidityChecks   call applyExperimentalValidityChecks
%   IncludeExcluded                   return excluded rows too
%   HardExcludeFileNames              filenames excluded before loading

% Parse the input arguments
% P = makeParser();
P = inputParser;
P.FunctionName = mfilename;
addRequired(P, 'dataFolder', @(x) ischar(x) || isstring(x));
addParameter(P, 'FilePattern', '*.mat', @(x) ischar(x) || isstring(x));
addParameter(P, 'FileSelectionArgs', {}, @(x) iscell(x));
addParameter(P, 'RequireSessionHeader', false, @(x) islogical(x) && isscalar(x));
addParameter(P, 'RequireSessionProbeHeader', false, @(x) islogical(x) && isscalar(x));
addParameter(P, 'ProbeDirDeg', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'ParentNProbeDirections', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'MinParentNProbeDirections', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'MaxParentNProbeDirections', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'ApplyExperimentalValidityChecks', false, @(x) islogical(x) && isscalar(x));
addParameter(P, 'IncludeExcluded', false, @(x) islogical(x) && isscalar(x));
addParameter(P, 'HardExcludeFileNames', ...
  { ...
  'IDReadout_Meetz_20260114.dat', ...
  'IDReadout_Meetz_20260114_2.dat', ...
  'IDReadout_Meetz_20260114_3.dat', ...
  'IDReadout_Meetz_20260304.dat', ...
  'IDReadout_Meetz_20260305.dat', ...
  'IDReadout_Meetz_20260306.dat', ...
  'IDReadout_Meetz_20260309.dat', ...
  'IDReadout_Meetz_20260310.dat', ...
  'IDReadout_Meetz_20260311.dat', ...
  'IDReadout_Meetz_20260312.dat', ...
  }, @(x) iscellstr(x) || isstring(x));
parse(P, dataFolder, varargin{:});
R0 = P.Results;
if ~isempty(R0.FileSelectionArgs)
  topArgs = removeParameterPair(varargin, 'FileSelectionArgs');
  nestedArgs = R0.FileSelectionArgs;
  P = makeParser();
  parse(P, dataFolder, topArgs{:}, nestedArgs{:});
  R = P.Results;
else
  R = R0;
end
% R = parseSelectionArgs(dataFolder, varargin{:});

dataFolder = char(R.dataFolder);
if ~exist(dataFolder, 'dir')
  error('selectAnalysisFiles:MissingFolder', 'Data folder not found: %s', dataFolder);
end

D = dir(fullfile(dataFolder, char(R.FilePattern)));
D = D(~[D.isdir]);

hardExcludeNames = cellstr(R.HardExcludeFileNames);
rows = {};

for k = 1:numel(D)
  fileName = D(k).name;
  if endsWith(fileName, '_fileInfo.mat')
    continue;
  end

  filePath = fullfile(D(k).folder, fileName);

  % Hard filename exclusions happen before loading metadata.
  if any(strcmp(fileName, hardExcludeNames))
    row = makeHardExcludedRow(filePath, 'hard filename exclusion');
    rows{end+1} = row; %#ok<AGROW>
    continue;
  end

  row = fileInfoFromAnalysisFile(filePath);

  if R.RequireSessionHeader && ~row.hasSessionHeader
    error('selectAnalysisFiles:MissingSessionHeader', ...
      'Selection requires sessionHeader, but %s does not contain it.', filePath);
  end

  if R.RequireSessionProbeHeader && ~row.hasSessionProbeHeader
    error('selectAnalysisFiles:MissingSessionProbeHeader', ...
      'Selection requires sessionProbeHeader, but %s does not contain it.', filePath);
  end

  if ~isempty(R.ProbeDirDeg)
    if ~row.hasSessionProbeHeader || ~isfinite(row.probeDirDeg)
      error('selectAnalysisFiles:MissingProbeDirDeg', ...
        'Selection by ProbeDirDeg requires sessionProbeHeader.probeDirDeg in %s.', filePath);
    end
    if abs(row.probeDirDeg - R.ProbeDirDeg) > 1e-9
      row = addExcludeReason(row, sprintf('probeDirDeg %.6g != %.6g', row.probeDirDeg, R.ProbeDirDeg));
    end
  end

  if ~isempty(R.ParentNProbeDirections)
    requireParentNProbeDirections(row, filePath, 'ParentNProbeDirections');
    if row.parentNProbeDirections ~= R.ParentNProbeDirections
      row = addExcludeReason(row, sprintf('parentNProbeDirections %.6g != %.6g', ...
        row.parentNProbeDirections, R.ParentNProbeDirections));
    end
  end

  if ~isempty(R.MinParentNProbeDirections)
    requireParentNProbeDirections(row, filePath, 'MinParentNProbeDirections');
    if row.parentNProbeDirections < R.MinParentNProbeDirections
      row = addExcludeReason(row, sprintf('parentNProbeDirections %.6g < %.6g', ...
        row.parentNProbeDirections, R.MinParentNProbeDirections));
    end
  end

  if ~isempty(R.MaxParentNProbeDirections)
    requireParentNProbeDirections(row, filePath, 'MaxParentNProbeDirections');
    if row.parentNProbeDirections > R.MaxParentNProbeDirections
      row = addExcludeReason(row, sprintf('parentNProbeDirections %.6g > %.6g', ...
        row.parentNProbeDirections, R.MaxParentNProbeDirections));
    end
  end

  if R.ApplyExperimentalValidityChecks
    S = load(filePath, 'sessionProbeHeader', 'sessionHeader');
    if isfield(S, 'sessionProbeHeader')
      metadata = S.sessionProbeHeader;
    elseif isfield(S, 'sessionHeader')
      metadata = S.sessionHeader;
    else
      error('selectAnalysisFiles:MissingValidityMetadata', ...
        'Experimental validity checks require sessionProbeHeader or sessionHeader in %s.', filePath);
    end
    [tfExclude, reasons] = applyExperimentalValidityChecks(metadata);
    if tfExclude
      row = addExcludeReasons(row, reasons);
    end
  end

  rows{end+1} = row; %#ok<AGROW>
end

if isempty(rows)
  fileInfo = emptyFileInfoTable();
else
  fileInfo = vertcat(rows{:});
end

if ~R.IncludeExcluded && ~isempty(fileInfo)
  fileInfo = fileInfo(~fileInfo.isExcluded, :);
end

if isempty(fileInfo)
  files = {};
else
  files = fileInfo.filePath;
end
end

% % -------------------------------------------------------------------------
% function R = parseSelectionArgs(dataFolder, varargin)
% % Parse once, then reparse with nested FileSelectionArgs appended. This lets
% % callers pass selection criteria through without needing to flatten them.
% P = makeParser();
% parse(P, dataFolder, varargin{:});
% R0 = P.Results;
% 
% if ~isempty(R0.FileSelectionArgs)
%   topArgs = removeParameterPair(varargin, 'FileSelectionArgs');
%   nestedArgs = R0.FileSelectionArgs;
%   P = makeParser();
%   parse(P, dataFolder, topArgs{:}, nestedArgs{:});
%   R = P.Results;
% else
%   R = R0;
% end
% end

% % -------------------------------------------------------------------------
% function P = makeParser()
% P = inputParser;
% P.FunctionName = mfilename;
% addRequired(P, 'dataFolder', @(x) ischar(x) || isstring(x));
% addParameter(P, 'FilePattern', '*.mat', @(x) ischar(x) || isstring(x));
% addParameter(P, 'FileSelectionArgs', {}, @(x) iscell(x));
% addParameter(P, 'RequireSessionHeader', false, @(x) islogical(x) && isscalar(x));
% addParameter(P, 'RequireSessionProbeHeader', false, @(x) islogical(x) && isscalar(x));
% addParameter(P, 'ProbeDirDeg', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
% addParameter(P, 'ParentNProbeDirections', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
% addParameter(P, 'MinParentNProbeDirections', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
% addParameter(P, 'MaxParentNProbeDirections', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
% addParameter(P, 'ApplyExperimentalValidityChecks', false, @(x) islogical(x) && isscalar(x));
% addParameter(P, 'IncludeExcluded', false, @(x) islogical(x) && isscalar(x));
% addParameter(P, 'HardExcludeFileNames', ...
%                       { ...
%                       'IDReadout_Meetz_20260114.dat', ...
%                       'IDReadout_Meetz_20260114_2.dat', ...
%                       'IDReadout_Meetz_20260114_3.dat', ...
%                       'IDReadout_Meetz_20260304.dat', ...
%                       'IDReadout_Meetz_20260305.dat', ...
%                       'IDReadout_Meetz_20260306.dat', ...
%                       'IDReadout_Meetz_20260309.dat', ...
%                       'IDReadout_Meetz_20260310.dat', ...
%                       'IDReadout_Meetz_20260311.dat', ...
%                       'IDReadout_Meetz_20260312.dat', ...
%                       }, @(x) iscellstr(x) || isstring(x));
% end

% -------------------------------------------------------------------------
function argsOut = removeParameterPair(argsIn, paramName)
argsOut = {};
k = 1;
while k <= numel(argsIn)
  if (ischar(argsIn{k}) || isstring(argsIn{k})) && strcmpi(char(argsIn{k}), paramName)
    k = k + 2;
  else
    argsOut{end+1} = argsIn{k}; %#ok<AGROW>
    k = k + 1;
  end
end
end

% -------------------------------------------------------------------------
function requireParentNProbeDirections(row, filePath, criterionName)
if ~isfinite(row.parentNProbeDirections)
  error('selectAnalysisFiles:MissingParentNProbeDirections', ...
    ['Selection by %s requires parentNProbeDirections metadata in %s. ' ...
     'For probe-session files this should come from sessionProbeHeader.parentNProbeDirections; ' ...
     'for parent/session files this should come from sessionHeader.nProbeDirections.'], ...
     criterionName, filePath);
end
end

% -------------------------------------------------------------------------
function row = makeHardExcludedRow(filePath, reason)
[folder, fileName, ext] = fileparts(filePath);
row = table( ...
  {filePath}, {folder}, {[fileName ext]}, ...
  false, false, ...
  NaN, {''}, NaN, {[]}, ...
  false, false, ...
  NaN, NaN, NaN, NaN, {''}, ...
  true, {{reason}}, ...
  'VariableNames', {'filePath','folder','fileName', ...
                    'hasSessionHeader','hasSessionProbeHeader', ...
                    'probeDirDeg','probeTag','parentNProbeDirections','parentProbeDirectionsDeg', ...
                    'parentIsSingleProbe','parentIsInterleavedProbe', ...
                    'prefCohNoisePC','probeCohNoisePC','nTrials','nNoiseTrials','parentFileName', ...
                    'isExcluded','excludeReasons'});
end

% -------------------------------------------------------------------------
function T = emptyFileInfoTable()
T = table( ...
  cell(0,1), cell(0,1), cell(0,1), ...
  false(0,1), false(0,1), ...
  zeros(0,1), cell(0,1), zeros(0,1), cell(0,1), ...
  false(0,1), false(0,1), ...
  zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), cell(0,1), ...
  false(0,1), cell(0,1), ...
  'VariableNames', {'filePath','folder','fileName', ...
                    'hasSessionHeader','hasSessionProbeHeader', ...
                    'probeDirDeg','probeTag','parentNProbeDirections','parentProbeDirectionsDeg', ...
                    'parentIsSingleProbe','parentIsInterleavedProbe', ...
                    'prefCohNoisePC','probeCohNoisePC','nTrials','nNoiseTrials','parentFileName', ...
                    'isExcluded','excludeReasons'});
end

% -------------------------------------------------------------------------
function row = addExcludeReason(row, reason)
row = addExcludeReasons(row, {reason});
end

% -------------------------------------------------------------------------
function row = addExcludeReasons(row, reasons)
if isempty(reasons)
  return;
end
oldReasons = row.excludeReasons{1};
row.excludeReasons{1} = [oldReasons(:)' reasons(:)'];
row.isExcluded = true;
end
