function [files, fileInfo] = selectAnalysisFiles(dataFolders, varargin)
% selectAnalysisFiles  Select MT readout analysis files by explicit metadata.
%
%   [files, fileInfo] = selectAnalysisFiles(dataFolders, ...)
%
% File structure is represented by hasSessionHeader / hasSessionProbeHeader.
% Parent acquisition context is represented by parentNProbeDirections and
% parentProbeDirectionsDeg. Filename conventions are used only for hard
% exclusions before loading metadata.
%
% Supported selection arguments:
%   Animal                            char array with the name of selection
%   Bin179With180                     count 179° and 180° as a single offset
%   FilePattern                       default '*.mat'
%   FileSelectionArgs                 cell array of additional args
%   HardExcludeFileNames              filenames excluded before loading
%   MaxParentNProbeDirections         maximum parent n-probe-directions
%   MinParentNProbeDirections         minimum parent n-probe-directions
%   MultipleProbeDirections           >1 parent probe direction
%   ParentNProbeDirections            exact parent n-probe-direction match
%   ProbeDirDeg                       require this derived probe direction
%   ReportExcluded                    print excluded files and reasons
%   SingleProbeDirection              maximum 1 parent probe direction

if isstring(dataFolders) || ischar(dataFolders)
  dataFolders = {char(dataFolders)};
end
P = makeParser();
parse(P, dataFolders, varargin{:});
R0 = P.Results;
% check for nested file selection arguments and include them if they exist
if ~isempty(R0.FileSelectionArgs)
    topArgs = removeParameterPair(varargin, 'FileSelectionArgs');
    nestedArgs = R0.FileSelectionArgs;
    P = makeParser();
    parse(P, dataFolders, topArgs{:}, nestedArgs{:});
    R = P.Results;
else
    R = R0;
end

% collect selected files from all folders
fileInfo = [];
excludedRows = [];
for dIndex = 1:numel(R.dataFolders)
  [dirFileInfo, dirExcludedRows] = doOneDataFolder(R.dataFolders{dIndex}, R);
  fileInfo = [fileInfo; dirFileInfo]; %#ok<AGROW>
  excludedRows = [excludedRows; dirExcludedRows]; %#ok<AGROW>
  % If binning 179 with 180, check whether either offset is selected
  if R.Bin179With180
    [path, name] = fileparts(R.dataFolders{dIndex});
    otherDataFolder = [];
    if contains(path, 'Probe179')
      otherDataFolder = [path(1:end-numel('Probe179')) 'Probe180/' name];
      otherProbeDir = 180;
    elseif contains(path, 'Probe180')
      otherDataFolder = [path(1:end-numel('Probe180')) 'Probe179/' name];
      otherProbeDir = 179;
    end
    % If 179 or 180 has been selected, concatenate entries for the other;
    if ~isempty(otherDataFolder) && exist(otherDataFolder, 'dir')
      modArgs = removeParameterPair(varargin, 'Bin179With180');
      modArgs = removeParameterPair(modArgs, 'ProbeDirDeg');
      if exist('nestedArgs')
        modNestArgs = removeParameterPair(nestedArgs, 'Bin179With180');
        modNestArgs = removeParameterPair(modNestArgs, 'ProbeDirDeg');
        modArgs = [modArgs, modNestArgs]; %#ok<AGROW>
      end
      [~, otherFileInfo] = selectAnalysisFiles(otherDataFolder, ...
                                    modArgs{:}, 'Bin179With180', false, 'ProbeDirDeg', otherProbeDir);
      fileInfo = [fileInfo; otherFileInfo]; %#ok<AGROW>
    end
  end
end

if isempty(fileInfo)
  files = {};
else
  files = fileInfo.filePath;
end
if R.ReportExcluded
  reportExcludedFiles(excludedRows);
end
end

%%-----------------------------------------------------------------
function [fileInfo, excludedRows] = doOneDataFolder(dataFolder, R)
% doOneDataFolder  Process the contents of one directory

if ~exist(dataFolder, 'dir')
    error('selectAnalysisFiles:MissingFolder', 'Data folder not found: %s', dataFolder);
end
D = dir(fullfile(dataFolder, char(R.FilePattern)));
D = D(~[D.isdir]);
hardExcludeNames = cellstr(R.HardExcludeFileNames);
selectedRows = {};
excludedRows = {};

for k = 1:numel(D)
    fileName = D(k).name;
    if endsWith(fileName, '_fileInfo.mat')
        continue;
    end
    filePath = fullfile(D(k).folder, fileName);

    % Hard filename exclusions happen before loading metadata.
    % Match the filename stem exactly, not just the beginning of the name.
    [~, fileStem, ~] = fileparts(fileName);
    if any(strcmpi(fileStem, hardExcludeNames))
      excludedRows{end+1} = makeExcludedReportRow(filePath, {'hard filename exclusion'}); %#ok<AGROW>
      continue;
    end

    % If this a .mat file, we will check for exclusions reasons beyond the
    % hard exclusions. Otherwise, only the hard exclusions will apply
    if ~strcmp(R.FilePattern, '*.mat')
      fileInfoRow = table({filePath}, {dataFolder}, {fileStem}, 'VariableNames', {'filePath','folder','fileName'});
      selectedRows{end+1} = fileInfoRow; %#ok<AGROW>
      continue;
    end

    row = fileInfoFromAnalysisFile(filePath);
    row = stripExclusionColumns(row);

    excludeReasons = {};

    if ~isempty(R.Animal)
      if isempty(row.animal)
          error('selectAnalysisFiles:MissingAnimal', ...
              'Selection by ProbeDirDeg requires sessionProbeHeader.probeDirDeg in %s.', filePath);
      end
      if ~strcmpi(R.Animal, 'All') & ~strcmpi(row.animal, R.Animal)
          excludeReasons{end+1} = sprintf('animal %s != %s', char(row.animal), R.Animal); %#ok<AGROW>
      end
    end

    if ~isempty(R.ProbeDirDeg)
        if ~isfinite(row.probeDirDeg)
            error('selectAnalysisFiles:MissingProbeDirDeg', ...
                'Selection by ProbeDirDeg requires sessionProbeHeader.probeDirDeg in %s.', filePath);
        end

        if abs(row.probeDirDeg - R.ProbeDirDeg) > 1e-9
            excludeReasons{end+1} = sprintf('probeDirDeg %.6g != %.6g', ...
                row.probeDirDeg, R.ProbeDirDeg); %#ok<AGROW>
        end
    end

    if ~isempty(R.ParentNProbeDirections)
        requireParentNProbeDirections(row, filePath, 'ParentNProbeDirections');

        if row.parentNProbeDirections ~= R.ParentNProbeDirections
            excludeReasons{end+1} = sprintf('parentNProbeDirections %.6g != %.6g', ...
                row.parentNProbeDirections, R.ParentNProbeDirections); %#ok<AGROW>
        end
    end

    if ~isempty(R.MinParentNProbeDirections)
        requireParentNProbeDirections(row, filePath, 'MinParentNProbeDirections');

        if row.parentNProbeDirections < R.MinParentNProbeDirections
            excludeReasons{end+1} = sprintf('parentNProbeDirections %.6g < %.6g', ...
                row.parentNProbeDirections, R.MinParentNProbeDirections); %#ok<AGROW>
        end
    end

    if ~isempty(R.MaxParentNProbeDirections)
        requireParentNProbeDirections(row, filePath, 'MaxParentNProbeDirections');

        if row.parentNProbeDirections > R.MaxParentNProbeDirections
            excludeReasons{end+1} = sprintf('parentNProbeDirections %.6g > %.6g', ...
                row.parentNProbeDirections, R.MaxParentNProbeDirections); %#ok<AGROW>
        end
    end

    if ~isempty(R.SingleProbeDirection) && R.SingleProbeDirection
      if (~isempty(R.MultipleProbeDirections) && R.MultipleProbeDirections)
        error('selectAnalysisFiles: cannot specify both singleProbeDireciton and multipleProbeDirections');
      end
      requireParentNProbeDirections(row, filePath, 'SingleProbeDirection');
      if row.parentNProbeDirections > 1
        excludeReasons{end+1} = sprintf('parentNProbeDirections %.6g > 1', row.parentNProbeDirections); %#ok<AGROW>
      end
    end

    if ~isempty(R.MultipleProbeDirections) && R.MultipleProbeDirections
      if (~isempty(R.SingleProbeDirection) && R.SingleProbeDirection)
        error('selectAnalysisFiles: cannot specify both singleProbeDireciton and multipleProbeDirections');
      end
      requireParentNProbeDirections(row, filePath, 'MultipleProbeDirections');
      if row.parentNProbeDirections < 2
        excludeReasons{end+1} = sprintf('parentNProbeDirections %.6g == 1', row.parentNProbeDirections); %#ok<AGROW>
      end
    end
    if isempty(excludeReasons)
        selectedRows{end+1} = row; %#ok<AGROW>
    else
        excludedRows{end+1} = makeExcludedReportRow(filePath, excludeReasons); %#ok<AGROW>
    end
end

if isempty(selectedRows)
    fileInfo = emptyFileInfoTable();
else
    fileInfo = vertcat(selectedRows{:});
end
end

%% -------------------------------------------------------------------------
function P = makeParser()
P = inputParser;
P.FunctionName = mfilename;

% addRequired(P,  'dataFolders', @(x) iscell(x) && ~isempty(x) && all(cellfun(@ischar, x(:))));
addRequired(P,  'dataFolders', @(x) iscell(x) && ~isempty(x));
addParameter(P, 'Animal', 'All', @(x) ischar(x) || isstring(x));
addParameter(P, 'Bin179With180', false, @(x) islogical(x) && isscalar(x));
addParameter(P, 'FilePattern', '*.mat', @(x) ischar(x) || isstring(x));
addParameter(P, 'FileSelectionArgs', {}, @(x) iscell(x));
addParameter(P, 'MaxParentNProbeDirections', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'MinParentNProbeDirections', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'MultipleProbeDirections', false, @(x) islogical(x) && isscalar(x));
addParameter(P, 'ParentNProbeDirections', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'ProbeDirDeg', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'ReportExcluded', false, @(x) islogical(x) && isscalar(x));
addParameter(P, 'SingleProbeDirection', false, @(x) islogical(x) && isscalar(x));

addParameter(P, 'HardExcludeFileNames', { ...
    'IDReadout_Meetz_20260114', ...
    'IDReadout_Meetz_20260114_2', ...
    'IDReadout_Meetz_20260114_3', ...
    'IDReadout_Meetz_20260303', ...
    'IDReadout_Meetz_20260304', ...
    'IDReadout_Meetz_20260305', ...
    'IDReadout_Meetz_20260306', ...
    'IDReadout_Meetz_20260309', ...
    'IDReadout_Meetz_20260310', ...
    'IDReadout_Meetz_20260311', ...
    'IDReadout_Meetz_20260312', ...
    }, @(x) iscellstr(x) || isstring(x));
end

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
function row = stripExclusionColumns(row)
% Protect the selector contract: returned fileInfo must describe selected
% files only, and must not carry stale exclusion bookkeeping columns.

varsToDrop = intersect({'isExcluded', 'excludeReasons'}, row.Properties.VariableNames);

if ~isempty(varsToDrop)
    row(:, varsToDrop) = [];
end
end

% -------------------------------------------------------------------------
function row = makeExcludedReportRow(filePath, reasons)
[~, fileName, ext] = fileparts(filePath);

row = table( ...
    {[fileName ext]}, ...
    {filePath}, ...
    {reasons(:)'}, ...
    'VariableNames', {'fileName', 'filePath', 'excludeReasons'});
end

% -------------------------------------------------------------------------
function reportExcludedFiles(excludedRows)
if isempty(excludedRows)
    fprintf('selectAnalysisFiles: no files were excluded.\n');
    return;
end

excludedInfo = vertcat(excludedRows{:});

fprintf('\nselectAnalysisFiles excluded %d file(s):\n', height(excludedInfo));

for ii = 1:height(excludedInfo)
    fprintf('  %s\n', excludedInfo.fileName{ii});

    reasons = excludedInfo.excludeReasons{ii};
    for rr = 1:numel(reasons)
        fprintf('    - %s\n', reasons{rr});
    end
end

fprintf('\n');
end

%% -------------------------------------------------------------------------
function T = emptyFileInfoTable()
T = table( ...
    cell(0,1), cell(0,1), cell(0,1), ...
    false(0,1), false(0,1), ...
    zeros(0,1), cell(0,1), zeros(0,1), cell(0,1), ...
    false(0,1), false(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), cell(0,1), ...
    'VariableNames', {'filePath','folder','fileName', ...
    'hasSessionHeader','hasSessionProbeHeader', ...
    'probeDirDeg','probeTag','parentNProbeDirections','parentProbeDirectionsDeg', ...
    'parentIsSingleProbe','parentIsInterleavedProbe', ...
    'prefCohNoisePC','probeCohNoisePC','nTrials','nNoiseTrials','parentFileName'});
end

%% -------------------------------------------------------------------------
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
{animal}, probeDirDeg, {probeTag}, parentNProbeDirections, {parentProbeDirectionsDeg}, ...
parentIsSingleProbe, parentIsInterleavedProbe, ...
prefCohNoisePC, probeCohNoisePC, nTrials, nNoiseTrials, {parentFileName}, ...
false, {excludeReasons}, ...
'VariableNames', {'filePath','folder','fileName', 'animal', 'probeDirDeg','probeTag','parentNProbeDirections', ...
        'parentProbeDirectionsDeg', 'parentIsSingleProbe','parentIsInterleavedProbe', 'prefCohNoisePC', ...
        'probeCohNoisePC','nTrials','nNoiseTrials','parentFileName', 'isExcluded','excludeReasons'});
end
