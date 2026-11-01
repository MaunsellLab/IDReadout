function [files, fileInfo] = selectAnalysisFiles(dataFolders, varargin)
% selectAnalysisFiles  Select MT readout analysis files by explicit metadata.
%
%   [files, fileInfo] = selectAnalysisFiles(dataFolders, ...)
%
% Selects .mat files with specified constraints.  Decisions are based
% primary on information in sessionHeader (which every .mat file must
% contain) or sometimes sessionProbeHeader (for more detailed selectivity).
%
% dataFolders can be a single data folder, or an array of folders that will
% all be scanned to return a single list
% 
% Top level exclusions are made based on a list of hard excluded file name
% and excluded probe directions, but otherwise all decisions are based on 
% header content
%
% Hard exclusions specified by a parent-session stem also apply to its
% ProbeSession derivatives (including the IDQ ProbeSession filename prefix).
%
% Typically this function is used by passing a dataFolder and an animal
% specification.
%
% Supported selection arguments:
%   Animal                            char array with the name of selection
%   Bin180Into179                     count 179° and 180° as a single offset
%   FilePattern                       default '*.mat'
%   HardExcludeFileNames              filenames excluded before loading
%   IncludeIDQ                        include IDQ Probe120 sessions when an
%                                     IDS Probe120/ProbeSessions folder is selected
%   LoadMetadata                      load headers and apply metadata filters
%                                     (default true); false applies only hard
%                                     filename/path exclusions
%   MaxParentNProbeDirections         maximum parent n-probe-directions
%   MinParentNProbeDirections         minimum parent n-probe-directions
%   MultipleProbeDirections           >1 parent probe direction
%   ParentNProbeDirections            exact parent n-probe-direction match
%   ProbeDirDeg                       require this derived probe direction
%   RandomBaseDir                     [] for either mode, or logical scalar
%   ReportExcluded                    print excluded files and reasons
%   SingleProbeDirection              maximum 1 parent probe direction
%   SourceDomain                      'All', 'IDQ', 'IDS', etc.
%   TwoAllocationMode                 [] for either IDL style; false/true to select

if isstring(dataFolders) || ischar(dataFolders)
  dataFolders = {char(dataFolders)};
end
P = makeParser();
parse(P, dataFolders, varargin{:});
R = P.Results;
validateLoadMetadataOptions(R);
R.dataFolders = includeRelatedIDQFolders(R.dataFolders, R.IncludeIDQ);

% collect selected files from all folders
fileInfo = emptyFileInfoTable();
excludedInfo = emptyExcludeTable();
for dIndex = 1:numel(R.dataFolders)

  % Special execution is implemented for Bin180Into179, which combines
  % probe directions 179° and 180°
  % If binning 179 with 180, check whether either offset is selected
  if R.Bin180Into179
    [path, name] = fileparts(R.dataFolders{dIndex});
    if contains(path, 'Probe180')   % do not process 180 when it is requested
      continue;
    elseif contains(path, 'Probe179')
      % if this is a 179° folder, pull in data for 180°
      dataFolder180 = [path(1:end-numel('Probe179')) 'Probe180/' name];
      if exist(dataFolder180, 'dir')
        % we need to modify args to get past this exclusion code on
        % recursive call
        modArgs = removeParameterPair(varargin, 'Bin180Into179');
        modArgs = removeParameterPair(modArgs, 'ProbeDirDeg');
        [~, fileInfo180] = selectAnalysisFiles(dataFolder180, ...
          modArgs{:}, 'Bin180Into179', false, 'ProbeDirDeg', 180);
        fileInfo = [fileInfo; fileInfo180]; %#ok<AGROW>
      end
    end
  end

  [dirFileInfo, dirExcludedInfo] = doOneDataFolder(R.dataFolders{dIndex}, R);
  
  fileInfo = [fileInfo; dirFileInfo]; %#ok<AGROW>
  excludedInfo = [excludedInfo; dirExcludedInfo]; %#ok<AGROW>
end

if isempty(fileInfo)
  files = {};
else
  files = fileInfo.filePath;
end
if R.ReportExcluded
  reportExcludedFiles(excludedInfo);
end
end

%%-----------------------------------------------------------------
function [fileInfo, excludedInfo] = doOneDataFolder(dataFolder, R)
% doOneDataFolder  Process the contents of one directory

if ~exist(dataFolder, 'dir')
    error('selectAnalysisFiles:MissingFolder', 'Data folder not found: %s', dataFolder);
end
D = dir(fullfile(dataFolder, char(R.FilePattern)));
D = D(~[D.isdir]);
hardExcludeNames = cellstr(R.HardExcludeFileNames);
hardExcludedStrings = cellstr(R.HardExcludeNameStrings);
selectedRows = {};
excludedRows = {};
fileInfo = emptyFileInfoTable();
excludedInfo = emptyExcludeTable();

for k = 1:numel(D)
    fileName = D(k).name;
    % if endsWith(fileName, '_fileInfo.mat')
    %     continue;
    % end
    filePath = fullfile(D(k).folder, fileName);


    % check for any path string exclusions -- used for excluded probe
    % directions and info files
    if any(contains(filePath, hardExcludedStrings))
      excludedRows{end+1} = makeExcludedReportRow(filePath, {'hard filename exclusion'}); %#ok<AGROW>
      continue;
    end

    % Hard filename exclusions happen before loading metadata. ProbeSession
    % derivatives also match the exact stem of their parent session.
    [~, fileStem, ~] = fileparts(fileName);
    if isHardExcludedFileStem(fileStem, hardExcludeNames)
      excludedRows{end+1} = makeExcludedReportRow(filePath, {'hard filename exclusion'}); %#ok<AGROW>
      continue;
    end

    % A filename-only scan is used by the incremental max-pooling cache. It
    % deliberately stops here so unchanged MAT files are never opened.
    if ~R.LoadMetadata
      fileInfoRow = fileInfoFromPathOnly(filePath,dataFolder,fileName);
      selectedRows{end+1} = fileInfoRow; %#ok<AGROW>
      continue;
    end

    % If this a .mat file, we will check for exclusions reasons beyond the
    % hard exclusions. Otherwise, only the hard exclusions will apply
    % Any pattern ending in .mat names analysis MAT files and should receive
    % metadata-based filtering (for example '*_sessionAnalysis.mat').
    if ~endsWith(char(R.FilePattern), '.mat', 'IgnoreCase', true)
      fileInfoRow = fileInfoFromPathOnly(filePath,dataFolder,fileName);
      selectedRows{end+1} = fileInfoRow; %#ok<AGROW>
      continue;
    end

    row = fileInfoFromAnalysisFile(filePath);
    excludeReasons = {};

    % IDQ-derived analyses are stored in the IDS analysis tree so that all
    % later IDS stages use one layout. Permit an IDS-only selection without
    % relying on filename conventions.
    if ~R.IncludeIDQ && isIDSDataFolder(dataFolder) && ...
        strcmpi(row.sourceDomain{1}, 'IDQ')
      excludeReasons{end+1} = 'IDQ-derived session excluded'; %#ok<AGROW>
    end

    if ~strcmpi(char(R.SourceDomain), 'All')
      if isempty(row.sourceDomain{1})
        error('selectAnalysisFiles:MissingSourceDomain', ...
          'Selection by SourceDomain requires sourceDomain metadata in %s.', filePath);
      end
      if ~strcmpi(row.sourceDomain{1}, char(R.SourceDomain))
        excludeReasons{end+1} = sprintf('sourceDomain %s != %s', ...
          row.sourceDomain{1}, char(R.SourceDomain)); %#ok<AGROW>
      end
    end

    if ~isempty(R.TwoAllocationMode)
      if ~strcmpi(row.sourceDomain{1},'IDL')
        excludeReasons{end+1} = 'TwoAllocationMode selects IDL sessions only'; %#ok<AGROW>
      elseif ~isfinite(row.twoAllocationMode)
        error('selectAnalysisFiles:MissingIDLAllocationMode', ...
          'IDL allocation header is missing in %s. Run makeIDLSessionAnalyses first.',filePath);
      elseif logical(row.twoAllocationMode)~=R.TwoAllocationMode
        excludeReasons{end+1} = 'different IDL allocation style'; %#ok<AGROW>
      end
    end

    % RandomBaseDir is an IDQ subgroup criterion in the current analysis.
    % Do not require that metadata from files already excluded by an
    % explicit SourceDomain selection (for example, native IDS files while
    % selecting randomized IDQ sessions).
    domainMatches = strcmpi(char(R.SourceDomain), 'All') || ...
      strcmpi(row.sourceDomain{1}, char(R.SourceDomain));
    if ~isempty(R.RandomBaseDir) && domainMatches
      if ~isfinite(row.randomBaseDir)
        error('selectAnalysisFiles:MissingRandomBaseDir', ...
          'Selection by RandomBaseDir requires randomBaseDir metadata in %s.', filePath);
      end
      if logical(row.randomBaseDir) ~= R.RandomBaseDir
        excludeReasons{end+1} = sprintf('randomBaseDir %d != %d', ...
          logical(row.randomBaseDir), R.RandomBaseDir); %#ok<AGROW>
      end
    end

    if ~isempty(R.Animal)
      if isempty(row.animal)
          error('selectAnalysisFiles:MissingAnimal', ...
              'Selection by Animal requires sessionHeader.animal or sessionHeader.subject in %s.', filePath);
      end
      if ~strcmpi(R.Animal, 'All') & ~strcmpi(row.animal, R.Animal)
          excludeReasons{end+1} = sprintf('animal %s != %s', char(row.animal), R.Animal); %#ok<AGROW>
      end
    end

    % if ~isempty(R.PrefCohNoisePC)
    %   if ~isfinite(row.prefCohNoisePC)
    %     error('selectAnalysisFiles:MissingPrefCohNoisePC', ...
    %       'Selection by PrefCohNoisePC requires sessionProbeHeader.PrefCohNoisePC in %s.', filePath);
    %   end
    % 
    %   if abs(row.prefCohNoisePC - R.PrefCohNoisePC) > 1e-9
    %     excludeReasons{end+1} = sprintf('prefCohNoisePC %.6g != %.6g', ...
    %       row.prefCohNoisePC, R.PrefCohNoisePC); %#ok<AGROW>
    %   end
    % end

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

if ~isempty(selectedRows)
    fileInfo = vertcat(selectedRows{:});
end
if ~isempty(excludedRows)
  excludedInfo = vertcat(excludedRows{:});
end
end

%% -------------------------------------------------------------------------
function P = makeParser()
P = inputParser;
P.FunctionName = mfilename;

addRequired(P,  'dataFolders', @(x) iscell(x) && ~isempty(x));
addParameter(P, 'Animal', 'All', @(x) ischar(x) || isstring(x));
addParameter(P, 'Bin180Into179', true, @(x) islogical(x) && isscalar(x));
addParameter(P, 'FilePattern', '*.mat', @(x) ischar(x) || isstring(x));
addParameter(P, 'IncludeIDQ', true, @(x) islogical(x) && isscalar(x));
addParameter(P, 'LoadMetadata', true, @(x) islogical(x) && isscalar(x));
addParameter(P, 'MaxParentNProbeDirections', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'MinParentNProbeDirections', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'MultipleProbeDirections', false, @(x) islogical(x) && isscalar(x));
addParameter(P, 'ParentNProbeDirections', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
% addParameter(P, 'PrefCohNoisePC', 5, @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'ProbeDirDeg', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'RandomBaseDir', [], ...
  @(x) isempty(x) || (islogical(x) && isscalar(x)));
addParameter(P, 'ReportExcluded', false, @(x) islogical(x) && isscalar(x));
addParameter(P, 'SingleProbeDirection', false, @(x) islogical(x) && isscalar(x));
addParameter(P, 'SourceDomain', 'All', @(x) ischar(x) || isstring(x));
addParameter(P, 'TwoAllocationMode', [], ...
  @(x) isempty(x) || (islogical(x) && isscalar(x)));

% addParameter(P, 'HardExcludeFileNames', { ...
%     'IDReadout_Meetz_20260114', ...
%     'IDReadout_Meetz_20260114_2', ...
%     'IDReadout_Meetz_20260114_3', ...
%     'IDReadout_Meetz_20260303', ...
%     'IDReadout_Meetz_20260304', ...
%     'IDReadout_Meetz_20260305', ...
%     'IDReadout_Meetz_20260306', ...
%     'IDReadout_Meetz_20260309', ...
%     'IDReadout_Meetz_20260310', ...
%     'IDReadout_Meetz_20260311', ...
%     'IDReadout_Meetz_20260312', ...
%     }, @(x) iscellstr(x) || isstring(x));
addParameter(P, 'HardExcludeFileNames', { ...
  'IDReadout_Meetz_20260114', ...
  'IDReadout_Meetz_20260114_2', ...
  'IDReadout_Meetz_20260114_3', ...
  }, @(x) iscellstr(x) || isstring(x));

addParameter(P, 'HardExcludeNameStrings', { ...
  '_fileInfo.mat', ...
  'Probe1/', ...
  }, @(x) iscellstr(x) || isstring(x));

end

% -------------------------------------------------------------------------
function validateLoadMetadataOptions(R)
% Filename-only mode supports hard exclusions, not header-based criteria.

if R.LoadMetadata
  return
end
hasMetadataCriterion = ~strcmpi(char(R.Animal),'All') || ...
  ~strcmpi(char(R.SourceDomain),'All') || ~isempty(R.RandomBaseDir) || ...
  ~isempty(R.TwoAllocationMode) || ...
  ~isempty(R.ProbeDirDeg) || ~isempty(R.ParentNProbeDirections) || ...
  ~isempty(R.MinParentNProbeDirections) || ...
  ~isempty(R.MaxParentNProbeDirections) || R.SingleProbeDirection || ...
  R.MultipleProbeDirections;
if hasMetadataCriterion
  error('selectAnalysisFiles:MetadataDisabled', ...
    ['LoadMetadata=false cannot be combined with Animal or header-based ' ...
     'selection criteria. Filter the filename-only result separately.']);
end
end

% -------------------------------------------------------------------------
function tf = isHardExcludedFileStem(fileStem, hardExcludeNames)
% Match either the derived filename or its parent-session filename stem.
% ProbeSessionIDR3_X_Probe120 -> IDR3_X
% IDReadout2_X_Probe25        -> IDReadout2_X

candidate = char(fileStem);
parentCandidate = regexprep(candidate, '^ProbeSession', '', 'ignorecase');
parentCandidate = regexprep(parentCandidate, ...
  '_Probe[-+]?\d+(?:\.\d+)?$', '', 'ignorecase');
tf = any(strcmpi(candidate,hardExcludeNames)) || ...
  any(strcmpi(parentCandidate,hardExcludeNames));
end

% -------------------------------------------------------------------------
function folders = includeRelatedIDQFolders(folders, includeIDQ)
% Add IDQ/Data/Probe120 for selections rooted at
% IDS/Data/Probe120/ProbeSessions. Other IDS folders and explicit IDQ
% selections are unchanged.

if ~includeIDQ
  return
end

related = {};
for k = 1:numel(folders)
  folder = char(folders{k});
  [probeFolder, leaf] = fileparts(stripTrailingFilesep(folder));
  [dataFolder, probeTag] = fileparts(probeFolder);
  [domainFolderPath, dataLeaf] = fileparts(dataFolder);
  [domainsFolder, domainLeaf] = fileparts(domainFolderPath);
  if strcmpi(leaf, 'ProbeSessions') && strcmpi(probeTag, 'Probe120') && ...
      strcmpi(dataLeaf, 'Data') && strcmpi(domainLeaf, 'IDS')
    idqFolder = fullfile(domainsFolder, 'IDQ', 'Data', 'Probe120');
    if isfolder(idqFolder)
      related{end + 1} = idqFolder; %#ok<AGROW>
    end
  end
end
folders = unique([folders(:); related(:)], 'stable');
end

% -------------------------------------------------------------------------
function folder = stripTrailingFilesep(folder)
while numel(folder) > 1 && (folder(end) == '/' || folder(end) == '\')
  folder(end) = [];
end
end

% -------------------------------------------------------------------------
function tf = isIDSDataFolder(folder)
folder = strrep(stripTrailingFilesep(char(folder)), '\', '/');
tf = contains(folder, '/IDS/Data/');
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

%% -------------------------------------------------------------------------
function T = emptyExcludeTable()

T = table(cell(0,1), cell(0,1), cell(0,1), 'VariableNames', {'fileName', 'filePath', 'excludeReasons'});
end

% -------------------------------------------------------------------------
function row = makeExcludedReportRow(filePath, reasons)

[~, fileName, ext] = fileparts(filePath);
row = table({[fileName ext]}, {filePath}, {reasons(:)'}, 'VariableNames', {'fileName', 'filePath', 'excludeReasons'});
end

% -------------------------------------------------------------------------
function reportExcludedFiles(excludedInfo)
if isempty(excludedInfo)
    fprintf('selectAnalysisFiles: no files were excluded.\n');
    return;
end
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
    cell(0,1), cell(0,1), cell(0,1), cell(0,1), cell(0,1), ...
    zeros(0,1), zeros(0,1), cell(0,1), zeros(0,1), cell(0,1), ...
    false(0,1), false(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), cell(0,1), ...
    'VariableNames', {'filePath', 'folder', 'fileName', 'animal', 'sourceDomain', ...
    'randomBaseDir','probeDirDeg','probeTag','parentNProbeDirections','parentProbeDirectionsDeg', ...
    'parentIsSingleProbe','parentIsInterleavedProbe', ...
    'prefCohNoisePC','probeCohNoisePC','nTrials','nNoiseTrials','parentFileName'});
T.twoAllocationMode = nan(0,1);
end

% -------------------------------------------------------------------------
function row = fileInfoFromPathOnly(filePath,dataFolder,fileName)
% Preserve the normal fileInfo schema without opening the selected file.

row = table( ...
    {filePath}, {dataFolder}, {fileName}, {''}, {''}, ...
    nan, nan, {''}, nan, {[]}, false, false, ...
    nan, nan, nan, nan, {''}, ...
    'VariableNames', {'filePath', 'folder', 'fileName', 'animal', 'sourceDomain', ...
    'randomBaseDir','probeDirDeg','probeTag','parentNProbeDirections','parentProbeDirectionsDeg', ...
    'parentIsSingleProbe','parentIsInterleavedProbe', ...
    'prefCohNoisePC','probeCohNoisePC','nTrials','nNoiseTrials','parentFileName'});
row.twoAllocationMode = NaN;
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
H = S.sessionHeader;

% IDR has historically supplied a rich sessionHeader. IDS initially used a
% deliberately lean parent header and stores offset-specific metadata in
% sessionProbeHeader. Accept both contracts so that the common file selector
% can perform its filename exclusions and later support Animal selection for
% either analysis without requiring task-specific directory scans.
nTrials = firstHeaderValue(H, {'numberOfTrials','nTrials'}, nan);
nNoiseTrials = firstHeaderValue(H, {'nNoiseTrials'}, nan);
animal = firstHeaderValue(H, {'animal','subject'}, '');
sourceDomain = firstHeaderValue(H, {'sourceDomain'}, '');
randomBaseDir = firstHeaderValue(H, {'randomBaseDir'}, nan);
parentFileName = firstHeaderValue(H, {'fileName'}, fileNameExt);
parentProbeDirectionsDeg = firstHeaderValue(H, {'probeDirectionsDeg','probeOffsetsDeg'}, []);
parentNProbeDirections = firstHeaderValue(H, {'nProbeDirections','nProbeOffsets'}, nan);
prefCohNoisePC = firstHeaderValue(H, {'prefCohNoisePC','noiseAmplitudePC','cohNoisePC'}, nan);

if isfield(S, 'sessionProbeHeader')
  P = S.sessionProbeHeader;
  sourceDomain = firstHeaderValue(P, {'sourceDomain'}, sourceDomain);
  randomBaseDir = firstHeaderValue(P, {'randomBaseDir'}, randomBaseDir);
  if isempty(parentProbeDirectionsDeg)
    parentProbeDirectionsDeg = firstHeaderValue(P, ...
      {'parentProbeDirectionsDeg','parentProbeOffsetsDeg'}, []);
  end
  if ~isfinite(parentNProbeDirections)
    parentNProbeDirections = firstHeaderValue(P, ...
      {'parentNProbeDirections','parentNProbeOffsets'}, nan);
  end
  if ~isfinite(nNoiseTrials)
    nNoiseTrials = firstHeaderValue(P, {'nNoiseTrials'}, nan);
  end
end
parentIsSingleProbe = parentNProbeDirections == 1;
parentIsInterleavedProbe = parentNProbeDirections > 1;
if isfield(S, 'sessionProbeHeader')
  probeDirDeg = firstHeaderValue(S.sessionProbeHeader, ...
    {'probeDirDeg','probeOffsetDeg'}, nan);
  probeTag = firstHeaderValue(S.sessionProbeHeader, {'probeTag'}, '');
  probeCohNoisePC = firstHeaderValue(S.sessionProbeHeader, ...
    {'probeCohNoisePC','noiseAmplitudePC'}, nan);
else
  probeCohNoisePC = nan;
  probeDirDeg = nan;
  probeTag = "";
end

fileInfoRow = table( ...
{filePath}, {folder}, {fileNameExt}, {animal}, {sourceDomain}, double(randomBaseDir), ...
probeDirDeg, {probeTag}, parentNProbeDirections, ...
{parentProbeDirectionsDeg}, parentIsSingleProbe, parentIsInterleavedProbe, ...
prefCohNoisePC, probeCohNoisePC, nTrials, nNoiseTrials, {parentFileName}, ...
'VariableNames', {'filePath','folder','fileName', 'animal', 'sourceDomain', 'randomBaseDir', ...
        'probeDirDeg','probeTag','parentNProbeDirections', ...
        'parentProbeDirectionsDeg', 'parentIsSingleProbe','parentIsInterleavedProbe', 'prefCohNoisePC', ...
        'probeCohNoisePC','nTrials','nNoiseTrials','parentFileName'});
fileInfoRow.twoAllocationMode = NaN;
if isfield(H,'twoAllocationMode')
  value = firstHeaderValue(H,{'twoAllocationMode'},NaN);
  assert(isscalar(value) && (isnumeric(value) || islogical(value)) && ...
    ismember(double(value),[0 1]), 'selectAnalysisFiles:IDLAllocationMode', ...
    'Invalid twoAllocationMode in %s.',filePath);
  fileInfoRow.twoAllocationMode = double(value);
end
end

%% -------------------------------------------------------------------------
function value = firstHeaderValue(H, candidateNames, defaultValue)
% Return the first available field, unwrapping Lablib's .data containers.

value = defaultValue;
for k = 1:numel(candidateNames)
  name = candidateNames{k};
  if ~isfield(H, name)
    continue
  end
  value = H.(name);
  while isstruct(value) && isfield(value, 'data')
    value = value.data;
  end
  if isstring(value) && isscalar(value)
    value = char(value);
  end
  return
end
end
