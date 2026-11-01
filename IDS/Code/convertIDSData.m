function convertIDSData()
%   Convert Lablib .dat files to .mat format
%   Sources are in .../IDS/Data/DatFiles
%   Output to .../IDS/Data/FullSessions
%
%   For each  .dat file in dataFolder, this function checks for a
%   corresponding *_fileInfo.mat file. If it does not exist, the .dat
%   file is read with readLLFile and the trials/header are saved.

[dataFolder, existed] = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'DatFiles'));
if ~existed 
  fprintf('  convertIDSData -- failed to find data to convert in %s', dataFolder);
  return;
end
[names, paths] = getDatFileList(dataFolder);
convertedFolder = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'FullSessions'));

% Convert any unconverted .dat files
for fi = 1:numel(paths)
  datPath = string(paths(fi));
  datName = string(names(fi));
  [~, baseName] = fileparts(datPath);
  tempInfoFileName = sprintf('%s/%s_fileInfo.mat', dataFolder, baseName);
  infoFileName = sprintf('%s/%s_fileInfo.mat', convertedFolder, baseName); % .mat headerfile name for .dat file
  outFileName = sprintf('%s/%s.mat', convertedFolder, baseName);           % .mat data file name for .dat file
  if isfile(infoFileName) && isfile(outFileName)                    % Skip if *_fileInfo.mat already exists
    continue;
  end
  if isfile(infoFileName) && ~isfile(outFileName)                   % Skip only *_fileInfo.mat exists, re-read
    delete(infoFileName);
  end
  fprintf('Reading %s\n', datName);
  header = readLLFile('i', datPath);
  movefile(tempInfoFileName, infoFileName);   % readLLFile leaves the info in working directory
  nTrials = header.numberOfTrials;
  trials  = cell(1, nTrials);
  fprintf('       converting %s (%d trials)\n', datName, nTrials);
  for t = 1:nTrials
    trials{t} = readLLFile('t', t);
    if t == nTrials || mod(t, 500) == 0                            % Occasionally update the waitbar text
      fprintf('       reading %s: trial %d of %d (%.0f%%)\n', datName, t, nTrials, 100 * t / nTrials);
    end
  end

  % correct mis-assignment of stepNoisePC that occurred in some early IDS dat files.
  d = datetime(strtrim(string(header.date)), 'InputFormat', 'MMMM d, yyyy');
  if d < datetime(2026, 8, 12) 
    [header, trials] = correctStepNoisePC(header, trials);
  end

  % Save header, sessionHeader, and trials; use -v7.3 if these can get large
  sessionHeader = makeSessionHeader(header);
  save(outFileName, 'trials', 'header', 'sessionHeader');
  fprintf('      Saved %s\n', outFileName);
end
end

% Early IDS versions added cohNoisePC to stepCohPC when constructing the
% displayed drift-direction step, even though stepCohPC was already an
% endpoint coherence. Correct the derived trial records to describe the
% coherences actually displayed; the primary .dat files remain unchanged.

function [header, trials] = correctStepNoisePC(header, trials)

correctionName = 'IDSStepCoherenceDoubleBase';
correctionVersion = 1;

if isfield(header, 'dataCorrections') && isfield(header.dataCorrections, correctionName)
  error('convertIDSData:CorrectionAlreadyApplied', ...
    'Correction %s has already been applied.', correctionName);
end
nCorrected = 0;
nSkipped = 0;
pedestalValuesPC = nan(1, numel(trials));
for t = 1:numel(trials)
  if ~isfield(trials{t}, 'trial')
    nSkipped = nSkipped + 1;
    continue;
  end

  % readLLFile normally wraps event values in a scalar struct with a data
  % field. Also accept an unwrapped trial descriptor for compatibility.
  trialEvent = trials{t}.trial;
  wrappedEvent = isstruct(trialEvent) && isscalar(trialEvent) && isfield(trialEvent, 'data');
  if wrappedEvent
    trialDesc = trialEvent.data;
  else
    trialDesc = trialEvent;
  end

  if ~isstruct(trialDesc) || ~isscalar(trialDesc) || ~isfield(trialDesc, 'cohNoisePC')
    error('convertIDSData:InvalidTrialDescriptor', ...
      'Trial %d does not contain a scalar trial descriptor with cohNoisePC.', t);
  end

  % Support both the current TrialDesc field names and the older analysis
  % aliases. A matched pair is required so that endpoint and baseline retain
  % the same convention.
  if isfield(trialDesc, 'baseNoisePC') && isfield(trialDesc, 'stepNoisePC')
    baseField = 'baseNoisePC';
    stepField = 'stepNoisePC';
  elseif isfield(trialDesc, 'baseCohPC') && isfield(trialDesc, 'stepCohPC')
    baseField = 'baseCohPC';
    stepField = 'stepCohPC';
  else
    error('convertIDSData:MissingCoherenceFields', ...
      'Trial %d lacks a matching base/step coherence field pair.', t);
  end

  % cohNoisePC is 7.071% on ordinary trials and zero on instruction trials.
  % Using the trial value avoids altering instruction trials and also handles
  % any affected sessions collected with a different noise amplitude.
  pedestalPC = trialDesc.cohNoisePC;
  trialDesc.(baseField) = pedestalPC;
  trialDesc.(stepField) = trialDesc.(stepField) + pedestalPC;
  pedestalValuesPC(t) = double(pedestalPC);
  nCorrected = nCorrected + 1;

  if wrappedEvent
    trials{t}.trial.data = trialDesc;
  else
    trials{t}.trial = trialDesc;
  end
end

correction = struct();
correction.applied = true;
correction.correctionVersion = correctionVersion;
correction.appliedDate = datetime('now');
correction.nTrialsCorrected = nCorrected;
correction.nTrialsSkippedWithoutTrialEvent = nSkipped;
correction.pedestalValuesPC = unique(pedestalValuesPC(~isnan(pedestalValuesPC)));
correction.baseFormula = 'basePC = trial.cohNoisePC';
correction.stepFormula = 'stepPC = rawStepPC + trial.cohNoisePC';
correction.reason = ['Affected IDS versions added cohNoisePC to an endpoint ' ...
  'that already included the coherence pedestal.'];

if ~isfield(header, 'dataCorrections') || isempty(header.dataCorrections)
  header.dataCorrections = struct();
end
header.dataCorrections.(correctionName) = correction;
end

%%------------------------------------------------
function sessionHeader = makeSessionHeader(header)
% makeSessionHeader  Build stable session-level metadata for analyses.
%
% sessionHeader is the analysis-facing metadata inherited from the converted
% Lablib header. It intentionally contains only session-level quantities. 

sessionHeader = struct();

% Stable fields used by downstream analyses or provenance displays.
copyFields = { ...
  'date', ...
  'subject', ...
  'taskName', ...
  'numberOfTrials', ...
  'frameRateHz', ...
  'preStepMS', ...
  'stepMS', ...
  'numNoiseStreams', ...
  'cohNoisePC', ...
  'numProbeOffsets', ...
  'dataCorrections' ...
  };

for k = 1:numel(copyFields)
  f = copyFields{k};
  if isfield(header, f)
    sessionHeader.(f) = localDataValue(header.(f));
  end
end
subStr = {'Meetz', 'Neesha'};
tf = cellfun(@(s) contains(header.fileName, s), subStr);
if sum(tf) == 1
  sessionHeader.animal = subStr{tf};
else
  sessionHeader.animal = 'unknown';
end
[~, sessionHeader.fileName] = fileparts(header.fileName);
sessionHeader.createdBy = mfilename;
sessionHeader.createdDate = datetime('now');
sessionHeader.metadataContract = 'session-level analysis metadata';
end

%% -------------------------------------------------------------------------
function v = localDataValue(x)

v = x;
while isstruct(v) && isfield(v, 'data')
  v = v.data;
end
if isnumeric(v) && ~isscalar(v)
  v = v(1);
end
end

%% -------------------------------------------------------------------------
function [names, paths] = getDatFileList(dataFolder)
S = dir(fullfile(dataFolder, '*.dat'));
names = {S.name}';
paths = fullfile(dataFolder, names);
[names, idx] = sort(names);
paths = paths(idx);
end
