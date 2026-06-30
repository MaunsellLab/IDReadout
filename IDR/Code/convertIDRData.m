function convertIDRData()
%   Convert Lablib .dat files to .mat format
%   convertIDRData() uses the default folder specified by domainFolder(mfilename('fullpath')).
%   convertIDRData(dataFolder) uses the specified folder.
%
%   For each  .dat file in dataFolder, this function checks for a
%   corresponding *_fileInfo.mat file. If it does not exist, the .dat
%   file is read with readLLFile and the trials/header are saved.

[dataFolder, existed] = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'DatFiles'));
if ~existed 
  fprintf('  convertIDRData -- failed to find data to convert in %s', dataFolder);
  return;
end

[paths, ~] = selectAnalysisFiles(dataFolder, 'FilePattern', '*.dat');
if isempty(paths)
  fprintf('No dat files found\n');
  return;
end

convertedFolder = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'FullSessions'));

% Convert any unconverted .dat files
for fi = 1:numel(paths)
  datPath = string(paths(fi));
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
  fprintf('       reading %s.mat\n', baseName);
  header = readLLFile('i', datPath);
  movefile(tempInfoFileName, infoFileName);   % readLLFile leaves the info with the dat -- move it to the mat
  nTrials = header.numberOfTrials;
  trials  = cell(1, nTrials);

  subStr = {'Meetz', 'Neesha'};
  tf = cellfun(@(s) contains(header.fileName, s), subStr);
  if sum(tf) == 1
    trialMeta.animal = subStr{tf};
  else
    trialMeta.animal = 'unknown';
  end
  
  trialMeta.probeDirectionsDeg = [];
  trialMeta.nProbeDirections = 0;
  trialMeta.nNoiseTrials = 0;
  fprintf('       converting %s.mat\n', baseName);
  for t = 1:nTrials
    trials{t} = readLLFile('t', t);
    D = trials{t}.trial.data;
    if D.cohNoise
      trialMeta.nNoiseTrials = trialMeta.nNoiseTrials + 1;
      if isfield(D, 'probeDirDeg')
        probeDirDeg = double(D.probeDirDeg);
        if isfinite(probeDirDeg) && probeDirDeg ~= -1
          trialMeta.probeDirectionsDeg(end+1) = probeDirDeg; 
        end
      end
    end
    if t == nTrials || mod(t, 500) == 0                            % Occasionally update the waitbar text
      fprintf('       reading %s: trial %d of %d (%.0f%%)\n', baseName, t, nTrials, 100 * t / nTrials);
    end
  end

  probeDirs = unique(round(trialMeta.probeDirectionsDeg, 6));
  % Old-format fallback: no per-trial probeDirDeg, but parent header has one.
  if isempty(probeDirs) && isfield(header, 'probeDirDeg') && isfield(header.probeDirDeg, 'data')
    probeDirs = double(header.probeDirDeg.data);
  end
  trialMeta.probeDirectionsDeg = probeDirs(:)';
  trialMeta.nProbeDirections = numel(probeDirs);
  trialMeta.probeTags = cell(1, numel(trialMeta.nProbeDirections));
  for p = 1:trialMeta.nProbeDirections
    trialMeta.probeTags{p} = sprintf('Probe%d', trialMeta.probeDirectionsDeg(p));
  end

  % correct indexing errors that existed in earlier version of IDR
  trials = correctIndices(header, trials);

  % Save header, sessionHeader, and trials; use -v7.3 if these can get large
  sessionHeader = makeSessionHeader(header, trialMeta);
  save(outFileName, 'trials', 'header', 'sessionHeader');
  fprintf('      Saved %s\n', outFileName);
end
end

% Several events were mislabeled as "change*" and "noChange*" when they
% were in fact "RF*" and "Opp*". These include "*Dots",
% "*PrefCohsPC", "*ProbeCohsPC", and "*TimesMS".  Specifically, the RF location
% was always assigned to be the "change" side and the opposite location was always
% assigned to be the noChange side. Thus, whenever the change was on the opposite side,
% the assignments were backwards.  Here, we patch the change/noChange assignments
% to correct this error by reassigning the values in the affected files.

function trials = correctIndices(header, trials)

d = datetime(strtrim(string(header.date)), 'InputFormat', 'MMMM d, yyyy');
if d < datetime(2026, 2, 13) || d > datetime(2026, 3, 12)
  return;
else
  fprintf('Data from %s. Correcting trial-shift error.\n', header.date);
  nTrials  = numel(trials);
  for k = 1:nTrials
    tr = trials{k};
    if ~isfield(tr, 'trialEnd')
      continue;                       % skip if missing outcome
    end
    if k < nTrials
      patchTrial = trials{k+1};
    else
      patchTrial = tr;
    end
    % Trial shift without swapping
    tr.changeDots =           patchTrial.changeDots;
    tr.noChangeDots =         patchTrial.noChangeDots;
    tr.changePrefCohsPC =     patchTrial.changePrefCohsPC;
    tr.changeProbeCohsPC =    patchTrial.changeProbeCohsPC;
    tr.noChangePrefCohsPC =   patchTrial.noChangePrefCohsPC;
    tr.noChangeProbeCohsPC =  patchTrial.noChangeProbeCohsPC;
    tr.changeTimesMS =        patchTrial.changeTimesMS;
    tr.noChangeTimesMS =      patchTrial.noChangeTimesMS;
    trials{k} = tr;
  end
end
end

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
sessionHeader.animal = getMetaField(trialMeta, 'animal', []);
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