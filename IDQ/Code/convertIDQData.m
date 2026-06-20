function convertIDQData()
%   Convert Lablib .dat files to .mat format
%   Sources are in .../IDQ/Data/DatFiles
%   Output to .../IDQ/Data/FullSessions
%
%   For each  .dat file in dataFolder, this function checks for a
%   corresponding *_fileInfo.mat file. If it does not exist, the .dat
%   file is read with readLLFile and the trials/header are saved.

[dataFolder, existed] = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'DatFiles'));
if ~existed 
  fprintf('  convertIDQData -- failed to find data to convert in %s', dataFolder);
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

  % Save header, sessionHeader, and trials; use -v7.3 if these can get large
  sessionHeader = makeSessionHeader(header);
  save(outFileName, 'trials', 'header', 'sessionHeader');
  fprintf('      Saved %s\n', outFileName);
end
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
  'baseDirDeg', ...
  'numNoiseStreams', ...
  'cohNoisePC' ...
  };

for k = 1:numel(copyFields)
  f = copyFields{k};
  if isfield(header, f)
    sessionHeader.(f) = localDataValue(header.(f));
  end
end
sessionHeader.dirsDeg = header.dirsDeg.data(1:sessionHeader.numNoiseStreams)';
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