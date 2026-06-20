function preprocessIDQSession(varargin)
% preprocessIDQSession
%
% Build processed IDQ session files from converted FullSessions files.
%
% Converted files live in:
%   Data/FullSessions
%
% Processed files are written to:
%   Data/ProcessedSessions
%
% Processed files contain:
%   header
%   sessionHeader
%   trialData
%   noiseBySideDir
%
% Raw 'trials' is not propagated.

% ---- File selection criteria: keep visible and centralized ----
selectArgs = { ...
  % Future examples:
  % 'MinTrials', 500, ...
  % 'TaskVersion', '...', ...
  };

% ---- Run options ----
p = inputParser;
p.addParameter('Replace', false, @islogical);
p.parse(varargin{:});

replace = p.Results.Replace;

domainPath = domainFolder(mfilename('fullpath'));

fullSessionFolder = fullfile(domainPath, 'Data', 'FullSessions');
processedSessionFolder = fullfile(domainPath, 'Data', 'ProcessedSessions');

if ~isfolder(processedSessionFolder)
  mkdir(processedSessionFolder);
end

files = selectIDQFiles(fullSessionFolder, selectArgs{:});

fprintf('preprocessIDQSession: %d files selected\n', numel(files));

for iFile = 1:numel(files)

  inPath = files(iFile).path;
  [~, baseName] = fileparts(inPath);
  outPath = fullfile(processedSessionFolder, [baseName '.mat']);

  if isfile(outPath) && ~replace
    fprintf('  %3d/%3d  exists, skipping: %s\n', ...
      iFile, numel(files), files(iFile).name);
    continue
  end

  fprintf('  %3d/%3d  processing: %s\n', ...
    iFile, numel(files), files(iFile).name);

  S = load(inPath, 'header', 'trials');

  header = S.header;
  trials = S.trials;

  sessionHeader = makeIDQSessionHeader(header, trials);
  trialData = extractIDQTrialData(trials);
  noiseBySideDir = extractIDQNoiseBySideDir(trials, trialData, sessionHeader);

  sessionHeader.nTrials = height(trialData);
  save(outPath, 'header', 'sessionHeader', 'trialData', 'noiseBySideDir', '-v7.3');
end

fprintf('preprocessIDQSession complete\n');

end