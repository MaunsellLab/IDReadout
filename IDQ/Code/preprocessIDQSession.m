function preprocessIDQSession(varargin)
% preprocessIDQSession
%
% Build processed IDQ session files from converted FullSessions files.
%
% Input files live in:
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

% ---- File selection criteria ----
selectArgs = { ...
  % 'MinTrials', 500, ...
  };

% ---- Run options ----
p = inputParser;
p.addParameter('Replace', false, @islogical);
p.addParameter('Verbose', false, @islogical);
p.parse(varargin{:});
opts = p.Results;

domainPath = domainFolder(mfilename('fullpath'));
fullSessionFolder = fullfile(domainPath, 'Data', 'FullSessions');
processedSessionFolder = validFolder(fullfile(domainPath, 'Data', 'ProcessedSessions'));
files = selectIDQFiles(fullSessionFolder, selectArgs{:});
if opts.Verbose
  fprintf('preprocessIDQSession: %d files selected\n', numel(files));
end
for iFile = 1:numel(files)

  inPath = files(iFile).path;
  [~, baseName] = fileparts(inPath);
  outPath = fullfile(processedSessionFolder, [baseName '.mat']);

  if isfile(outPath) && ~opts.Replace
    if opts.Verbose
      fprintf('  %3d/%3d  exists, skipping: %s\n', iFile, numel(files), files(iFile).name);
    end
    continue
  end
  if opts.Verbose
    fprintf('  %3d/%3d  processing: %s\n', iFile, numel(files), files(iFile).name);
  end
  S = load(inPath, 'header', 'trials');
  header = S.header;
  trials = S.trials;
  sessionHeader = makeIDQSessionHeader(header, trials);
  trialData = extractIDQTrialData(trials);
  noiseBySideDir = extractIDQNoiseBySideDir(trials, trialData, sessionHeader);
  sessionHeader.nTrials = height(trialData);
  save(outPath, 'header', 'sessionHeader', 'trialData', 'noiseBySideDir', '-v7.3');
end

makeIDQSessionAnalyses();
plotIDQSessionSummaries();

fprintf('preprocessIDQSession complete\n');

end

function noiseBySideDir = extractIDQNoiseBySideDir(trials, trialData, sessionHeader)
% extractIDQNoiseBySideDir
%
% Extract full six-stream IDQ noise for valid completed trials.
%
% Output:
%   noiseBySideDir: 2 sides x 3 directions x nFrames x nTrials
%
% Side with coherence increment step:
%   1 = RF, 2 = Opp
% Direction dimension:
%   1, 2, 3 corresponding to sessionHeader.dirDeg
%
% Noise is kept in raw task sign convention. DEC trials are not flipped.


msPerVFrame = 1000.0 / sessionHeader.frameRateHz;
m = round((sessionHeader.preStepMS + sessionHeader.stepMS) / msPerVFrame);
nTrials = numel(trialData.validIdx);
noiseBySideDir = nan(2, 3, m, nTrials);
t = trials(trialData.validIdx);

for kk = 1:nTrials
  tr = t{kk};

  changeNoise = fillFromTimes('change', tr, m, msPerVFrame);
  noChangeNoise = fillFromTimes('noChange', tr, m, msPerVFrame);

  if tr.trial.data.changeSide == 0
    % RF changed, Opp noChange
    noiseBySideDir(1, :, :, kk) = changeNoise;
    noiseBySideDir(2, :, :, kk) = noChangeNoise;
  else
    % Opp changed, RF noChange
    noiseBySideDir(1, :, :, kk) = noChangeNoise;
    noiseBySideDir(2, :, :, kk) = changeNoise;
  end
end
end

%%=========================================================================
function noise = fillFromTimes(fieldName, tr, m, msPerVFrame)

noise0 = tr.([fieldName, 'CohsPC0']).data;
noise1 = tr.([fieldName, 'CohsPC1']).data;
noise2 = tr.([fieldName, 'CohsPC2']).data;

% In a few of the early data files collected with IDQ, the noiseCohPC was
% saved as an integer, rounding the value down to 7 or -7.  We correct that
% here. It is a small correction, and had negligible effect on NLL of fits

% This correction was wrong -- the noise was actually presented with 
% integer roundings

% if abs(noise0(1)) == 7
%   correction = 10 / sqrt(2) / 7;
%   noise0 = noise0 * correction;
%   noise1 = noise1 * correction;
%   noise2 = noise2 * correction;
% end

timesMS = tr.([fieldName, 'TimesMS']).data;
nTimes = min(numel(timesMS), numel(noise0));
timesMS = timesMS(1:nTimes);
noise0  = noise0(1:nTimes);
noise1  = noise1(1:nTimes);
noise2  = noise2(1:nTimes);

noise = nan(1, 3, m);

for tIndex = 1:nTimes
  t0 = timesMS(tIndex);
  theVFrame = max(1, floor(t0 / msPerVFrame) + 1);
  if theVFrame > m
    continue;
  end
  if tIndex < nTimes
    t1 = timesMS(tIndex + 1);
    nextVFrame = floor(t1 / msPerVFrame) + 1;
  else
    nextVFrame = m + 1;
  end
  if nextVFrame <= theVFrame
    continue;
  end
  if nextVFrame > m + 1
    nextVFrame = m + 1;
  end
  noise(:, 1, theVFrame:nextVFrame - 1) = noise0(tIndex);
  noise(:, 2, theVFrame:nextVFrame - 1) = noise1(tIndex);
  noise(:, 3, theVFrame:nextVFrame - 1) = noise2(tIndex);
end
end

function trialData = extractIDQTrialData(trials)
% extractIDQTrialData
%
% Extract analysis-facing trialData from raw converted IDQ trials.
%
% Raw trials may use task-native 0-based indexing.
% trialData uses MATLAB/IDR-facing 1-based indexing.
%
% Conventions:
%   trialData.correct         logical
%   trialData.stepSignIndex   1 = DEC, 2 = INC
%   trialData.sideIndex       1 = RF, 2 = Opp
%   trialData.chosenSideIndex 1 = RF, 2 = Opp
%   trialData.dirIndex        1, 2, 3

eotCode = getTrialField(trials, 'extendedEOT', 'data');
certify = getTrialField(trials, 'trialCertify', 'data');
validMask = (eotCode == 0 | eotCode == 1) & (certify == 0);
trialData.validIdx = find(validMask);
t = trials(trialData.validIdx);
trialData.correct = eotCode(trialData.validIdx) == 0;
trialData.trialCertify = getTrialField(t, 'trialCertify', 'data');
trialData.chosenSideIndex = getTrialField(t, 'targetChosen', 'data') + 1;
trialData.stepSignIndex = getTrialField(t, 'trial', 'data', 'stepIndex') + 1;
trialData.sideIndex = getTrialField(t, 'trial', 'data', 'changeSide') + 1;
trialData.dirIndex = getTrialField(t, 'trial', 'data', 'dirIndex') + 1;
trialData.stepCoh = getTrialField(t, 'trial', 'data', 'stepCohPC');
trialData.hasStepNoise = logical(getTrialField(t, 'trial', 'data', 'stepHasNoise'));
trialData.sideBias = getTrialField(t, 'RFBias', 'data');
trialData.dirBias = getTrialField(t, 'dirBias', 'data');

end

