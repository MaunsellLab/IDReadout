
function NEWpreprocessIDQSession(varargin)
% preprocessIDQSession
%
% Build processed IDQ session files from converted FullSessions files.
%
% Converted files read from: Data/FullSessions
% Processed files written to: Data/ProcessedSessions
%
% Processed files contain:
%   header
%   sessionHeader
%   trialData
%   noiseBySideDir

% ---- Run options ----
p = inputParser;
p.addParameter('Replace', true, @islogical);
p.parse(varargin{:});

replace = p.Results.Replace;
domainPath = domainFolder(mfilename('fullpath'));
selectArgs = {};      % Examples: 'MinTrials', 500,  'TaskVersion', ...
fullSessionFolder = fullfile(domainPath, 'Data', 'NEWFullSessions');
processedSessionFolder = validFolder(fullfile(domainPath, 'Data', 'NEWProcessedSessions'));
files = selectIDQFiles(fullSessionFolder, selectArgs{:});
fprintf('preprocessIDQSession: %d files selected\n', numel(files));

for fIndex = 1:numel(files)
  inPath = files(fIndex).path;
  [~, baseName] = fileparts(inPath);
  outPath = fullfile(processedSessionFolder, [baseName '.mat']);
  if isfile(outPath) && ~replace
    continue
  end

  fprintf('       processing %s\n', baseName);
  load(inPath, 'header', 'trials');

  [sessionHeader, trialData] = makeIDQSessionHeaders(header, trials);
  % trialData = extractIDQTrialData(trials);
  noiseBySideDir = extractIDQNoiseBySideDir(trials, trialData, sessionHeader);


  % MakeSessionHeader already gets the number of trials?
  % sessionHeader.nTrials = height(trialData);
  save(outPath, 'header', 'sessionHeader', 'trialData', 'noiseBySideDir', '-v7.3');
end

% makeIDQSessionAnalyses();
% makeIDQSessionSummaries();

fprintf('preprocessIDQSession complete\n');

end

%%=========================================================================
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

noiseBySideDir = nan(2, 3, sessionHeader.nNoiseVF, sessionHeader.nTrials);
t = trials(trialData.validIdx);
for tIndex = 1:sessionHeader.nTrials
  tr = t{tIndex};
  changeNoise = fillFromTimes('change', tr, sessionHeader);
  noChangeNoise = fillFromTimes('noChange', tr, sessionHeader);
  if tr.trial.data.changeSide == 0
    % RF changed, Opp noChange
    noiseBySideDir(1, :, :, tIndex) = changeNoise;
    noiseBySideDir(2, :, :, tIndex) = noChangeNoise;
  else
    % Opp changed, RF noChange
    noiseBySideDir(1, :, :, tIndex) = noChangeNoise;
    noiseBySideDir(2, :, :, tIndex) = changeNoise;
  end
end
end

%%=========================================================================
function noise = fillFromTimes(fieldName, tr, sessionHeader)
% return the 3 noise streams for one RDK (fieldName) on one trial (tr).
% Noises are specified by video frame (nVF) after converting from trial time
% in ms using msPerVF.

nVF = sessionHeader.nNoiseVF;
msPerVF = sessionHeader.msPerVF;
timesMS = tr.([fieldName, 'TimesMS']).data;
nTimes = numel(timesMS);
timeStepMS = mean(diff(timesMS));              % size of steps in the time vector
noisePeriodMS = 3 * timeStepMS;                % duration of 
nPhaseVF = round(noisePeriodMS / msPerVF);
cohStr = {'CohsPC0', 'CohsPC1', 'CohsPC2'};
inNoise = nan(1, 3, nTimes);
for nIndex = 1:3
  inNoise(1, nIndex, :) = tr.([fieldName, cohStr{nIndex}]).data;
  % inNoise(1, 2, :) = tr.([fieldName, 'CohsPC1']).data;
  % inNoise(1, 3, :) = tr.([fieldName, 'CohsPC2']).data;
  nPhase(nIndex) = noisePhase(inNoise(:, nIndex, :), timesMS, nPhaseVF, msPerVF);
end
nanIdx = isnan(nPhase);
if sum(nanIdx) >= 2
  error('Cannot process to interdeterminate noise phases');
elseif sum(nanIdx == 1) 
  currentPhase = nPhase(~nanIdx);
  if all(ismember(currentPhase, [0 2 4]))
    validVals = [0 2 4];
  elseif all(ismember(currentPhase, [1 3 5]))
    validVals = [1 3 5];
  else
    error('fillMissingPhaseValue:InvalidPhaseValues','Phase values must belong to either [0 2 4] or [1 3 5].');
  end
  missingVal = setdiff(validVals, currentPhase);
  if numel(missingVal) ~= 1
    error('fillMissingPhaseValue:CannotInferMissingValue', 'Could not infer a unique missing value.');
  end
  fprintf('nan value on %s %s side\n', sessionHeader.fileName, fieldName);
  nPhase(nanIdx) = missingVal;
end
noise = nan(1, 3, nVF);
for tIndex = 1:nTimes
  VFIndex = round(timesMS(tIndex) / msPerVF) + 1;
  if VFIndex > nVF
    continue;
  end
  VFPhase = mod(round(timesMS(tIndex) / msPerVF), nPhaseVF);
  for nIndex = 1:3
    if nPhase(nIndex) == VFPhase
      noise(:, nIndex, VFIndex) = inNoise(:, nIndex, tIndex);
      % fprintf('time %ld (%.1f ms), VFIndex %d noise %d (phase %d), noise(time) %.0f\n',...
      %   tIndex, timesMS(tIndex), VFIndex, nIndex, nPhase(nIndex), inNoise(:, nIndex, tIndex));
    end
  end
end
end

%%=========================================================================
function phase = noisePhase(noise, timesMS, nPhaseVF, msPerVF)

changeIndex = find(diff(noise) ~= 0, 1, 'first') + 1;
if isempty(changeIndex)
  phase = nan;
  return;
end
changeTimeMS = timesMS(changeIndex);
changeVFIndex = round(changeTimeMS / msPerVF);
phase = mod(changeVFIndex, nPhaseVF);
end

%%=========================================================================
function [sessionHeader, trialData] = makeIDQSessionHeaders(header, trials)
% makeIDQSessionHeader
%
% Build lightweight sessionHeader for processed IDQ session files.
% Fast to load: does not include noise arrays here.

eotCode = getTrialField(trials, 'extendedEOT', 'data');
certify = getTrialField(trials, 'trialCertify', 'data');
validMask = (eotCode == 0 | eotCode == 1) & (certify == 0);

% construct the trialData struc
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

%construct the sessionHeader struct
sessionHeader = struct();
% ---- Identity / provenance ----
if isfield(header, 'fileName')
  sessionHeader.fileName = header.fileName;
end

if isfield(header, 'taskVersion')
  sessionHeader.taskVersion = header.taskVersion;
end

% ---- Trial counts ----
sessionHeader.nTrialsTotal = numel(trials);
sessionHeader.nTrials = numel(t);

% ---- Timing ----
sessionHeader.frameRateHz = header.frameRateHz.data(1);
sessionHeader.preStepMS = header.preStepMS.data(1);
sessionHeader.stepMS = header.stepMS.data(1);
sessionHeader.stepDurationMS = header.stepMS.data(1);

% Use IDR-style naming if useful downstream.
msPerVF = 1000.0 / sessionHeader.frameRateHz;
nNoiseVF = round((sessionHeader.preStepMS + sessionHeader.stepMS) / msPerVF);

sessionHeader.msPerVF = msPerVF;
sessionHeader.nNoiseVF = nNoiseVF;
sessionHeader.tMS = (0:nNoiseVF-1) * msPerVF;
cohNoiseFrameMS = [header.blockStatus.data.cohNoiseFrameMS];
sessionHeader.noiseFrameMS = cohNoiseFrameMS(1);

stepStartMS = sessionHeader.preStepMS;
stepStopMS = sessionHeader.preStepMS + sessionHeader.stepMS;

sessionHeader.stepWindowFrameIndex = find(sessionHeader.tMS >= stepStartMS & sessionHeader.tMS < stepStopMS);

sessionHeader.nDirs = 3;
sessionHeader.dirsDeg = header.dirsDeg.data(1:sessionHeader.nDirs);
sessionHeader.noiseAmplitudePC = header.cohNoisePC.data(1);
sessionHeader.nCohSteps = header.blockStatus.data(1).numSteps;
sessionHeader.stepCohPC = header.blockStatus.data(1).stepCohPC(1:sessionHeader.nCohSteps)';

% ---- Bias summaries ----
sideBias = getTrialField(t, 'RFBias', 'data');
sessionHeader.sideBiasMedian = median(sideBias, 'omitnan');
[~, sessionHeader.sideBiasIQR] = iqr(sideBias);
mag = getTrialField(t, 'dirBias', 'data', 'magnitude');
sessionHeader.dirBiasMagMean = mean(mag, 'omitnan');
sessionHeader.dirBiasMagSD = std(mag, 'omitnan');
dirDeg = getTrialField(t, 'dirBias', 'data', 'dirDeg');
sessionHeader.dirBiasDirDegMean = mean(dirDeg, 'omitnan');
sessionHeader.dirBiasDirDegSD = std(dirDeg, 'omitnan');
end
