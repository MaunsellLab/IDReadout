function sessionHeader = makeIDQSessionHeader(header, trials)
% makeIDQSessionHeader
%
% Build lightweight sessionHeader for processed IDQ session files.
%
% Keep this fast to load. Do not place trialwise noise arrays here.

eotCode = getTrialField(trials, 'extendedEOT', 'data');
keep = eotCode == 0 | eotCode == 1;
t = trials(keep);
sessionHeader = struct();

% ---- Identity / provenance ----
if isfield(header, 'fileName')
  sessionHeader.fileName = header.fileName;
end

if isfield(header, 'animal')
  sessionHeader.animal = header.animal;
end

if isfield(header, 'taskVersion')
  sessionHeader.taskVersion = header.taskVersion;
end

% ---- Trial counts ----
sessionHeader.nTrialsTotal = numel(trials);
sessionHeader.nTrials = numel(t);   % valid completed trials stored in processed file

% ---- Timing ----
sessionHeader.frameRateHz = header.frameRateHz.data(1);
sessionHeader.preStepMS = header.preStepMS.data(1);
sessionHeader.stepMS = header.stepMS.data(1);
sessionHeader.stepDurationMS = header.stepMS.data(1);

% Use IDR-style naming if useful downstream.
frameRateHz = sessionHeader.frameRateHz;
msPerVFrame = 1000.0 / frameRateHz;
nNoiseFrames = round((sessionHeader.preStepMS + sessionHeader.stepMS) / msPerVFrame);

sessionHeader.msPerVFrame = msPerVFrame;
sessionHeader.nNoiseFrames = nNoiseFrames;
sessionHeader.tMS = (0:nNoiseFrames-1) * msPerVFrame;

stepStartMS = sessionHeader.preStepMS;
stepStopMS = sessionHeader.preStepMS + sessionHeader.stepMS;

sessionHeader.stepWindowFrameIndex = ...
  find(sessionHeader.tMS >= stepStartMS & sessionHeader.tMS < stepStopMS);

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
