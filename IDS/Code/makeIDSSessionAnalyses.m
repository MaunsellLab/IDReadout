function sessionAnalysisAll = makeIDSSessionAnalyses(varargin)
% makeIDSSessionAnalyses  Build compact IDS analyses from probe sessions.

p = inputParser;
p.addParameter('Replace', false, @islogical);
p.addParameter('Verbose', true, @islogical);
p.parse(varargin{:});
opts = p.Results;

root = domainFolder(mfilename('fullpath'));
probeDirs = dir(fullfile(root, 'Data', 'Probe*'));
probeDirs = probeDirs([probeDirs.isdir]);
sessionAnalysisAll = {};
for iDir = 1:numel(probeDirs)
  probeTag = probeDirs(iDir).name;
  sourceFolder = fullfile(probeDirs(iDir).folder, probeTag, 'ProbeSessions');
  if ~isfolder(sourceFolder)
    continue
  end
  dataFilePaths = selectAnalysisFiles(sourceFolder);
  analysisFolder = validFolder(fullfile(root, 'Data', probeTag, 'SessionAnalysis'));
  for iFile = 1:numel(dataFilePaths)
    sourcePath = dataFilePaths{iFile};
    [~, baseName] = fileparts(sourcePath);
    outPath = fullfile(analysisFolder, [baseName '_sessionAnalysis.mat']);
    if isfile(outPath) && ~opts.Replace
      vars = who('-file', outPath);
      if all(ismember({'sessionAnalysis','sessionHeader','sessionProbeHeader'}, vars))
        continue
      end
    end
    sessionAnalysis = makeOneAnalysis(sourcePath);
    % Top-level headers make compact analysis files selectable by the common
    % selectAnalysisFiles utility without loading the larger analysis struct.
    sessionHeader = sessionAnalysis.sessionHeader; %#ok<NASGU>
    sessionProbeHeader = sessionAnalysis.sessionProbeHeader; %#ok<NASGU>
    save(outPath, 'sessionAnalysis', 'sessionHeader', 'sessionProbeHeader', '-v7.3');
    sessionAnalysisAll{end+1,1} = sessionAnalysis; %#ok<AGROW>
    if opts.Verbose
      fprintf('makeIDSSessionAnalyses: saved %s\n', outPath);
    end
  end
end
if isscalar(sessionAnalysisAll)
  sessionAnalysisAll = sessionAnalysisAll{1};
end
end

function A = makeOneAnalysis(sourcePath)
S = load(sourcePath, 'header', 'sessionHeader', 'sessionProbeHeader', ...
  'probeTrials', 'trialIdx');
P = S.sessionProbeHeader;
[trialData, validLocalIdx] = extractTrialData(S.probeTrials, S.trialIdx, P.probeOffsetDeg);
[noiseBySideStream, tMS, stepFrames] = extractNoise(S.probeTrials(validLocalIdx), trialData, S.header);

[changeStream, noChangeStream] = sideAlignedNoise(noiseBySideStream, trialData.sideIndex);
changeSum = sum(changeStream, 3, 'omitnan');
noChangeSum = sum(noChangeStream, 3, 'omitnan');
fullFrames = (1:numel(tMS))';

T = makeTrialTable(trialData, P);
names = {'Drift','Plus','Minus'};
for k = 1:3
  T.(['rectStepChangeNoisePred' names{k}]) = mean(changeStream(stepFrames,:,k), 1, 'omitnan')';
  T.(['rectFullChangeNoisePred' names{k}]) = mean(changeStream(fullFrames,:,k), 1, 'omitnan')';
  T.(['rectStepNoChangeNoisePred' names{k}]) = mean(noChangeStream(stepFrames,:,k), 1, 'omitnan')';
  T.(['rectFullNoChangeNoisePred' names{k}]) = mean(noChangeStream(fullFrames,:,k), 1, 'omitnan')';
end

noisyStepCoh = unique(T.stepCoh(T.hasStepNoise));
if numel(noisyStepCoh) ~= 1
  error('makeIDSSessionAnalyses:NoisyStepMismatch', ...
    'Expected one noisy step magnitude in %s; found %s.', P.sessionID, mat2str(noisyStepCoh'));
end

A = struct();
A.fileName = sprintf('%s_%s', P.sessionID, P.probeTag);
A.sourceSessionID = P.sourceSessionID;
A.sourceProbeSession = sourcePath;
A.createdBy = mfilename;
A.createdAt = datetime('now');
A.header = S.header;
A.sessionHeader = S.sessionHeader;
A.sessionProbeHeader = P;
A.nTrials = height(T);
A.nStepNoiseTrials = sum(T.hasStepNoise);
A.noisyStepCoh = noisyStepCoh;
A.tMS = tMS;
A.stepFrames = stepFrames;
A.trialTable = T;
A.noiseBySideStream = noiseBySideStream;
A.changeStreamNoiseByFrameTrial = changeStream;
A.noChangeStreamNoiseByFrameTrial = noChangeStream;
A.changeSumNoiseByFrameTrial = changeSum;
A.noChangeSumNoiseByFrameTrial = noChangeSum;
A.sumNoiseByFrameTrial = changeSum;       % IDQ-compatible aliases
A.dirNoiseByFrameTrial = changeStream;
A.sumNoiseKernel = computeKernel(changeSum, T.correct, T.hasStepNoise);
A.noChangeSumNoiseKernel = computeKernel(noChangeSum, T.correct, T.hasStepNoise);
A.streamNames = {'Drift','Drift+Probe','Drift-Probe'};
A.streamDirectionDefinition = P.streamDirectionDefinition;
end

function [D, validLocalIdx] = extractTrialData(probeTrials, parentTrialIdx, expectedOffset)
n = numel(probeTrials);
keep = false(n,1);
for k = 1:n
  tr = probeTrials{k};
  if ~isfield(tr, 'trial') || ~isfield(tr.trial, 'data') || ...
      ~isfield(tr, 'extendedEOT') || ~isfield(tr, 'trialCertify')
    continue
  end
  td = tr.trial.data;
  keep(k) = ~logical(td.instructTrial) && double(tr.trialCertify.data) == 0 && ...
    ismember(double(tr.extendedEOT.data), [0 1]);
end
validLocalIdx = find(keep);
t = probeTrials(validLocalIdx);
m = numel(t);
D = struct();
D.parentTrialIdx = parentTrialIdx(validLocalIdx)';
D.correct = false(m,1); D.chosenSideIndex = nan(m,1); D.stepSignIndex = nan(m,1);
D.sideIndex = nan(m,1); D.stepIndex = nan(m,1); D.baseCoh = nan(m,1);
D.stepEndpointCoh = nan(m,1); D.stepCoh = nan(m,1); D.hasStepNoise = false(m,1);
D.driftDirDeg = nan(m,1); D.probeOffsetDeg = nan(m,1); D.offsetIndex = nan(m,1);
D.cohNoisePC = nan(m,1); D.streamDirDeg = nan(m,3);
for k = 1:m
  tr = t{k}; td = tr.trial.data;
  D.correct(k) = double(tr.extendedEOT.data) == 0;
  D.chosenSideIndex(k) = double(tr.targetChosen.data) + 1;
  D.stepSignIndex(k) = double(td.stepSign) + 1;
  D.sideIndex(k) = double(td.stepSide) + 1;
  D.stepIndex(k) = double(td.stepIndex) + 1;
  D.baseCoh(k) = double(td.baseCohPC);
  D.stepEndpointCoh(k) = double(td.stepCohPC);
  D.stepCoh(k) = abs(D.stepEndpointCoh(k) - D.baseCoh(k));
  D.hasStepNoise(k) = logical(td.stepHasNoise);
  D.driftDirDeg(k) = mod(double(td.driftDirDeg), 360);
  D.probeOffsetDeg(k) = double(td.probeOffsetDeg);
  D.offsetIndex(k) = double(td.offsetIndex) + 1;
  D.cohNoisePC(k) = double(td.cohNoisePC);
  D.streamDirDeg(k,:) = [D.driftDirDeg(k), ...
    mod(D.driftDirDeg(k) + D.probeOffsetDeg(k), 360), ...
    mod(D.driftDirDeg(k) - D.probeOffsetDeg(k), 360)];
end
assert(all(D.probeOffsetDeg == expectedOffset), ...
  'makeIDSSessionAnalyses:MixedProbeOffsets', 'Probe session contains mixed offsets.');
assert(all(abs(D.cohNoisePC - 10/sqrt(2)) < 1e-3), ...
  'makeIDSSessionAnalyses:UnexpectedNoiseAmplitude', 'Expected 7.071%% noise amplitude.');
end

function [noise, tMS, stepFrames] = extractNoise(trials, D, header)
frameRateHz = header.frameRateHz.data(1);
preStepMS = header.preStepMS.data(1);
stepMS = header.stepMS.data(1);
frameMS = 1000 / frameRateHz;
nFrames = round((preStepMS + stepMS) / frameMS);
tMS = (0:nFrames-1)' * frameMS;
stepFrames = find(tMS >= preStepMS & tMS < preStepMS + stepMS);
noise = nan(2, 3, nFrames, numel(trials));
for k = 1:numel(trials)
  tr = trials{k};
  change = fillSide(tr, 'change', nFrames, frameMS);
  noChange = fillSide(tr, 'noChange', nFrames, frameMS);
  if D.sideIndex(k) == 1
    noise(1,:,:,k) = change; noise(2,:,:,k) = noChange;
  else
    noise(1,:,:,k) = noChange; noise(2,:,:,k) = change;
  end
end
end

function values = fillSide(tr, prefix, nFrames, frameMS)
fields = {[prefix 'DriftCohsPC'], [prefix 'Probe0CohsPC'], [prefix 'Probe1CohsPC']};
times = double(tr.([prefix 'TimesMS']).data(:));
values = nan(1, 3, nFrames);
for stream = 1:3
  cohs = double(tr.(fields{stream}).data(:));
  n = min(numel(times), numel(cohs));
  for j = 1:n
    firstFrame = max(1, floor(times(j)/frameMS) + 1);
    if j < n
      nextFrame = floor(times(j+1)/frameMS) + 1;
    else
      nextFrame = nFrames + 1;
    end
    nextFrame = min(nextFrame, nFrames + 1);
    if firstFrame <= nFrames && nextFrame > firstFrame
      values(1,stream,firstFrame:nextFrame-1) = cohs(j);
    end
  end
end
end

function [change, noChange] = sideAlignedNoise(noise, sideIndex)
[~, nStreams, nFrames, nTrials] = size(noise);
change = nan(nFrames, nTrials, nStreams);
noChange = nan(nFrames, nTrials, nStreams);
for k = 1:nTrials
  c = sideIndex(k); nc = 3-c;
  change(:,k,:) = reshape(squeeze(noise(c,:,:,k))', nFrames, 1, nStreams);
  noChange(:,k,:) = reshape(squeeze(noise(nc,:,:,k))', nFrames, 1, nStreams);
end
end

function T = makeTrialTable(D, P)
n = numel(D.correct);
T = table();
T.sessionName = repmat(string(sprintf('%s_%s', P.sessionID, P.probeTag)), n, 1);
T.sourceSessionID = repmat(string(P.sourceSessionID), n, 1);
T.trialIndex = (1:n)';
T.parentTrialIdx = D.parentTrialIdx(:);
for f = {'correct','chosenSideIndex','stepSignIndex','sideIndex','stepIndex', ...
    'baseCoh','stepEndpointCoh','stepCoh','hasStepNoise','driftDirDeg', ...
    'probeOffsetDeg','offsetIndex','cohNoisePC'}
  T.(f{1}) = D.(f{1})(:);
end
T.plusProbeDirDeg = D.streamDirDeg(:,2);
T.minusProbeDirDeg = D.streamDirDeg(:,3);
end

function K = computeKernel(noise, correct, hasNoise)
use = logical(hasNoise(:)); correct = logical(correct(:));
K = struct();
K.meanCorrect = mean(noise(:,use & correct), 2, 'omitnan');
K.meanError = mean(noise(:,use & ~correct), 2, 'omitnan');
K.kernel = K.meanCorrect - K.meanError;
K.nCorrect = sum(use & correct); K.nError = sum(use & ~correct);
K.nTrials = K.nCorrect + K.nError;
end
