function sessionAnalysis = makeIDLSessionAnalysis(sourcePath)
% makeIDLSessionAnalysis  Extract one compact IDL FullSession analysis.
%
% sessionAnalysis = makeIDLSessionAnalysis(sourcePath)
%
% Loads sessionHeader and trials from a converted IDL FullSession. The
% returned structure contains all trial-level quantities currently needed
% by downstream IDL analyses, together with frame-resolved coherence noise
% for the full 750-ms prestep and 250-ms step. It intentionally does not
% retain the source trials.

S = load(sourcePath, 'sessionHeader', 'trials');
sessionHeader = S.sessionHeader;
trials = S.trials;
% Normalize the compact metadata contract for sessions converted before the
% current makeIDLSessionHeader implementation was installed.
sessionHeader = idlSetAllocationHeader(sessionHeader,trials);
sessionHeader.sourceDomain = 'IDL';
assert(double(sessionHeader.numberOfTrials) == numel(trials), ...
  'makeIDLSessionAnalysis:TrialCountMismatch', ...
  'sessionHeader.numberOfTrials does not match the number of recorded trials.');

[trialTable, validTrialIndex, configurationTable, stepTable] = ...
  extractTrialTable(trials, sessionHeader);
[changeNoise, noChangeNoise, stepTimeMS, fullChangeNoise, ...
  fullNoChangeNoise, fullTimeMS] = ...
  extractStepNoise(trials(validTrialIndex), sessionHeader);

sessionAnalysis = struct();
sessionAnalysis.version = 2;
sessionAnalysis.sourceSessionPath = sourcePath;
sessionAnalysis.sourceSessionID = string(sessionHeader.fileName);
sessionAnalysis.createdBy = mfilename;
sessionAnalysis.createdAt = datetime('now');
sessionAnalysis.sessionHeader = sessionHeader;
sessionAnalysis.nTrialsTotal = numel(trials);
sessionAnalysis.nAnalysisTrials = height(trialTable);
sessionAnalysis.nExcludedTrials = numel(trials) - height(trialTable);
sessionAnalysis.nCorrect = sum(trialTable.correct);
sessionAnalysis.nWrong = sum(~trialTable.correct);
sessionAnalysis.stepTimeMS = stepTimeMS;
sessionAnalysis.trialTable = trialTable;
sessionAnalysis.configurationTable = configurationTable;
sessionAnalysis.stepTable = stepTable;
sessionAnalysis.changeNoiseByFrameTrial = changeNoise;
sessionAnalysis.noChangeNoiseByFrameTrial = noChangeNoise;
sessionAnalysis.fullTimeMS = fullTimeMS;
sessionAnalysis.fullChangeNoiseByFrameTrial = fullChangeNoise;
sessionAnalysis.fullNoChangeNoiseByFrameTrial = fullNoChangeNoise;
sessionAnalysis.noiseArrayContract = ...
  '[step frame, analysis trial, signal stream]; coherence percentage points';
sessionAnalysis.fullNoiseArrayContract = ...
  '[prestep/step frame, analysis trial, signal stream]; coherence percentage points';
end

%% ------------------------------------------------------------------------
function [T, validTrialIndex, configurationTable, stepTable] = ...
    extractTrialTable(trials, sessionHeader)

nTotal = numel(trials);
analysisTrial = false(nTotal, 1);
for k = 1:nTotal
  certified = scalarData(trials{k}.trialCertify) == 0;
  eot = scalarData(trials{k}.extendedEOT);
  D = trials{k}.trial.data;
  topUp = isfield(D,'topUpTrial') && logical(D.topUpTrial);
  analysisTrial(k) = certified && ismember(eot, [0 1]) && ~topUp;
end
validTrialIndex = find(analysisTrial);
n = numel(validTrialIndex);
assert(n>0,'makeIDLSessionAnalysis:NoAnalysisTrials', ...
  'No certified, non-top-up correct/error trials remain.');

sourceTrialIndex = validTrialIndex(:);
sourceTrialNumber = nan(n, 1);
correct = false(n, 1);
chosenSideIndex = nan(n, 1);
stepSideIndex = nan(n, 1);
reactionTimeMS = nan(n, 1);
allocIndex = nan(n, 1);
alloc0 = nan(n, 1);
alloc1 = nan(n, 1);
unsignedAllocPC = nan(n, 1);
highStreamIndex = nan(n, 1);
totalStepIndex = nan(n, 1);
totalStepCohPC = nan(n, 1);
prestepCohPC = nan(n, 1);
stepEndCoh0PC = nan(n, 1);
stepEndCoh1PC = nan(n, 1);
signalDir0Deg = nan(n, 1);
signalDir1Deg = nan(n, 1);
stepHasNoise = false(n, 1);
noiseCohPC = nan(n, 1);
requestedRepsPerBlock = nan(n, 1);
configurationEpoch = nan(n, 1);
splitUnsplitRatio = nan(n, 1);

previousConfiguration = [];
epoch = 0;
for row = 1:n
  tr = trials{validTrialIndex(row)};
  D = tr.trial.data;
  B = tr.blockStatus.data(1);
  info = idlAllocationConfiguration(B);
  splitUnsplitRatio(row) = info.splitUnsplitRatio;

  configuration = blockStatusSignature(B);
  if row == 1 || ~isequaln(configuration, previousConfiguration)
    epoch = epoch + 1;
    previousConfiguration = configuration;
  end
  configurationEpoch(row) = epoch;

  sourceTrialNumber(row) = double(tr.trialNumber);
  correct(row) = scalarData(tr.extendedEOT) == 0;
  chosenSideIndex(row) = scalarData(tr.targetChosen) + 1;
  stepSideIndex(row) = double(D.stepSide) + 1;
  reactionTimeMS(row) = scalarData(tr.reactTimeMS);
  allocIndex(row) = double(D.allocIndex) + 1;
  alloc0(row) = double(D.alloc0);
  alloc1(row) = double(D.alloc1);
  unsignedAllocPC(row) = 100 * max(alloc0(row), alloc1(row));
  if alloc0(row) > alloc1(row)
    highStreamIndex(row) = 1;
  elseif alloc1(row) > alloc0(row)
    highStreamIndex(row) = 2;
  else
    highStreamIndex(row) = 0;
  end
  totalStepIndex(row) = double(D.totalStepIndex) + 1;
  totalStepCohPC(row) = double(D.totalStepCohPC);
  prestepCohPC(row) = double(D.prestepCohPC);
  stepEndCoh0PC(row) = double(D.stepEndCoh0PC);
  stepEndCoh1PC(row) = double(D.stepEndCoh1PC);
  signalDir0Deg(row) = mod(double(D.signalDir0Deg), 360);
  signalDir1Deg(row) = mod(double(D.signalDir1Deg), 360);
  stepHasNoise(row) = logical(D.stepHasNoise);
  noiseCohPC(row) = double(D.noiseCohPC);

  nAllocationSlots = numel(B.allocEnabled);
  nStepSlots = numel(B.totalStepCohPC) / nAllocationSlots;
  assert(nStepSlots == round(nStepSlots), ...
    'makeIDLSessionAnalysis:BlockArrayShape', ...
    'Could not infer the number of step slots in blockStatus.');
  blockIndex = double(D.allocIndex) * nStepSlots + ...
    double(D.totalStepIndex) + 1;
  assert(abs(double(B.totalStepCohPC(blockIndex)) - totalStepCohPC(row)) < 1e-5, ...
    'makeIDLSessionAnalysis:TrialBlockStepMismatch', ...
    'The trial step does not match its blockStatus configuration.');
  requestedRepsPerBlock(row) = double(B.numRepsEachBlock(blockIndex));
end

tol = 1e-5;
assert(all(abs(alloc0 + alloc1 - 1) < tol), ...
  'makeIDLSessionAnalysis:AllocationSum', ...
  'Signal allocations do not sum to one on every analysis trial.');
assert(all(abs((stepEndCoh0PC - prestepCohPC) - alloc0 .* totalStepCohPC) < tol), ...
  'makeIDLSessionAnalysis:Stream0Endpoint', ...
  'Stream-0 endpoint coherence is inconsistent with allocation and total step.');
assert(all(abs((stepEndCoh1PC - prestepCohPC) - alloc1 .* totalStepCohPC) < tol), ...
  'makeIDLSessionAnalysis:Stream1Endpoint', ...
  'Stream-1 endpoint coherence is inconsistent with allocation and total step.');
assert(all(mod(signalDir1Deg - signalDir0Deg, 360) == 90), ...
  'makeIDLSessionAnalysis:SignalSeparation', ...
  'Signal streams are not separated by 90 degrees on every analysis trial.');

observedAllocations = unique(round(unsignedAllocPC, 6))';
expectedAllocations = sort(sessionHeader.enabledAllocationPC);
assert(isequal(observedAllocations, expectedAllocations), ...
  'makeIDLSessionAnalysis:UnexpectedAllocations', ...
  'Expected unsigned allocations %s; found %s.', ...
  mat2str(expectedAllocations), mat2str(observedAllocations));

T = table(sourceTrialIndex, sourceTrialNumber, correct, ...
  chosenSideIndex, stepSideIndex, reactionTimeMS, allocIndex, alloc0, alloc1, ...
  unsignedAllocPC, highStreamIndex, totalStepIndex, totalStepCohPC, ...
  prestepCohPC, stepEndCoh0PC, stepEndCoh1PC, signalDir0Deg, signalDir1Deg, ...
  stepHasNoise, noiseCohPC, requestedRepsPerBlock, configurationEpoch,splitUnsplitRatio);

configurationTable = makeConfigurationTable(T);
stepTable = makeStepTable(trials(validTrialIndex), T);
end

%% ------------------------------------------------------------------------
function signature = blockStatusSignature(B)

info = idlAllocationConfiguration(B);
signature = [ ...
  info.splitUnsplitRatio, ...
  double(B.randomizeDriftDir), double(B.stepHasNoise), ...
  double(B.numSignalStreams), double(B.cohNoiseFrameMS), ...
  double(B.signalSepDeg), double(B.noiseCohPC), ...
  double(B.numSteps), double(B.totalStepCohPC(:))', ...
  double(B.numAllocs), double(B.allocEnabled(:))', ...
  double(B.alloc0(:))', double(B.alloc1(:))', ...
  double(B.prestepCohPC), double(B.stepEndCoh0PC(:))', ...
  double(B.stepEndCoh1PC(:))', double(B.numRepsEachBlock(:))' ...
  ];
signature = round(signature, 6);
end

%% ------------------------------------------------------------------------
function C = makeConfigurationTable(T)

epochs = unique(T.configurationEpoch, 'stable');
n = numel(epochs);
configurationEpoch = epochs(:);
firstSourceTrial = nan(n, 1);
lastSourceTrial = nan(n, 1);
nAnalysisTrials = nan(n, 1);
for k = 1:n
  idx = T.configurationEpoch == epochs(k);
  firstSourceTrial(k) = min(T.sourceTrialNumber(idx));
  lastSourceTrial(k) = max(T.sourceTrialNumber(idx));
  nAnalysisTrials(k) = sum(idx);
end
C = table(configurationEpoch, firstSourceTrial, lastSourceTrial, nAnalysisTrials);
C.splitUnsplitRatio = arrayfun(@(e) T.splitUnsplitRatio( ...
  find(T.configurationEpoch==e,1)),epochs(:));
end

%% ------------------------------------------------------------------------
function S = makeStepTable(validTrials, T)

configuredRows = [];
epochs = unique(T.configurationEpoch, 'stable');
for epochIndex = 1:numel(epochs)
  trialRow = find(T.configurationEpoch == epochs(epochIndex), 1, 'first');
  B = validTrials{trialRow}.blockStatus.data(1);
  nAllocationSlots = numel(B.allocEnabled);
  nStepSlots = numel(B.totalStepCohPC) / nAllocationSlots;
  assert(nStepSlots == round(nStepSlots), ...
    'makeIDLSessionAnalysis:BlockArrayShape', ...
    'Could not infer the number of step slots in blockStatus.');

  for allocationIndex = 1:double(B.numAllocs)
    if ~logical(B.allocEnabled(allocationIndex))
      continue
    end
    for stepIndex = 1:double(B.numSteps)
      blockIndex = (allocationIndex - 1) * nStepSlots + stepIndex;
      a0 = double(B.alloc0(allocationIndex));
      a1 = double(B.alloc1(allocationIndex));
      end0 = double(B.stepEndCoh0PC(blockIndex));
      end1 = double(B.stepEndCoh1PC(blockIndex));
      configuredRows(end + 1,:) = [ ... %#ok<AGROW>
        epochs(epochIndex), 100 * max(a0, a1), stepIndex, ...
        double(B.totalStepCohPC(blockIndex)), max(end0, end1), ...
        min(end0, end1), double(B.numRepsEachBlock(blockIndex))];
    end
  end
end

[U, ~, group] = unique(round(configuredRows, 6), 'rows');
nSignedAllocations = accumarray(group, 1);
assert(all(nSignedAllocations == 2), ...
  'makeIDLSessionAnalysis:AsymmetricSignedConfigurations', ...
  'Each unsigned allocation-step configuration must have two signed forms.');
nAnalysisTrials = zeros(size(U,1), 1);
for row = 1:size(U,1)
  trialMask = T.configurationEpoch == U(row,1) & ...
    abs(T.unsignedAllocPC - U(row,2)) < 1e-6 & ...
    T.totalStepIndex == U(row,3) & ...
    abs(T.totalStepCohPC - U(row,4)) < 1e-6;
  nAnalysisTrials(row) = sum(trialMask);
end

S = array2table([U nSignedAllocations nAnalysisTrials], 'VariableNames', { ...
  'configurationEpoch', 'unsignedAllocPC', 'totalStepIndex', ...
  'totalStepCohPC', 'highStepEndCohPC', 'lowStepEndCohPC', ...
  'requestedRepsPerBlock', 'nSignedAllocations', 'nAnalysisTrials'});
S = sortrows(S, {'configurationEpoch','unsignedAllocPC','totalStepIndex'});
end

%% ------------------------------------------------------------------------
function [changeNoise, noChangeNoise, stepTimeMS, fullChangeNoise, ...
    fullNoChangeNoise, relativeFullTimeMS] = ...
    extractStepNoise(trials, sessionHeader)

frameRateHz = double(sessionHeader.frameRateHz);
preStepMS = double(sessionHeader.preStepMS);
stepMS = double(sessionHeader.stepMS);
frameMS = 1000 / frameRateHz;
nFrames = round((preStepMS + stepMS) / frameMS);
fullTimeMS = (0:nFrames - 1)' * frameMS;
stepFrames = find(fullTimeMS >= preStepMS & fullTimeMS < preStepMS + stepMS);
stepTimeMS = fullTimeMS(stepFrames) - preStepMS;
assert(numel(stepFrames) == round(stepMS / frameMS), ...
  'makeIDLSessionAnalysis:StepFrameCount', ...
  'The reconstructed step does not contain the expected number of video frames.');

nTrials = numel(trials);
nStepFrames = numel(stepFrames);
fullChangeNoise = nan(nFrames, nTrials, 2, 'single');
fullNoChangeNoise = nan(nFrames, nTrials, 2, 'single');
for k = 1:nTrials
  fullChangeNoise(:,k,:) = fillSide(trials{k}, 'change', fullTimeMS);
  fullNoChangeNoise(:,k,:) = fillSide(trials{k}, 'noChange', fullTimeMS);
end
changeNoise = fullChangeNoise(stepFrames,:,:);
noChangeNoise = fullNoChangeNoise(stepFrames,:,:);
relativeFullTimeMS = fullTimeMS - preStepMS;
end

%% ------------------------------------------------------------------------
function values = fillSide(tr, prefix, stepTimeMS)

times = double(tr.([prefix 'TimesMS']).data(:));
fields = {[prefix 'Noise0CohsPC'], [prefix 'Noise1CohsPC']};
values = nan(numel(stepTimeMS), 1, 2, 'single');
timeToleranceMS = 1e-6;
for stream = 1:2
  cohs = double(tr.(fields{stream}).data(:));
  assert(numel(cohs) == numel(times), ...
    'makeIDLSessionAnalysis:NoiseEventCount', ...
    '%s times and stream-%d coherence events have different lengths.', ...
    prefix, stream - 1);
  for frame = 1:numel(stepTimeMS)
    eventIndex = find(times <= stepTimeMS(frame) + timeToleranceMS, 1, 'last');
    assert(~isempty(eventIndex), 'makeIDLSessionAnalysis:MissingNoiseValue', ...
      'No %s stream-%d noise value precedes step frame %d.', ...
      prefix, stream - 1, frame);
    values(frame,1,stream) = single(cohs(eventIndex));
  end
end
end

%% ------------------------------------------------------------------------
function value = scalarData(event)

value = double(event.data);
assert(isscalar(value), 'makeIDLSessionAnalysis:ExpectedScalarEvent', ...
  'Expected a scalar trial event.');
end
