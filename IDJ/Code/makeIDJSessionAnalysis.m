function sessionAnalysis = makeIDJSessionAnalysis(sourcePath)
% makeIDJSessionAnalysis  Extract one compact IDJ FullSession analysis.
%
% sessionAnalysis = makeIDJSessionAnalysis(sourcePath)
%
% Loads sessionHeader and trials from a converted IDJ FullSession. The
% returned structure contains all trial-level quantities currently needed
% by downstream IDJ analyses, together with frame-resolved coherence noise
% for the full 750-ms prestep and 250-ms step. It intentionally does not
% retain the source trials. Flagged top-up trials are always excluded
% before extracting trial fields, configurations, or any stream noise.

S = load(sourcePath, 'sessionHeader', 'trials');
sessionHeader = S.sessionHeader;
trials = S.trials;
% Normalize the compact metadata contract for sessions converted before the
% current makeIDJSessionHeader implementation was installed.
sessionHeader.version = 1;
sessionHeader.sourceDomain = 'IDJ';
assert(double(sessionHeader.numberOfTrials) == numel(trials), ...
  'makeIDJSessionAnalysis:TrialCountMismatch', ...
  'sessionHeader.numberOfTrials does not match the number of recorded trials.');

[trialTable, validTrialIndex, configurationTable, stepTable, nExcludedTopUpTrials] = ...
  extractTrialTable(trials, sessionHeader);
[changeNoise, noChangeNoise, stepTimeMS, fullChangeNoise, ...
  fullNoChangeNoise, fullTimeMS] = ...
  extractStepNoise(trials(validTrialIndex), sessionHeader);

sessionAnalysis = struct();
sessionAnalysis.version = 2;
sessionAnalysis.topUpTrialsExcluded = true;
sessionAnalysis.nExcludedTopUpTrials = nExcludedTopUpTrials;
sessionAnalysis.sourceSessionPath = sourcePath;
sessionAnalysis.sourceSessionID = string(sessionHeader.fileName);
sessionAnalysis.createdBy = mfilename;
sessionAnalysis.createdAt = datetime('now');
sessionAnalysis.sessionHeader = sessionHeader;
sessionAnalysis.nSignalStreams = 3;
sessionAnalysis.streamTable = makeStreamTable(trialTable);
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
function [T, validTrialIndex, configurationTable, stepTable, nExcludedTopUpTrials] = ...
    extractTrialTable(trials, sessionHeader)

nTotal = numel(trials);
analysisTrial = false(nTotal, 1);
nExcludedTopUpTrials = 0;
for k = 1:nTotal
  certified = scalarData(trials{k}.trialCertify) == 0;
  eot = scalarData(trials{k}.extendedEOT);
  if certified && ismember(eot, [0 1])
    D = trials{k}.trial.data;
    assert(isfield(D,'topUpTrial') && isscalar(D.topUpTrial) && ...
      ismember(double(D.topUpTrial),[0 1]), ...
      'makeIDJSessionAnalysis:TopUpFlag', ...
      'Analysis-eligible trial %d lacks a valid topUpTrial flag.',k);
    topUp = logical(D.topUpTrial);
    nExcludedTopUpTrials = nExcludedTopUpTrials + double(topUp);
    analysisTrial(k) = ~topUp;
  end
end
validTrialIndex = find(analysisTrial);
n = numel(validTrialIndex);

sourceTrialIndex = validTrialIndex(:);
sourceTrialNumber = nan(n, 1);
correct = false(n, 1);
chosenSideIndex = nan(n, 1);
stepSideIndex = nan(n, 1);
reactionTimeMS = nan(n, 1);
allocIndex = nan(n, 1);
alloc0 = nan(n, 1);
alloc1 = nan(n, 1);
alloc2 = nan(n, 1);
minorAllocPC = nan(n, 1);
topUpTrial = false(n, 1);
unsignedAllocPC = nan(n, 1);
highStreamIndex = nan(n, 1);
totalStepIndex = nan(n, 1);
totalStepCohPC = nan(n, 1);
prestepCohPC = nan(n, 1);
stepEndCoh0PC = nan(n, 1);
stepEndCoh1PC = nan(n, 1);
stepEndCoh2PC = nan(n, 1);
signalDir0Deg = nan(n, 1);
signalDir1Deg = nan(n, 1);
signalDir2Deg = nan(n, 1);
stepHasNoise = false(n, 1);
noiseCohPC = nan(n, 1);
requestedRepsPerBlock = nan(n, 1);
configurationEpoch = nan(n, 1);

previousConfiguration = [];
epoch = 0;
for row = 1:n
  tr = trials{validTrialIndex(row)};
  D = tr.trial.data;
  B = tr.blockStatus.data(1);

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
  alloc2(row) = double(D.alloc2);
  [nominal, highStreamIndex(row)] = canonicalAllocation( ...
    [alloc0(row) alloc1(row) alloc2(row)]);
  unsignedAllocPC(row) = 100 * max(nominal);
  minorAllocPC(row) = 100 * min(nominal);
  topUpTrial(row) = logical(D.topUpTrial);
  totalStepIndex(row) = double(D.totalStepIndex) + 1;
  totalStepCohPC(row) = double(D.totalStepCohPC);
  prestepCohPC(row) = double(D.prestepCohPC);
  stepEndCoh0PC(row) = double(D.stepEndCoh0PC);
  stepEndCoh1PC(row) = double(D.stepEndCoh1PC);
  stepEndCoh2PC(row) = double(D.stepEndCoh2PC);
  signalDir0Deg(row) = mod(double(D.signalDir0Deg), 360);
  signalDir1Deg(row) = mod(double(D.signalDir1Deg), 360);
  signalDir2Deg(row) = mod(double(D.signalDir2Deg), 360);
  assert(double(B.numSignalStreams) == 3, ...
    'makeIDJSessionAnalysis:StreamCount', 'Expected three signal streams.');
  stepHasNoise(row) = logical(D.stepHasNoise);
  noiseCohPC(row) = double(D.noiseCohPC);

  nAllocationSlots = numel(B.allocEnabled);
  nStepSlots = numel(B.totalStepCohPC) / nAllocationSlots;
  assert(nStepSlots == round(nStepSlots), ...
    'makeIDJSessionAnalysis:BlockArrayShape', ...
    'Could not infer the number of step slots in blockStatus.');
  blockIndex = double(D.allocIndex) * nStepSlots + ...
    double(D.totalStepIndex) + 1;
  assert(abs(double(B.totalStepCohPC(blockIndex)) - totalStepCohPC(row)) < 1e-5, ...
    'makeIDJSessionAnalysis:TrialBlockStepMismatch', ...
    'The trial step does not match its blockStatus configuration.');
  requestedRepsPerBlock(row) = double(B.numRepsEachBlock(blockIndex));
end

tol = 1e-5;
assert(all(abs(alloc0 + alloc1 + alloc2 - 1) < tol), ...
  'makeIDJSessionAnalysis:AllocationSum', ...
  'Signal allocations do not sum to one on every analysis trial.');
assert(all(abs((stepEndCoh0PC - prestepCohPC) - alloc0 .* totalStepCohPC) < tol), ...
  'makeIDJSessionAnalysis:Stream0Endpoint', ...
  'Stream-0 endpoint coherence is inconsistent with allocation and total step.');
assert(all(abs((stepEndCoh1PC - prestepCohPC) - alloc1 .* totalStepCohPC) < tol), ...
  'makeIDJSessionAnalysis:Stream1Endpoint', ...
  'Stream-1 endpoint coherence is inconsistent with allocation and total step.');
assert(all(abs((stepEndCoh2PC - prestepCohPC) - alloc2 .* totalStepCohPC) < tol), ...
  'makeIDJSessionAnalysis:Stream2Endpoint', ...
  'Stream-2 endpoint coherence is inconsistent with allocation and total step.');
separations = sort([mod(signalDir1Deg - signalDir0Deg, 360), ...
  mod(signalDir2Deg - signalDir0Deg, 360)], 2);
assert(all(abs(separations - [120 240]) < tol, 'all'), ...
  'makeIDJSessionAnalysis:SignalSeparation', ...
  'Signal streams must be separated by 120 degrees.');
assert(n > 0, 'makeIDJSessionAnalysis:NoAnalysisTrials', ...
  'No certified correct/wrong non-top-up trials remain in this session.');
% Validate allocations trial by trial; partial sessions need not contain all four.

T = table(sourceTrialIndex, sourceTrialNumber, correct, ...
  chosenSideIndex, stepSideIndex, reactionTimeMS, allocIndex, alloc0, alloc1, alloc2, ...
  unsignedAllocPC, minorAllocPC, highStreamIndex, topUpTrial, totalStepIndex, totalStepCohPC, ...
  prestepCohPC, stepEndCoh0PC, stepEndCoh1PC, stepEndCoh2PC, ...
  signalDir0Deg, signalDir1Deg, signalDir2Deg, ...
  stepHasNoise, noiseCohPC, requestedRepsPerBlock, configurationEpoch);

configurationTable = makeConfigurationTable(T);
stepTable = makeStepTable(trials(validTrialIndex), T);
end

%% ------------------------------------------------------------------------
function signature = blockStatusSignature(B)

signature = [ ...
  double(B.randomizeDriftDir), double(B.stepHasNoise), ...
  double(B.numSignalStreams), double(B.cohNoiseFrameMS), ...
  double(B.signalSepDeg), double(B.noiseCohPC), ...
  double(B.numSteps), double(B.totalStepCohPC(:))', ...
  double(B.numAllocs), double(B.allocEnabled(:))', ...
  double(B.alloc0(:))', double(B.alloc1(:))', double(B.alloc2(:))', ...
  double(B.prestepCohPC), double(B.stepEndCoh0PC(:))', ...
  double(B.stepEndCoh1PC(:))', double(B.stepEndCoh2PC(:))', ...
  double(B.topUpFraction), double(B.numRepsEachBlock(:))' ...
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
    'makeIDJSessionAnalysis:BlockArrayShape', ...
    'Could not infer the number of step slots in blockStatus.');

  for allocationIndex = 1:double(B.numAllocs)
    if ~logical(B.allocEnabled(allocationIndex))
      continue
    end
    for stepIndex = 1:double(B.numSteps)
      blockIndex = (allocationIndex - 1) * nStepSlots + stepIndex;
      a = [double(B.alloc0(allocationIndex)), ...
        double(B.alloc1(allocationIndex)), double(B.alloc2(allocationIndex))];
      nominal = canonicalAllocation(a);
      endpoints = [double(B.stepEndCoh0PC(blockIndex)), ...
        double(B.stepEndCoh1PC(blockIndex)), double(B.stepEndCoh2PC(blockIndex))];
      % Endpoint "high" follows the nominal major stream. Equal thirds have
      % no major stream, so average their tiny recorded rounding differences.
      if max(nominal) - min(nominal) < 1e-6
        highEnd = mean(endpoints);
        lowEnd = highEnd;
      else
        highEnd = endpoints(nominal == max(nominal));
        lowEnd = mean(endpoints(nominal < max(nominal)));
      end
      configuredRows(end + 1,:) = [ ... %#ok<AGROW>
        epochs(epochIndex), 100 * max(nominal), stepIndex, ...
        double(B.totalStepCohPC(blockIndex)), highEnd, lowEnd, ...
        double(B.numRepsEachBlock(blockIndex))];
    end
  end
end

[U, ~, group] = unique(round(configuredRows, 6), 'rows');
nAllocationConfigurations = accumarray(group, 1);
% IDJ has one configured form per allocation, unlike IDL signed pairs.
nAnalysisTrials = zeros(size(U,1), 1);
for row = 1:size(U,1)
  trialMask = T.configurationEpoch == U(row,1) & ...
    abs(T.unsignedAllocPC - U(row,2)) < 1e-6 & ...
    T.totalStepIndex == U(row,3) & ...
    abs(T.totalStepCohPC - U(row,4)) < 1e-6;
  nAnalysisTrials(row) = sum(trialMask);
end

S = array2table([U nAllocationConfigurations nAnalysisTrials], 'VariableNames', { ...
  'configurationEpoch', 'unsignedAllocPC', 'totalStepIndex', ...
  'totalStepCohPC', 'highStepEndCohPC', 'lowStepEndCohPC', ...
  'requestedRepsPerBlock', 'nAllocationConfigurations', 'nAnalysisTrials'});
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
  'makeIDJSessionAnalysis:StepFrameCount', ...
  'The reconstructed step does not contain the expected number of video frames.');

nTrials = numel(trials);
nStepFrames = numel(stepFrames);
fullChangeNoise = nan(nFrames, nTrials, 3, 'single');
fullNoChangeNoise = nan(nFrames, nTrials, 3, 'single');
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
fields = {[prefix 'Noise0CohsPC'], [prefix 'Noise1CohsPC'], ...
  [prefix 'Noise2CohsPC']};
values = nan(numel(stepTimeMS), 1, 3, 'single');
timeToleranceMS = 1e-6;
for stream = 1:3
  cohs = double(tr.(fields{stream}).data(:));
  assert(numel(cohs) == numel(times), ...
    'makeIDJSessionAnalysis:NoiseEventCount', ...
    '%s times and stream-%d coherence events have different lengths.', ...
    prefix, stream - 1);
  for frame = 1:numel(stepTimeMS)
    eventIndex = find(times <= stepTimeMS(frame) + timeToleranceMS, 1, 'last');
    assert(~isempty(eventIndex), 'makeIDJSessionAnalysis:MissingNoiseValue', ...
      'No %s stream-%d noise value precedes step frame %d.', ...
      prefix, stream - 1, frame);
    values(frame,1,stream) = single(cohs(eventIndex));
  end
end
end

%% ------------------------------------------------------------------------
function value = scalarData(event)

value = double(event.data);
assert(isscalar(value), 'makeIDJSessionAnalysis:ExpectedScalarEvent', ...
  'Expected a scalar trial event.');
end

%% ------------------------------------------------------------------------
function [nominal, highStreamIndex] = canonicalAllocation(raw)
% Keep recorded fractions in trialTable; use nominal fractions for grouping.
patterns = [1 0 0; 0.7 0.15 0.15; 0.5 0.25 0.25; 1/3 1/3 1/3];
[distance, row] = min(max(abs(sort(patterns, 2) - sort(raw)), [], 2));
assert(all(isfinite(raw)) && distance < 1e-4, ...
  'makeIDJSessionAnalysis:UnexpectedAllocation', ...
  'Unexpected recorded allocation: %s.', mat2str(raw));
if row == 4
  nominal = [1/3 1/3 1/3];
  highStreamIndex = 0;
else
  [~, highStreamIndex] = max(raw);
  nominal = repmat(patterns(row, 2), 1, 3);
  nominal(highStreamIndex) = patterns(row, 1);
end
end

%% ------------------------------------------------------------------------
function S = makeStreamTable(T)
% Each stream contributes one row. Use this table for stream-noise pooling,
% never to duplicate observations in behavioral likelihood fits.
n = height(T);
analysisTrialRow = repelem((1:n)', 3);
sourceTrialIndex = T.sourceTrialIndex(analysisTrialRow);
streamIndex = repmat((1:3)', n, 1);
recordedAllocation = reshape([T.alloc0 T.alloc1 T.alloc2]', [], 1);
allocationPC = nan(3*n, 1);
role = strings(3*n, 1);
for trial = 1:n
  rows = (trial - 1) * 3 + (1:3);
  [a, major] = canonicalAllocation([T.alloc0(trial) T.alloc1(trial) T.alloc2(trial)]);
  allocationPC(rows) = 100 * a(:);
  if major == 0
    role(rows) = "equal";
  else
    role(rows) = "minor";
    role(rows(major)) = "major";
  end
end
S = table(analysisTrialRow, sourceTrialIndex, streamIndex, ...
  recordedAllocation, allocationPC, role);
end
