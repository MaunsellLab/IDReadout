function correction = runIDQRebalancingCorrection(varargin)
% runIDQRebalancingCorrection
%
% Reanalyse IDQ gains using conditional-expected direction occupancies under
% the historical LLRDKMixture rebalancing code.
%
% Place this file in the IDQ folder (beside plotIDQSummary.m) and run:
%
%   correction = runIDQRebalancingCorrection;
%
% Or specify folders explicitly:
%
%   correction = runIDQRebalancingCorrection( ...
%       'FullSessionFolder', '/path/to/IDQ/Data/FullSessions', ...
%       'OutputFolder', '/path/to/output');
%
% The function reads the raw FullSessions locally. It writes compact files:
%
%   IDQ_RebalancingCorrection.mat
%   IDQ_RebalancingCorrection_Gains.csv
%   IDQ_RebalancingCorrection_Predictors.csv
%
% The MAT output contains nominal and corrected gain fits, paired gain
% changes, nominal-versus-corrected predictor diagnostics, an approximate
% residual finite-dot uncertainty, and the six trial-level predictors. It
% does not contain the raw trials or frame-by-frame arrays.
%
% Required analysis functions already used by plotIDQSummary:
%   fitIDQInitialSessionThresholds
%   fitIDQAcrossAlignedWeibull
%   fitIDQNoiseGain
%   idqWeibullP
%   fitIDQSessionWeibullFixedBeta
%   idqWeibullAlphaForThreshold
%
% Historical behaviors represented here:
%   * probability groups, rather than individual directions, are balanced;
%   * target/source counter updates have the original short-circuit order;
%   * step frame indices are truncated and at most one due step is advanced
%     after each displayed frame;
%   * newborn directions use the old 1000-level LLRandom draw;
%   * dots retain age when transferred and expire at the recorded lifetime.
%
% The state propagation is a conditional-mean calculation. Transferred dots
% are removed proportionally from the source direction/age cells and assigned
% uniformly among the target group's directions. The unobserved departures
% from this conditional mean are summarized as residual uncertainty.

p = inputParser;
thisFolder = fileparts(mfilename('fullpath'));
p.addParameter('FullSessionFolder', fullfile(thisFolder, '../Data', 'FullSessions'), ...
    @(x) ischar(x) || isstring(x));
p.addParameter('OutputFolder', fullfile(thisFolder, '../Data', 'AcrossSessionSummaries'), ...
    @(x) ischar(x) || isstring(x));
p.addParameter('FilePattern', '*.mat', @(x) ischar(x) || isstring(x));
p.addParameter('Verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;

fullSessionFolder = char(opts.FullSessionFolder);
outputFolder = char(opts.OutputFolder);
if ~isfolder(fullSessionFolder)
    error('runIDQRebalancingCorrection:MissingFolder', ...
        'FullSessions folder does not exist: %s', fullSessionFolder);
end
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

files = dir(fullfile(fullSessionFolder, char(opts.FilePattern)));
files = files(~endsWith({files.name}, '_fileInfo.mat'));
if isempty(files)
    error('runIDQRebalancingCorrection:NoFiles', ...
        'No FullSession MAT files found in %s.', fullSessionFolder);
end
[~, order] = sort({files.name});
files = files(order);

if opts.Verbose
    fprintf('IDQ rebalancing correction: %d sessions\n', numel(files));
end

sessionData = cell(numel(files), 1);
sessionDiagnostics = repmat(emptySessionDiagnostic(), numel(files), 1);
for iSession = 1:numel(files)
    filePath = fullfile(files(iSession).folder, files(iSession).name);
    if opts.Verbose
        fprintf('  %2d/%2d  %s\n', iSession, numel(files), files(iSession).name);
    end
    S = load(filePath, 'header', 'trials');
    [sessionData{iSession}, sessionDiagnostics(iSession)] = ...
        reconstructOneSession(S.header, S.trials, iSession, files(iSession).name, opts.Verbose);
end

% The psychometric calibration is identical for the two gain fits.
trialTables = cellfun(@(s) s.trialTable, sessionData, 'UniformOutput', false);
trialTable = vertcat(trialTables{:});
initialBetaWeibull = 2;
initialLapse = 0.02;
targetPerformance = 0.75;
lapseBounds = [0 0.05];

initialPsych = fitIDQInitialSessionThresholds(trialTable, ...
    initialBetaWeibull, initialLapse, targetPerformance);
alignedWeibull = fitIDQAcrossAlignedWeibull(initialPsych.alignedCoh, ...
    trialTable.correct, targetPerformance, lapseBounds);

nominalLOO = computeLOOPredictors(sessionData, 'nominal');
correctedLOO = computeLOOPredictors(sessionData, 'corrected');

trialTable.nominalNoisePredDir1 = nominalLOO.noisePredDir1;
trialTable.nominalNoisePredDir2 = nominalLOO.noisePredDir2;
trialTable.nominalNoisePredDir3 = nominalLOO.noisePredDir3;
trialTable.correctedNoisePredDir1 = correctedLOO.noisePredDir1;
trialTable.correctedNoisePredDir2 = correctedLOO.noisePredDir2;
trialTable.correctedNoisePredDir3 = correctedLOO.noisePredDir3;

nominalFitTable = selectPredictors(trialTable, 'nominal');
correctedFitTable = selectPredictors(trialTable, 'corrected');
nominalGains = fitIDQNoiseGain(nominalFitTable, initialPsych.sessionFits, ...
    alignedWeibull, targetPerformance);
correctedGains = fitIDQNoiseGain(correctedFitTable, initialPsych.sessionFits, ...
    alignedWeibull, targetPerformance);

pairedGains = makePairedGainTable(nominalGains, correctedGains);
predictorDiagnostics = makePredictorDiagnostics(trialTable);
residualUncertainty = estimateResidualUncertainty(sessionDiagnostics, ...
    correctedLOO, trialTable);

% Retain a compact trial table sufficient to inspect or refit the results.
keepNames = { ...
    'sessionIndex', 'sessionName', 'trialIndex', 'validIdx', 'correct', ...
    'stepCoh', 'hasStepNoise', 'dirIndex', ...
    'nominalNoisePredDir1', 'nominalNoisePredDir2', 'nominalNoisePredDir3', ...
    'correctedNoisePredDir1', 'correctedNoisePredDir2', 'correctedNoisePredDir3'};
compactTrialTable = trialTable(:, keepNames);

correction = struct();
correction.createdAt = datetime('now');
correction.createdBy = mfilename;
correction.methodVersion = 1;
correction.fullSessionFolder = fullSessionFolder;
correction.nSessions = numel(sessionData);
correction.nTrials = height(trialTable);
correction.nStepNoiseTrials = sum(trialTable.hasStepNoise);
correction.predictorType = 'looStep';
correction.initialPsych = initialPsych;
correction.alignedWeibull = alignedWeibull;
correction.nominalGains = nominalGains;
correction.correctedGains = correctedGains;
correction.pairedGains = pairedGains;
correction.predictorDiagnostics = predictorDiagnostics;
correction.residualUncertainty = residualUncertainty;
correction.nominalWeightDiagnostics = nominalLOO.weightDiagnostics;
correction.correctedWeightDiagnostics = correctedLOO.weightDiagnostics;
correction.sessionDiagnostics = struct2table(sessionDiagnostics);
correction.trialTable = compactTrialTable;
correction.notes = { ...
    'Corrected predictors are conditional expected realized signed coherences.'; ...
    'LOO temporal kernels were recomputed from the corrected framewise streams.'; ...
    'Individual historical dot assignments are not reconstructed.'; ...
    'Residual uncertainty quantifies unobserved finite-dot allocation approximately.'};

matPath = fullfile(outputFolder, 'IDQ_RebalancingCorrection.mat');
gainCSVPath = fullfile(outputFolder, 'IDQ_RebalancingCorrection_Gains.csv');
predictorCSVPath = fullfile(outputFolder, 'IDQ_RebalancingCorrection_Predictors.csv');
save(matPath, 'correction', '-v7.3');
writetable(pairedGains, gainCSVPath);
writetable(predictorDiagnostics, predictorCSVPath);

if opts.Verbose
    fprintf('\nPaired gain changes (corrected - nominal):\n');
    disp(pairedGains(:, {'model', 'predictor', 'nominalGain', ...
        'correctedGain', 'gainChange'}));
    fprintf('Saved:\n  %s\n  %s\n  %s\n', matPath, gainCSVPath, predictorCSVPath);
end
end

%% Session reconstruction -------------------------------------------------
function [S, diagnostic] = reconstructOneSession(header, trials, sessionIndex, fileName, verbose)

nAll = numel(trials);
eot = nan(nAll, 1);
certify = nan(nAll, 1);
for i = 1:nAll
    eot(i) = scalarData(trials{i}, 'extendedEOT');
    certify(i) = scalarData(trials{i}, 'trialCertify');
end
validIdx = find((eot == 0 | eot == 1) & certify == 0);
t = trials(validIdx);
nTrials = numel(t);

frameRateHz = firstData(header.frameRateHz);
preStepMS = firstData(header.preStepMS);
stepMS = firstData(header.stepMS);
nFrames = round((preStepMS + stepMS) * frameRateHz / 1000);
tMS = (0:nFrames - 1)' * 1000 / frameRateHz;
stepFrames = find(tMS >= preStepMS & tMS < preStepMS + stepMS);

nominal = nan(nFrames, nTrials, 3);
corrected = nan(nFrames, nTrials, 3);

trialIndex = (1:nTrials)';
correct = false(nTrials, 1);
stepCoh = nan(nTrials, 1);
hasStepNoise = false(nTrials, 1);
dirIndex = nan(nTrials, 1);

numDots = nan(nTrials, 1);
lifeFrames = nan(nTrials, 1);
noiseAmplitude = nan(nTrials, 1);

for iTrial = 1:nTrials
    tr = t{iTrial};
    td = tr.trial.data;
    correct(iTrial) = eot(validIdx(iTrial)) == 0;
    stepCoh(iTrial) = double(td.stepCohPC);
    hasStepNoise(iTrial) = logical(td.stepHasNoise);
    dirIndex(iTrial) = double(td.dirIndex) + 1;

    nominal(:, iTrial, :) = fillNominalFrames(tr, 'change', nFrames, frameRateHz);
    [corrected(:, iTrial, :), oneDiag] = ...
        expectedRealizedFrames(tr, 'change', nFrames, frameRateHz);
    numDots(iTrial) = oneDiag.numDots;
    lifeFrames(iTrial) = oneDiag.lifeFrames;
    noiseAmplitude(iTrial) = oneDiag.noiseAmplitudePC;
    if verbose && (mod(iTrial, 100) == 0 || iTrial == nTrials)
        fprintf('           reconstructed %d/%d valid trials\n', iTrial, nTrials);
    end
end

[~, baseName] = fileparts(fileName);
noisyCoh = unique(stepCoh(hasStepNoise));
if numel(noisyCoh) ~= 1
    error('runIDQRebalancingCorrection:NoisyCohMismatch', ...
        'Expected one noisy step coherence in %s.', fileName);
end

trialTable = table();
trialTable.sessionIndex = repmat(sessionIndex, nTrials, 1);
trialTable.sessionName = repmat(string(baseName), nTrials, 1);
trialTable.trialIndex = trialIndex;
trialTable.validIdx = validIdx(:);
trialTable.correct = correct;
trialTable.stepCoh = stepCoh;
trialTable.hasStepNoise = hasStepNoise;
trialTable.dirIndex = dirIndex;
trialTable.noisyStepCoh = repmat(noisyCoh, nTrials, 1);

S = struct();
S.fileName = baseName;
S.trialTable = trialTable;
S.tMS = tMS;
S.stepFrames = stepFrames;
S.nominalDirNoise = nominal;
S.correctedDirNoise = corrected;
S.nominalSumNoise = sum(nominal, 3);
S.correctedSumNoise = sum(corrected, 3);

diagnostic = emptySessionDiagnostic();
diagnostic.fileName = string(baseName);
diagnostic.nTrials = nTrials;
diagnostic.nStepNoiseTrials = sum(hasStepNoise);
diagnostic.frameRateHz = frameRateHz;
diagnostic.preStepMS = preStepMS;
diagnostic.stepMS = stepMS;
diagnostic.numDotsMin = min(numDots);
diagnostic.numDotsMax = max(numDots);
diagnostic.lifeFramesMin = min(lifeFrames);
diagnostic.lifeFramesMax = max(lifeFrames);
diagnostic.noiseAmplitudePCMin = min(noiseAmplitude);
diagnostic.noiseAmplitudePCMax = max(noiseAmplitude);
end

function d = emptySessionDiagnostic()
d = struct('fileName', "", 'nTrials', NaN, 'nStepNoiseTrials', NaN, ...
    'frameRateHz', NaN, 'preStepMS', NaN, 'stepMS', NaN, ...
    'numDotsMin', NaN, 'numDotsMax', NaN, ...
    'lifeFramesMin', NaN, 'lifeFramesMax', NaN, ...
    'noiseAmplitudePCMin', NaN, 'noiseAmplitudePCMax', NaN);
end

%% Nominal frame expansion ------------------------------------------------
function noise = fillNominalFrames(tr, fieldName, nFrames, frameRateHz)

timesMS = vectorData(tr, [fieldName 'TimesMS']);
coh = nan(numel(timesMS), 3);
for iDir = 1:3
    coh(:, iDir) = vectorData(tr, sprintf('%sCohsPC%d', fieldName, iDir - 1));
end
nSteps = min(size(coh, 1), numel(timesMS));
timesMS = timesMS(1:nSteps);
coh = coh(1:nSteps, :);

noise = nan(nFrames, 1, 3);
msPerFrame = 1000 / frameRateHz;
for iStep = 1:nSteps
    firstFrame = max(1, floor(timesMS(iStep) / msPerFrame) + 1);
    if firstFrame > nFrames
        continue
    end
    if iStep < nSteps
        nextFrame = floor(timesMS(iStep + 1) / msPerFrame) + 1;
    else
        nextFrame = nFrames + 1;
    end
    if nextFrame <= firstFrame
        continue
    end
    nextFrame = min(nextFrame, nFrames + 1);
    for iDir = 1:3
        noise(firstFrame:nextFrame - 1, 1, iDir) = coh(iStep, iDir);
    end
end
end

%% Expected historical movie ---------------------------------------------
function [realizedNoise, diagnostic] = expectedRealizedFrames( ...
    tr, fieldName, nFrames, frameRateHz)

td = tr.trial.data;
dots = tr.([fieldName 'Dots']).data;
timesMS = double(vectorData(tr, [fieldName 'TimesMS']));
coh = nan(numel(timesMS), 3);
for iDir = 1:3
    coh(:, iDir) = double(vectorData(tr, ...
        sprintf('%sCohsPC%d', fieldName, iDir - 1)));
end
nSteps = min(size(coh, 1), numel(timesMS));
timesMS = timesMS(1:nSteps);
coh = coh(1:nSteps, :);

baseDirDeg = double(td.baseDirDeg);
dirsDeg = mod(baseDirDeg + [0 120 240], 360);
dirBins = dirsDeg + 1;
noiseAmplitudePC = abs(double(td.cohNoisePC));
preStepMS = double(td.preStepMS);
stepCohPC = double(td.stepCohPC);
stepDir = double(td.dirIndex) + 1;

fieldWidthDeg = 2 * (double(dots.radiusDeg) + double(dots.speedDPS) / frameRateHz);
numDots = fix(fieldWidthDeg * fieldWidthDeg * double(dots.density));
lifeFrames = double(dots.lifeFrames);
if numDots < 1
    error('runIDQRebalancingCorrection:BadDotCount', 'Computed fewer than one dot.');
end

targets = repmat(struct('probCat', [], 'genCat', [], 'multiplicity', [], ...
    'noisePedestalPC', [], 'stepIncrementPC', [], 'noiseCohPC', []), nSteps, 1);
for iStep = 1:nSteps
    noiseInc = zeros(1, 3);
    nonzero = abs(coh(iStep, :)) > 1e-7;
    noiseInc(nonzero) = (noiseAmplitudePC + coh(iStep, nonzero)) / 100;

    stepInc = zeros(1, 3);
    if strcmp(fieldName, 'change') && timesMS(iStep) >= preStepMS
        stepInc(stepDir) = stepCohPC / 100;
    end
    coherentInc = noiseInc + stepInc;
    bkgProb = (1 - sum(coherentInc)) / 360;
    if bkgProb < 0
        error('runIDQRebalancingCorrection:NegativeBackground', ...
            'Coherent increments exceed total probability.');
    end

    prob360 = repmat(single(bkgProb), 1, 360);
    for iDir = 1:3
        prob360(dirBins(iDir)) = single(bkgProb + coherentInc(iDir));
    end
    % Match LLDotFuncStep's float accumulation and forced final CDF value.
    cdf = zeros(size(prob360), 'single');
    probSum = single(0);
    for iBin = 1:numel(prob360)
      probSum = probSum + prob360(iBin);
      cdf(iBin) = probSum;
    end
    cdf(end) = single(1);
    % cdf = cumsum(prob360, 'native');
    % cdf(end) = single(1);
    edgeCount = ceil(1000 * double(cdf));
    gen360 = diff([0 edgeCount]) / 1000;

    otherMask = true(1, 360);
    otherMask(dirBins) = false;
    targets(iStep).probCat = double([prob360(dirBins), prob360(find(otherMask, 1))]);
    targets(iStep).genCat = [gen360(dirBins), sum(gen360(otherMask))];
    targets(iStep).multiplicity = [1 1 1 357];
    if any(nonzero)
        targets(iStep).noisePedestalPC = noiseAmplitudePC;
    else
        targets(iStep).noisePedestalPC = 0;
    end
    targets(iStep).stepIncrementPC = 100 * stepInc;
    targets(iStep).noiseCohPC = coh(iStep, :);
end

stepFirstFrames = fix(timesMS * frameRateHz / 1000);
currentStep = 1;

if lifeFrames > 0
    age = mod(0:numDots - 1, lifeFrames);
    ageCounts = accumarray(age(:) + 1, 1, [lifeFrames 1]);
    state = ageCounts * targets(1).genCat;
else
    ageCounts = numDots;
    state = numDots * targets(1).genCat;
end

realizedNoise = nan(nFrames, 1, 3);
for iFrame = 1:nFrames
    realizedNoise(iFrame, 1, :) = noiseEvidenceFromState( ...
        state, targets(currentStep), numDots);
    if iFrame == nFrames
        continue
    end

    if lifeFrames > 0
        state = renewExpiredCohort(state, targets(currentStep).genCat);
    end

    % Historical code advanced at most one due function step per frame.
    nextDisplayedFrame0 = iFrame; % zero-based frame following MATLAB frame iFrame
    if currentStep < nSteps && ...
            nextDisplayedFrame0 >= stepFirstFrames(currentStep + 1)
        currentStep = currentStep + 1;
        state = rebalanceExpectedState(state, targets(currentStep), ...
            numDots, dirsDeg);
    end
end

diagnostic = struct('numDots', numDots, 'lifeFrames', lifeFrames, ...
    'noiseAmplitudePC', noiseAmplitudePC);
end

function state = renewExpiredCohort(state, genCat)

nAges = size(state, 1);
expiringCount = sum(state(end, :));
newState = zeros(size(state));
if nAges > 1
    newState(2:end, :) = state(1:end - 1, :);
end
newState(1, :) = expiringCount * genCat;
state = newState;
end

function evidence = noiseEvidenceFromState(state, target, numDots)

countCat = sum(state, 1);
minProb = min(target.probCat);
bkgMask = target.probCat == minProb;
nBkgDirs = sum(target.multiplicity(bkgMask));
bkgCountPerDir = sum(countCat(bkgMask)) / nBkgDirs;

realizedCohPC = 100 * (countCat(1:3) - bkgCountPerDir) / numDots;
evidence = realizedCohPC - target.stepIncrementPC - target.noisePedestalPC;
evidence = reshape(evidence, 1, 1, 3);
end

%% Historical probability-group rebalancing ------------------------------
function state = rebalanceExpectedState(state, target, numDots, dirsDeg)

[groupCats, groupProb, groupNumDirs] = makeProbabilityGroups(target);
nGroups = numel(groupCats);

% LLRDKProbData redistributes a direction whose occupancy is more than ten
% times the target expectation. This usually catches a coherent direction
% that has just joined the background group.
[~, dirOrder] = sort(dirsDeg);
for iOrder = 1:3
    iCat = dirOrder(iOrder);
    iGroup = find(cellfun(@(x) any(x == iCat), groupCats), 1);
    threshold = 10 * numDots * groupProb(iGroup);
    if sum(state(:, iCat)) > threshold
        byAge = state(:, iCat);
        state(:, iCat) = 0;
        weights = target.multiplicity(groupCats{iGroup});
        weights = weights / sum(weights);
        state(:, groupCats{iGroup}) = state(:, groupCats{iGroup}) + ...
            byAge * weights;
    end
end

numDotsGroup = zeros(nGroups, 1);
numToMove = zeros(nGroups, 1);
for iGroup = 1:nGroups
    numDotsGroup(iGroup) = round(sum(state(:, groupCats{iGroup}), 'all'));
    targetCount = round(numDots * groupProb(iGroup) * groupNumDirs(iGroup));
    numToMove(iGroup) = numDotsGroup(iGroup) - targetCount;
end

% Sort by descending excess, matching NSSortDescriptor in the old code.
originalOrder = (1:nGroups)';
[~, sortOrder] = sortrows([-numToMove originalOrder], [1 2]);
groupCats = groupCats(sortOrder);
groupProb = groupProb(sortOrder);
groupNumDirs = groupNumDirs(sortOrder);
numToMove = numToMove(sortOrder);

% Match the integer-rounding correction applied to numToMove.
sumExcess = sum(numToMove);
guard = 0;
while sumExcess ~= 0
    guard = guard + 1;
    if guard > 1000
        error('runIDQRebalancingCorrection:RoundingCorrectionFailed', ...
            'Could not balance integer rounding residual.');
    end
    signExcess = sign(sumExcess);
    changed = false;
    for iGroup = 1:nGroups
        if sumExcess == 0
            break
        end
        if groupProb(iGroup) == 0 || groupProb(iGroup) == 1
            continue
        end
        if sign(numToMove(iGroup)) == signExcess
            while sumExcess ~= 0 && numToMove(iGroup) ~= 0
                sumExcess = sumExcess - signExcess;
                numToMove(iGroup) = numToMove(iGroup) - signExcess;
                changed = true;
            end
        end
    end
    if ~changed
        error('runIDQRebalancingCorrection:RoundingCorrectionStalled', ...
            'Integer rounding correction stalled.');
    end
end

% Move dots with the original do-while short-circuit semantics. In
% particular, the source counter is not decremented when incrementing the
% target counter makes the first condition false.
for iSource = 1:nGroups - 1
    if numToMove(iSource) <= 0
        continue
    end
    for iTarget = iSource + 1:nGroups
        if numToMove(iTarget) >= 0
            continue
        end
        while true
            state = transferExpectedDot(state, groupCats{iSource}, ...
                groupCats{iTarget}, target.multiplicity);
            numToMove(iTarget) = numToMove(iTarget) + 1;
            if numToMove(iTarget) < 0
                numToMove(iSource) = numToMove(iSource) - 1;
                if numToMove(iSource) > 0
                    continue
                end
            end
            break
        end
    end
end
end

function state = transferExpectedDot(state, sourceCats, targetCats, multiplicity)

sourceMass = sum(state(:, sourceCats), 'all');
if sourceMass <= 0
    return
end
removed = state(:, sourceCats) / sourceMass;
state(:, sourceCats) = state(:, sourceCats) - removed;
removedByAge = sum(removed, 2);
targetWeights = multiplicity(targetCats);
targetWeights = targetWeights / sum(targetWeights);
state(:, targetCats) = state(:, targetCats) + removedByAge * targetWeights;
end

function [groupCats, groupProb, groupNumDirs] = makeProbabilityGroups(target)

probCat = target.probCat;
[uniqueProb, ~, groupIndex] = unique(probCat, 'sorted');
nGroups = numel(uniqueProb);
groupCats = cell(nGroups, 1);
groupProb = nan(nGroups, 1);
groupNumDirs = nan(nGroups, 1);
for iGroup = 1:nGroups
    groupCats{iGroup} = find(groupIndex == iGroup);
    groupProb(iGroup) = uniqueProb(iGroup);
    groupNumDirs(iGroup) = sum(target.multiplicity(groupCats{iGroup}));
end
end

%% Leave-one-session-out predictors --------------------------------------
function out = computeLOOPredictors(sessionData, kind)

nSessions = numel(sessionData);
nFrames = numel(sessionData{1}.tMS);
predCell = cell(nSessions, 3);
weightsBySession = nan(nFrames, nSessions);

sessionIndex = (1:nSessions)';
fileName = strings(nSessions, 1);
kernelStepSum = nan(nSessions, 1);
kernelStepMean = nan(nSessions, 1);
kernelStepMin = nan(nSessions, 1);
kernelStepMax = nan(nSessions, 1);
fracNegativeWeights = nan(nSessions, 1);
effectiveNFrames = nan(nSessions, 1);

for iSession = 1:nSessions
    SA = sessionData{iSession};
    useFrames = SA.stepFrames;
    allCorrect = [];
    allError = [];
    for jSession = 1:nSessions
        if jSession == iSession && nSessions > 1
            continue
        end
        SB = sessionData{jSession};
        T = SB.trialTable;
        idxCorrect = T.hasStepNoise & T.correct;
        idxError = T.hasStepNoise & ~T.correct;
        sumNoise = getNoiseArray(SB, kind, true);
        allCorrect = [allCorrect sumNoise(:, idxCorrect)]; %#ok<AGROW>
        allError = [allError sumNoise(:, idxError)]; %#ok<AGROW>
    end
    kernel = mean(allCorrect, 2, 'omitnan') - mean(allError, 2, 'omitnan');
    kUse = kernel(useFrames);
    kSum = sum(kUse, 'omitnan');
    if ~isfinite(kSum) || kSum == 0
        error('runIDQRebalancingCorrection:BadKernel', ...
            'LOO kernel is zero or nonfinite for %s.', SA.fileName);
    end
    w = kUse(:) / kSum;
    dirNoise = getNoiseArray(SA, kind, false);
    for iDir = 1:3
        predCell{iSession, iDir} = ...
            sum(dirNoise(useFrames, :, iDir) .* w, 1, 'omitnan')';
    end
    weightsBySession(useFrames, iSession) = w;
    fileName(iSession) = string(SA.fileName);
    kernelStepSum(iSession) = kSum;
    kernelStepMean(iSession) = mean(kUse, 'omitnan');
    kernelStepMin(iSession) = min(kUse, [], 'omitnan');
    kernelStepMax(iSession) = max(kUse, [], 'omitnan');
    fracNegativeWeights(iSession) = mean(w < 0);
    effectiveNFrames(iSession) = 1 / sum(w .^ 2);
end

out = struct();
out.noisePredDir1 = vertcat(predCell{:, 1});
out.noisePredDir2 = vertcat(predCell{:, 2});
out.noisePredDir3 = vertcat(predCell{:, 3});
out.weightsBySession = weightsBySession;
out.weightDiagnostics = table(sessionIndex, fileName, kernelStepSum, ...
    kernelStepMean, kernelStepMin, kernelStepMax, fracNegativeWeights, ...
    effectiveNFrames);
end

function x = getNoiseArray(S, kind, summed)
switch kind
    case 'nominal'
        if summed
            x = S.nominalSumNoise;
        else
            x = S.nominalDirNoise;
        end
    case 'corrected'
        if summed
            x = S.correctedSumNoise;
        else
            x = S.correctedDirNoise;
        end
    otherwise
        error('Unknown predictor kind: %s', kind);
end
end

function T = selectPredictors(T, prefix)
for iDir = 1:3
    T.(sprintf('noisePredDir%d', iDir)) = ...
        T.(sprintf('%sNoisePredDir%d', prefix, iDir));
end
end

%% Compact results --------------------------------------------------------
function paired = makePairedGainTable(nominal, corrected)

models = {'combined', 'absolute', 'driftRelative', 'driftNonDrift'};
model = strings(0, 1);
predictor = strings(0, 1);
nominalGain = [];
nominalSE = [];
nominalCILow = [];
nominalCIHigh = [];
correctedGain = [];
correctedSE = [];
correctedCILow = [];
correctedCIHigh = [];
for iModel = 1:numel(models)
    a = nominal.(models{iModel});
    b = corrected.(models{iModel});
    for iGain = 1:numel(a.gain)
        model(end + 1, 1) = string(models{iModel}); %#ok<AGROW>
        predictor(end + 1, 1) = string(a.predictorNames{iGain}); %#ok<AGROW>
        nominalGain(end + 1, 1) = a.gain(iGain); %#ok<AGROW>
        nominalSE(end + 1, 1) = a.SE(iGain); %#ok<AGROW>
        nominalCILow(end + 1, 1) = a.CI95(iGain, 1); %#ok<AGROW>
        nominalCIHigh(end + 1, 1) = a.CI95(iGain, 2); %#ok<AGROW>
        correctedGain(end + 1, 1) = b.gain(iGain); %#ok<AGROW>
        correctedSE(end + 1, 1) = b.SE(iGain); %#ok<AGROW>
        correctedCILow(end + 1, 1) = b.CI95(iGain, 1); %#ok<AGROW>
        correctedCIHigh(end + 1, 1) = b.CI95(iGain, 2); %#ok<AGROW>
    end
end
gainChange = correctedGain - nominalGain;
percentChange = 100 * gainChange ./ abs(nominalGain);
paired = table(model, predictor, nominalGain, nominalSE, nominalCILow, ...
    nominalCIHigh, correctedGain, correctedSE, correctedCILow, ...
    correctedCIHigh, gainChange, percentChange);
end

function diagnostics = makePredictorDiagnostics(T)

idx = T.hasStepNoise;
Xn = [T.nominalNoisePredDir1 T.nominalNoisePredDir2 T.nominalNoisePredDir3];
Xc = [T.correctedNoisePredDir1 T.correctedNoisePredDir2 T.correctedNoisePredDir3];

label = strings(0, 1);
n = [];
correlation = [];
intercept = [];
calibrationSlope = [];
throughOriginSlope = [];
rmse = [];
meanCorrection = [];
sdCorrection = [];

for iDir = 1:3
    addOne(sprintf('Absolute direction %d', iDir), Xn(idx, iDir), Xc(idx, iDir));
end

row = (1:height(T))';
drift = T.dirIndex;
plusDir = mod(drift, 3) + 1;
minusDir = mod(drift - 2, 3) + 1;
xnDrift = Xn(sub2ind(size(Xn), row, drift));
xcDrift = Xc(sub2ind(size(Xc), row, drift));
xnPlus = Xn(sub2ind(size(Xn), row, plusDir));
xcPlus = Xc(sub2ind(size(Xc), row, plusDir));
xnMinus = Xn(sub2ind(size(Xn), row, minusDir));
xcMinus = Xc(sub2ind(size(Xc), row, minusDir));
addOne('Drift direction', xnDrift(idx), xcDrift(idx));
addOne('+120 degrees', xnPlus(idx), xcPlus(idx));
addOne('-120 degrees', xnMinus(idx), xcMinus(idx));
addOne('Pooled non-drift', [xnPlus(idx); xnMinus(idx)], ...
    [xcPlus(idx); xcMinus(idx)]);
addOne('Sum of directions', sum(Xn(idx, :), 2), sum(Xc(idx, :), 2));

diagnostics = table(label, n, correlation, intercept, calibrationSlope, ...
    throughOriginSlope, rmse, meanCorrection, sdCorrection);

    function addOne(theLabel, x, y)
        valid = isfinite(x) & isfinite(y);
        x = double(x(valid));
        y = double(y(valid));
        label(end + 1, 1) = string(theLabel);
        n(end + 1, 1) = numel(x);
        if numel(x) > 1
            R = corrcoef(x, y);
            correlation(end + 1, 1) = R(1, 2);
            b = [ones(size(x)) x] \ y;
            intercept(end + 1, 1) = b(1);
            calibrationSlope(end + 1, 1) = b(2);
            throughOriginSlope(end + 1, 1) = (x' * y) / (x' * x);
            rmse(end + 1, 1) = sqrt(mean((y - x).^2));
            meanCorrection(end + 1, 1) = mean(y - x);
            sdCorrection(end + 1, 1) = std(y - x);
        else
            correlation(end + 1, 1) = NaN;
            intercept(end + 1, 1) = NaN;
            calibrationSlope(end + 1, 1) = NaN;
            throughOriginSlope(end + 1, 1) = NaN;
            rmse(end + 1, 1) = NaN;
            meanCorrection(end + 1, 1) = NaN;
            sdCorrection(end + 1, 1) = NaN;
        end
    end
end

function residual = estimateResidualUncertainty(sessionDiagnostics, loo, T)

N = median([sessionDiagnostics.numDotsMin sessionDiagnostics.numDotsMax], 'all');
L = median([sessionDiagnostics.lifeFramesMin sessionDiagnostics.lifeFramesMax], 'all');
A = median([sessionDiagnostics.noiseAmplitudePCMin ...
    sessionDiagnostics.noiseAmplitudePCMax], 'all');
movedDots = N * (2 * A / 100);
onsetSDPC = 100 * sqrt(movedDots * 0.25) / N;

% Approximate decay weights for an isolated redistribution error.
decay = (L:-1:1) / L;
rectangularMeanFactor = sum(decay) / 6; % six 120-Hz frames per 50-ms sample
fiveSampleMeanSDPC = onsetSDPC * rectangularMeanFactor / sqrt(5);

% Scale the reference estimate by each LOO temporal weight vector. This is
% deliberately a simple estimate, not a full covariance propagation.
w = loo.weightsBySession;
weightedScale = nan(size(w, 2), 1);
for iSession = 1:size(w, 2)
    wi = w(:, iSession);
    wi = wi(isfinite(wi));
    weightedScale(iSession) = sqrt(sum(wi.^2)) / sqrt(1 / 30);
end

residual = struct();
residual.description = [ ...
    'Approximate SD from the unknown allocation of a two-way redistributed ' ...
    'high-coherence dot pool; it excludes uncertainty in fitted gains.'];
residual.numDots = N;
residual.lifeFrames = L;
residual.noiseAmplitudePC = A;
residual.approxMovedDots = movedDots;
residual.onsetAllocationSDPC = onsetSDPC;
residual.rectangularFiveNoiseSampleMeanSDPC = fiveSampleMeanSDPC;
residual.looWeightedScaleBySession = weightedScale;
residual.approxLOOPredictorSDPCBySession = fiveSampleMeanSDPC * weightedScale;
residual.nStepNoiseTrials = sum(T.hasStepNoise);
residual.caveat = [ ...
    'Aperture sampling and dependence across adjacent transitions can increase ' ...
    'this value; use roughly 0.2-0.4 percentage points as a conservative range.'];
end

%% Raw Lablib helpers -----------------------------------------------------
function x = scalarData(s, fieldName)
x = firstData(s.(fieldName));
end

function x = firstData(x)
while isstruct(x) && isfield(x, 'data')
    x = x.data;
end
x = double(x(1));
end

function x = vectorData(s, fieldName)
x = s.(fieldName);
while isstruct(x) && isfield(x, 'data')
    x = x.data;
end
x = double(x(:));
end
