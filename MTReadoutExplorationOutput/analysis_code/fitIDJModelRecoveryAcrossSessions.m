function [result, fig] = fitIDJModelRecoveryAcrossSessions(varargin)
% fitIDJModelRecoveryAcrossSessions  Leave-one-session-out IDJ model recovery.
%
% [result, fig] = fitIDJModelRecoveryAcrossSessions(Name, Value, ...)
%
% SessionAnalysis files are selected from IDJ/Data/SessionAnalysis using
% selectAnalysisFiles. Name-value arguments not recognized as analysis
% options below are passed through to selectAnalysisFiles, for example:
%
%   fitIDJModelRecoveryAcrossSessions('Animal', 'Meetz')
%   fitIDJModelRecoveryAcrossSessions('Animal', 'All')
%   fitIDJModelRecoveryAcrossSessions('Animal','All','SigmaFitMode','reference')
%   fitIDJModelRecoveryAcrossSessions('Animal','All','IncludeNoiseFreeSessions',true)
%
% Primary analysis excludes sessions with no recorded step noise. Set
% IncludeNoiseFreeSessions=true to include them as a secondary analysis.
% Exclusions are retained in result.excludedNoiseFreeSessions.
%
% Each leave-one-session-out kernel uses noisy 100:0:0 trials from OTHER
% sessions of the same animal. Sessions without informative training noise
% use a uniform 250-ms kernel. Zero-noise trials still enter behavioral fits.
% SigmaFitMode='all' (default) fits one sigma per session/model using ALL
% allocations and scores those same trials: these are fitted likelihoods.
% SigmaFitMode='reference' fits sigma on 100:0:0 and scores divided allocations
% as held-out predictions. Kernel training remains LOO in both modes.
% Three figures are exported when both animals are present; the returned
% figure handle is the combined figure.
%
% The primary analysis:
%   1. estimates one temporal kernel from other sessions' 100:0:0 trials;
%   2. normalizes that kernel to unit L1 mass;
%   3. temporally integrates absolute coherence within each stream;
%   4. pools the three streams with L1, L2, L-infinity, or Square;
%   5. fits one post-pooling Gaussian sigma per session/model, shared across
%      allocations (all trials by default; 100:0 in reference mode); and
%   6. scores all trials (all mode) or divided trials (reference mode).
%
% Square sums squared integrated stream responses WITHOUT the L2 square
% root: D = sum(changeStream.^2,2) - sum(noChangeStream.^2,2).
% The coherence pedestal and recorded noise are retained before squaring;
% no additional baseline subtraction is introduced. Integration precedes
% squaring. Square has its own fitted Gaussian sigma, in squared
% coherence units, and the same parameter count as the other models.
% Model order remains L1, L2, Linf, with Square appended as column four.
% Primary scores are scoredLogLikelihood and totalScoredLogLikelihood.
% Legacy heldOutLogLikelihood/totalHeldOutLogLikelihood are EMPTY in all
% mode; they are populated only in reference mode. noiseAuditTable reports
% actual recorded noise and pedestal ranges, not just the enable flag.
%
% The default correct-minus-error kernel is deliberately simple. Probit and
% logit alternatives use five 50-ms predictors and can be selected with
% 'KernelMethod'. The same leave-one-session-out kernel is used for both
% streams, both stimulus patches, and all allocations in the held-out
% session.

p = inputParser;
p.KeepUnmatched = true;
p.PartialMatching = false;
p.addParameter('KernelMethod', 'correct-error', ...
  @(x) any(strcmpi(string(x), ["correct-error" "probit" "logit"])));
p.addParameter('ReferenceAllocationPC', 100, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('SigmaFitMode', 'all', ...
  @(x) any(strcmpi(string(x), ["all" "reference"])));
p.addParameter('IncludeNoiseFreeSessions', false, ...
  @(x) islogical(x) && isscalar(x));
% p.addParameter('PDFPath', '', @(x) ischar(x) || isstring(x));
% p.addParameter('MATPath', '', @(x) ischar(x) || isstring(x));
p.addParameter('ExportPDF', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MakePlots', true, @(x) islogical(x) && isscalar(x));
p.addParameter('SaveResult', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Visible', 'on', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opts = p.Results;
opts.KernelMethod = lower(string(opts.KernelMethod));
opts.SigmaFitMode = lower(string(opts.SigmaFitMode));

[dataFolder, existed] = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'SessionAnalysis'));
assert(existed, 'fitIDJModelRecoveryAcrossSessions:MissingDataFolder', ...
  'The IDJ SessionAnalysis folder does not exist: %s', dataFolder);
selectionArguments = unmatchedNameValuePairs(p.Unmatched);
sourceList = selectAnalysisFiles(dataFolder, selectionArguments{:});
if ischar(sourceList)
  sourceList = {sourceList};
elseif isstring(sourceList)
  sourceList = cellstr(sourceList(:));
else
  sourceList = sourceList(:);
end
assert(iscell(sourceList), ...
  'fitIDJModelRecoveryAcrossSessions:UnexpectedFileList', ...
  'selectAnalysisFiles must return file names as a cell or string array.');
nSessions = numel(sourceList);
assert(nSessions >= 1, 'fitIDJModelRecoveryAcrossSessions:TooFewSessions', ...
  'At least one session is required.');

A = cell(nSessions, 1);
for k = 1:nSessions
  A{k} = loadSessionAnalysis(sourceList{k});
end
validateSessions(A, opts.ReferenceAllocationPC);
allIDs = string(cellfun(@(x) char(x.sourceSessionID), A, 'UniformOutput',false));
allAnimals = string(cellfun(@(x) char(x.sessionHeader.animal), A, 'UniformOutput',false));
audit = makeNoiseAuditTable(A,allIDs(:),allAnimals(:));
selected = opts.IncludeNoiseFreeSessions | audit.nNoiseTrials>0;
excludedNoiseFreeSessions = audit(~selected,:);
A = A(selected);
sourceList = sourceList(selected);
audit = audit(selected,:);
nSessions = numel(A);
assert(nSessions>=1, 'fitIDJModelRecoveryAcrossSessions:NoNoisySessions', ...
  ['No sessions with recorded step noise were selected. To analyze noise-free ' ...
   'sessions, set IncludeNoiseFreeSessions=true.']);
animal = string(cellfun(@(x) char(x.sessionHeader.animal), A, ...
  'UniformOutput', false));
animal = animal(:);

modelNames = ["L1" "L2" "Linf" "Square"];
nModels = numel(modelNames);
nFrames = numel(A{1}.stepTimeMS);
kernelByFrameSession = nan(nFrames, nSessions);
kernelTrainingSessionCount = zeros(nSessions, 1);
kernelSourceBySession = strings(nSessions, 1);
sigmaBySessionModel = nan(nSessions, nModels);
fitLogLikelihood = nan(nSessions, nModels);
referenceLogLikelihood = nan(nSessions, nModels);
dividedLogLikelihood = nan(nSessions, nModels);
scoredLogLikelihood = nan(nSessions, nModels);
referenceTrialCount = zeros(nSessions, 1);
dividedTrialCount = zeros(nSessions, 1);
fitTrialCount = zeros(nSessions, 1);
scoredTrialCount = zeros(nSessions, 1);
predictionTables = cell(nSessions, 1);

kernelSummaries = makeKernelSummaries(A, opts.ReferenceAllocationPC);
for heldOut = 1:nSessions
  included = find(animal == animal(heldOut) & [kernelSummaries.canTrain]');
  included(included == heldOut) = [];
  kernelTrainingSessionCount(heldOut) = numel(included);
  if isempty(included)
    kernel = ones(nFrames, 1) / nFrames;
    kernelSourceBySession(heldOut) = "uniform boxcar (no informative training noise)";
  else
    [kernel, informativeKernel] = estimateLOOKernel(A, kernelSummaries, included, ...
      opts.ReferenceAllocationPC, opts.KernelMethod);
    kernelSourceBySession(heldOut) = "same-animal LOO";
    if ~informativeKernel
      kernelSourceBySession(heldOut) = "uniform boxcar (zero training kernel)";
    end
  end
  kernelByFrameSession(:,heldOut) = kernel;

  T = A{heldOut}.trialTable;
  referenceMask = abs(double(T.unsignedAllocPC) - ...
    opts.ReferenceAllocationPC) < 1e-6;
  dividedMask = ~referenceMask;
  if opts.SigmaFitMode == "all"
    fitMask = true(height(T),1);
    scoreMask = fitMask;
  else
    fitMask = referenceMask;
    scoreMask = dividedMask;
  end
  referenceTrialCount(heldOut) = sum(referenceMask);
  dividedTrialCount(heldOut) = sum(dividedMask);
  fitTrialCount(heldOut) = sum(fitMask);
  scoredTrialCount(heldOut) = sum(scoreMask);

  probabilities = nan(height(T), nModels);
  for modelIndex = 1:nModels
    decisionEvidence = modelEvidence(A{heldOut}, kernel, modelNames(modelIndex));
    sigma = fitDecisionSigma(decisionEvidence(fitMask), T.correct(fitMask));
    probability = gaussianCDF(decisionEvidence ./ sigma);

    sigmaBySessionModel(heldOut,modelIndex) = sigma;
    probabilities(:,modelIndex) = probability;
    fitLogLikelihood(heldOut,modelIndex) = bernoulliLogLikelihood( ...
      probability(fitMask), T.correct(fitMask));
    referenceLogLikelihood(heldOut,modelIndex) = bernoulliLogLikelihood( ...
      probability(referenceMask), T.correct(referenceMask));
    dividedLogLikelihood(heldOut,modelIndex) = bernoulliLogLikelihood( ...
      probability(dividedMask), T.correct(dividedMask));
    scoredLogLikelihood(heldOut,modelIndex) = bernoulliLogLikelihood( ...
      probability(scoreMask), T.correct(scoreMask));
  end

  n = height(T);
  sourceSessionID = repmat(string(A{heldOut}.sourceSessionID), n, 1);
  P = table(sourceSessionID, double(T.sourceTrialNumber), logical(T.correct), ...
    double(T.unsignedAllocPC), double(T.totalStepCohPC), fitMask, ...
    probabilities(:,1), probabilities(:,2), probabilities(:,3), ...
    probabilities(:,4), ...
    'VariableNames', {'sourceSessionID','sourceTrialNumber','correct', ...
    'unsignedAllocPC','totalStepCohPC','isCalibrationTrial', ...
    'pL1','pL2','pLinf','pSquare'});
  P.alloc0 = double(T.alloc0);
  P.alloc1 = double(T.alloc1);
  P.alloc2 = double(T.alloc2);
  P.isReferenceTrial = referenceMask;
  P.isFitTrial = fitMask;
  P.isScoredTrial = scoreMask;
  P.prestepCohPC = double(T.prestepCohPC);
  P.hasRecordedStepNoise = hasStepNoise(A{heldOut});
  P.noiseCohPC = nan(n,1);
  if ismember('noiseCohPC',T.Properties.VariableNames)
    P.noiseCohPC = double(T.noiseCohPC);
  end
  predictionTables{heldOut} = P;
end

predictionTable = vertcat(predictionTables{:});
psychometricTable = makePsychometricTable(predictionTable);
sessionIDs = string(cellfun(@(x) char(x.sourceSessionID), A, ...
  'UniformOutput', false));
sessionTable = table(sessionIDs(:), animal(:), ...
  kernelSourceBySession, kernelTrainingSessionCount, ...
  referenceTrialCount, dividedTrialCount, fitTrialCount, scoredTrialCount, ...
  'VariableNames', {'sourceSessionID','animal','kernelSource', ...
  'nKernelTrainingSessions','nReferenceTrials','nDividedTrials', ...
  'nFitTrials','nScoredTrials'});
for modelIndex = 1:nModels
  suffix = char(modelNames(modelIndex));
  sessionTable.(['sigma' suffix]) = sigmaBySessionModel(:,modelIndex);
  sessionTable.(['scoredLogLikelihood' suffix]) = scoredLogLikelihood(:,modelIndex);
end
sessionTable.nCalibrationTrials = fitTrialCount;
sessionTable.nHeldOutTrials = zeros(nSessions,1);
if opts.SigmaFitMode == "reference"
  sessionTable.nHeldOutTrials = dividedTrialCount;
  for modelIndex = 1:nModels
    sessionTable.(['heldOutLogLikelihood' char(modelNames(modelIndex))]) = ...
      dividedLogLikelihood(:,modelIndex);
  end
end

result = struct();
result.version = 2;
result.topUpTrialsExcluded = true;
result.sourceDomain = 'IDJ';
result.nSignalStreams = 3;
result.createdBy = mfilename;
result.createdAt = datetime('now');
result.sourceFolder = string(dataFolder);
result.selectionOptions = p.Unmatched;
result.sourceFiles = string(sourceList);
result.kernelMethod = opts.KernelMethod;
result.sigmaFitMode = opts.SigmaFitMode;
result.includeNoiseFreeSessions = opts.IncludeNoiseFreeSessions;
result.excludedNoiseFreeSessions = excludedNoiseFreeSessions;
if opts.SigmaFitMode == "all"
  result.scoreDescription = "Fitted likelihood, all allocations";
else
  result.scoreDescription = "Held-out likelihood, divided allocations";
end
result.kernelContract = ...
  ['Same-animal LOO using noisy reference trials only; boxcar if no ' ...
   'informative training session; unit L1 mass; one kernel for all streams'];
result.modelContract = ...
  ['absolute coherence including pedestal and recorded noise -> temporal ' ...
   'integration within stream -> L1/L2/Linf or sum of squared responses ' ...
   'without root -> Change minus No-change -> Gaussian decision noise'];
result.referenceAllocationPC = opts.ReferenceAllocationPC;
result.modelNames = modelNames;
result.stepTimeMS = double(A{1}.stepTimeMS(:));
result.kernelByFrameSession = kernelByFrameSession;
result.kernelSourceBySession = kernelSourceBySession;
result.kernelTrainingSessionCount = kernelTrainingSessionCount;
result.meanLOOKernel = mean(kernelByFrameSession, 2);
result.sigmaBySessionModel = sigmaBySessionModel;
result.fitLogLikelihood = fitLogLikelihood;
result.calibrationLogLikelihood = fitLogLikelihood;
result.referenceLogLikelihood = referenceLogLikelihood;
result.dividedLogLikelihood = dividedLogLikelihood;
result.scoredLogLikelihood = scoredLogLikelihood;
result.totalScoredLogLikelihood = sum(scoredLogLikelihood, 1);
result.heldOutLogLikelihood = [];
result.totalHeldOutLogLikelihood = [];
if opts.SigmaFitMode == "reference"
  result.heldOutLogLikelihood = dividedLogLikelihood;
  result.totalHeldOutLogLikelihood = sum(dividedLogLikelihood,1);
end
result.deltaLogLikelihoodVsLinf = result.totalScoredLogLikelihood - ...
  result.totalScoredLogLikelihood(3);
result.sessionTable = sessionTable;
result.psychometricTable = psychometricTable;
result.predictionTable = predictionTable;
result.noiseAuditTable = audit;
result.kernelTrainingScope = ...
  'Same-animal LOO, noisy reference trials only; boxcar without informative training sessions';
animalNames = unique(animal, 'stable');
if numel(animalNames) == 1
  result.animal = animalNames(1);
else
  result.animal = "Combined";
end
result.groupTable = makeModelComparisonTable(result);

result.animalResults = struct();
for k = 1:numel(animalNames)
  field = matlab.lang.makeValidName(char(animalNames(k)));
  result.animalResults.(field) = subsetRecoveryResult(result, animalNames(k));
end

fig = gobjects(0);
if opts.MakePlots
  fig = plotRecoverySummary(result, char(opts.Visible));
  [plotFolder] = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Plots', 'AcrossSessions'));
  if numel(animalNames) > 1
    plotPath = fullfile(plotFolder, 'ModelRecovery_Combined.pdf');
  else
    plotPath = fullfile(plotFolder, sprintf('ModelRecovery_%s.pdf', char(animalNames(1))));
  end
  if opts.ExportPDF
    exportgraphics(fig, plotPath, 'ContentType', 'vector');
  end
  if numel(animalNames) > 1
    for k = 1:numel(animalNames)
      field = matlab.lang.makeValidName(char(animalNames(k)));
      animalFig = plotRecoverySummary(result.animalResults.(field), char(opts.Visible));
      animalPath = fullfile(plotFolder, sprintf('ModelRecovery_%s.pdf', field));
      if opts.ExportPDF
        exportgraphics(animalFig, animalPath, 'ContentType', 'vector');
      end
      if strcmpi(opts.Visible,'off'), close(animalFig); end
    end
  end
end % MakePlots
if opts.SaveResult
  [dataFolder] = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'AcrossSessions'));
  dataPath = fullfile(dataFolder, 'ModelRecovery.mat');
  save(char(dataPath), 'result', '-v7.3');
end % SaveResult
end

%% ------------------------------------------------------------------------
function nameValueArgs = unmatchedNameValuePairs(unmatched)

names = fieldnames(unmatched);
nameValueArgs = cell(1, 2 * numel(names));
for k = 1:numel(names)
  nameValueArgs{2 * k - 1} = names{k};
  nameValueArgs{2 * k} = unmatched.(names{k});
end
end

%% ------------------------------------------------------------------------
function A = loadSessionAnalysis(source)

S = load(source);
if isfield(S, 'sessionAnalysis')
  A = S.sessionAnalysis;
  assert(idjHasTopUpExclusion(A),'fitIDJModelRecoveryAcrossSessions:RebuildAnalysis', ...
    'Rebuild compact sessions with makeIDJSessionAnalyses(''Replace'',true); top-up exclusion is required.');
else
  error('fitIDJModelRecoveryAcrossSessions:MissingSessionAnalysis', ...
    '%s does not contain sessionAnalysis.', source);
end
end

%% ------------------------------------------------------------------------
function validateSessions(A, referenceAllocationPC)

stepTimeMS = double(A{1}.stepTimeMS(:));
for k = 1:numel(A)
  required = {'sourceSessionID','trialTable','stepTimeMS', ...
    'changeNoiseByFrameTrial','noChangeNoiseByFrameTrial'};
  assert(all(isfield(A{k}, required)), ...
    'fitIDJModelRecoveryAcrossSessions:MissingFields', ...
    'Session %d lacks one or more required sessionAnalysis fields.', k);
  assert(isequal(double(A{k}.stepTimeMS(:)), stepTimeMS), ...
    'fitIDJModelRecoveryAcrossSessions:TimeGridMismatch', ...
    'All sessions must use the same reconstructed video-frame grid.');
  n = height(A{k}.trialTable);
  requiredTrial = {'alloc0','alloc1','alloc2','highStreamIndex', ...
    'unsignedAllocPC','totalStepCohPC','prestepCohPC','correct','sourceTrialNumber'};
  assert(all(ismember(requiredTrial,A{k}.trialTable.Properties.VariableNames)), ...
    'fitIDJModelRecoveryAcrossSessions:TrialSchema', ...
    'Session %d is missing IDJ trial fields; rebuild its compact analysis.',k);

  assert(isequal(size(A{k}.changeNoiseByFrameTrial), ...
    [numel(stepTimeMS) n 3]), ...
    'fitIDJModelRecoveryAcrossSessions:ChangeNoiseShape', ...
    'Unexpected change-noise array shape in session %d.', k);
  assert(isequal(size(A{k}.noChangeNoiseByFrameTrial), ...
    [numel(stepTimeMS) n 3]), ...
    'fitIDJModelRecoveryAcrossSessions:NoChangeNoiseShape', ...
    'Unexpected no-change-noise array shape in session %d.', k);
  T = A{k}.trialTable;
  assert(all(isfinite(A{k}.changeNoiseByFrameTrial(:))) && ...
    all(isfinite(A{k}.noChangeNoiseByFrameTrial(:))) && ...
    all(isfinite(T.prestepCohPC)) && all(double(T.prestepCohPC)>=0), ...
    'fitIDJModelRecoveryAcrossSessions:NoiseValues', ...
    'Noise arrays and trial pedestals must be finite in session %d.', k);
  if ismember('noiseCohPC',T.Properties.VariableNames)
    declaredZero = abs(double(T.noiseCohPC))<1e-8;
    actualNoise = hasStepNoise(A{k});
    assert(~any(declaredZero & actualNoise), ...
      'fitIDJModelRecoveryAcrossSessions:NoiseMetadataMismatch', ...
      ['Session %d declares zero coherence noise on trials with nonzero ' ...
       'recorded step noise; check the compact-session reconstruction.'], k);
  end
  assert(any(abs(double(A{k}.trialTable.unsignedAllocPC) - ...
    referenceAllocationPC) < 1e-6), ...
    'fitIDJModelRecoveryAcrossSessions:MissingReferenceAllocation', ...
    'Session %d has no %s trials.', k, idjAllocationLabel(referenceAllocationPC));
end
end

%% ------------------------------------------------------------------------
function present = hasStepNoise(A)
% Use recorded amplitudes, not the stepHasNoise enable flag.
C = any(any(abs(double(A.changeNoiseByFrameTrial))>1e-8,1),3);
N = any(any(abs(double(A.noChangeNoiseByFrameTrial))>1e-8,1),3);
present = reshape(C | N,[],1);
end

function audit = makeNoiseAuditTable(A, sessionIDs, animal)
n = numel(A);
nTrials = zeros(n,1); nNoiseTrials = zeros(n,1);
nZeroNoiseTrials = zeros(n,1); minPedestalPC = zeros(n,1);
maxPedestalPC = zeros(n,1); rmsNoisePC = zeros(n,1);
for k = 1:n
  T = A{k}.trialTable;
  nTrials(k) = height(T);
  nNoiseTrials(k) = sum(hasStepNoise(A{k}));
  nZeroNoiseTrials(k) = nTrials(k)-nNoiseTrials(k);
  minPedestalPC(k) = min(double(T.prestepCohPC));
  maxPedestalPC(k) = max(double(T.prestepCohPC));
  C = double(A{k}.changeNoiseByFrameTrial);
  N = double(A{k}.noChangeNoiseByFrameTrial);
  rmsNoisePC(k) = sqrt((sum(C(:).^2)+sum(N(:).^2))/(numel(C)+numel(N)));
end
audit = table(sessionIDs,animal,nTrials,nNoiseTrials,nZeroNoiseTrials, ...
  minPedestalPC,maxPedestalPC,rmsNoisePC, ...
  'VariableNames',{'sourceSessionID','animal','nTrials','nNoiseTrials', ...
  'nZeroNoiseTrials','minPedestalPC','maxPedestalPC','rmsNoisePC'});
end

%% ------------------------------------------------------------------------
function summaries = makeKernelSummaries(A, referenceAllocationPC)

nSessions = numel(A);
nFrames = numel(A{1}.stepTimeMS);
summaries = repmat(struct('meanCorrect', zeros(nFrames,1), ...
  'meanError', zeros(nFrames,1), 'nCorrect', 0, 'nError', 0, ...
  'canTrain', false), nSessions, 1);
for sessionIndex = 1:nSessions
  T = A{sessionIndex}.trialTable;
  mask = abs(double(T.unsignedAllocPC) - referenceAllocationPC) < 1e-6;
  rows = find(mask);
  X = nan(nFrames, numel(rows));
  for j = 1:numel(rows)
    stream = double(T.highStreamIndex(rows(j)));
    assert(ismember(stream, [1 2 3]), ...
      'fitIDJModelRecoveryAcrossSessions:ReferenceStream', ...
      'Every reference-allocation trial must have one high stream.');
    X(:,j) = double(A{sessionIndex}.changeNoiseByFrameTrial(:,rows(j),stream));
  end
  informative = any(abs(X)>1e-8,1);
  X = X(:,informative);
  rows = rows(informative);
  isCorrect = logical(T.correct(rows));
  if ~any(isCorrect) || ~any(~isCorrect)
    continue
  end
  summaries(sessionIndex).meanCorrect = mean(X(:,isCorrect), 2);
  summaries(sessionIndex).meanError = mean(X(:,~isCorrect), 2);
  summaries(sessionIndex).nCorrect = sum(isCorrect);
  summaries(sessionIndex).nError = sum(~isCorrect);
  summaries(sessionIndex).canTrain = true;
end
end

%% ------------------------------------------------------------------------
function [kernel, informative] = estimateLOOKernel(A, summaries, included, ...
    referenceAllocationPC, method)

switch method
  case "correct-error"
    nCorrect = [summaries(included).nCorrect];
    nError = [summaries(included).nError];
    correctSum = zeros(size(summaries(1).meanCorrect));
    errorSum = zeros(size(summaries(1).meanError));
    for k = included(:)'
      correctSum = correctSum + summaries(k).nCorrect .* summaries(k).meanCorrect;
      errorSum = errorSum + summaries(k).nError .* summaries(k).meanError;
    end
    kernel = correctSum ./ sum(nCorrect) - errorSum ./ sum(nError);
  case {"probit", "logit"}
    kernel = fitBinnedGLMKernel(A(included), referenceAllocationPC, method);
  otherwise
    error('fitIDJModelRecoveryAcrossSessions:UnknownKernelMethod', ...
      'Unknown kernel method %s.', method);
end

if sum(kernel) < 0
  kernel = -kernel;
end
mass = sum(abs(kernel));
informative = isfinite(mass) && mass > 0;
if ~informative
  kernel = ones(numel(A{1}.stepTimeMS),1)/numel(A{1}.stepTimeMS);
  return
end
kernel = double(kernel(:)) ./ mass;
end

%% ------------------------------------------------------------------------
function frameKernel = fitBinnedGLMKernel(A, referenceAllocationPC, link)

nSessions = numel(A);
nFrames = numel(A{1}.stepTimeMS);
nBins = 5;
assert(mod(nFrames, nBins) == 0, ...
  'fitIDJModelRecoveryAcrossSessions:KernelBinCount', ...
  'The number of step frames must divide evenly into five 50-ms bins.');
framesPerBin = nFrames / nBins;

X = [];
y = [];
sessionIndex = [];
for k = 1:nSessions
  T = A{k}.trialTable;
  rows = find(abs(double(T.unsignedAllocPC) - referenceAllocationPC) < 1e-6);
  noise = nan(numel(rows), nFrames);
  for j = 1:numel(rows)
    stream = double(T.highStreamIndex(rows(j)));
    noise(j,:) = double(A{k}.changeNoiseByFrameTrial(:,rows(j),stream))';
  end
  informative = any(abs(noise)>1e-8,2);
  noise = noise(informative,:);
  rows = rows(informative);
  binned = nan(numel(rows), nBins);
  for bin = 1:nBins
    frames = (bin - 1) * framesPerBin + (1:framesPerBin);
    binned(:,bin) = mean(noise(:,frames), 2);
  end
  X = [X; double(T.totalStepCohPC(rows)) binned]; %#ok<AGROW>
  y = [y; double(T.correct(rows))]; %#ok<AGROW>
  sessionIndex = [sessionIndex; repmat(k, numel(rows), 1)]; %#ok<AGROW>
end

sessionTerms = double(sessionIndex == (1:nSessions));
design = [sessionTerms X];
scale = std(design, 0, 1);
scale(1:nSessions) = 1;
scale(scale == 0) = 1;
Z = design ./ scale;

theta0 = zeros(size(Z,2), 1);
theta0(nSessions + 1) = 0.5;
objective = @(theta) glmNLL(theta, Z, y, link);
optimOptions = optimset('Display', 'off', 'MaxFunEvals', 50000, ...
  'MaxIter', 10000, 'TolFun', 1e-8, 'TolX', 1e-8);
[theta, nll, exitflag] = fminsearch(objective, theta0, optimOptions); %#ok<ASGLU>
assert(isfinite(nll) && exitflag ~= -1, ...
  'fitIDJModelRecoveryAcrossSessions:KernelGLMFailure', ...
  'The %s temporal-kernel fit failed.', link);
beta = theta ./ scale(:);
binKernel = beta(nSessions + 2:end);

frameKernel = repelem(binKernel ./ framesPerBin, framesPerBin);
end

%% ------------------------------------------------------------------------
function nll = glmNLL(theta, X, y, link)

eta = X * theta;
if link == "probit"
  probability = gaussianCDF(eta);
else
  eta = min(max(eta, -40), 40);
  probability = 1 ./ (1 + exp(-eta));
end
nll = -bernoulliLogLikelihood(probability, y);
end

%% ------------------------------------------------------------------------
function decisionEvidence = modelEvidence(A, kernel, modelName)

T = A.trialTable;
n = height(T);
allocation = [double(T.alloc0) double(T.alloc1) double(T.alloc2)];
stepByStream = allocation .* double(T.totalStepCohPC);
pedestal = double(T.prestepCohPC);

changeStream = nan(n, 3);
noChangeStream = nan(n, 3);
for stream = 1:3
  changeAbsolute = double(A.changeNoiseByFrameTrial(:,:,stream)) + ...
    reshape(pedestal + stepByStream(:,stream), 1, n);
  noChangeAbsolute = double(A.noChangeNoiseByFrameTrial(:,:,stream)) + ...
    reshape(pedestal, 1, n);
  changeStream(:,stream) = sum(changeAbsolute .* kernel, 1)';
  noChangeStream(:,stream) = sum(noChangeAbsolute .* kernel, 1)';
end

switch modelName
  case "L1"
    changePooled = sum(abs(changeStream), 2);
    noChangePooled = sum(abs(noChangeStream), 2);
  case "L2"
    changePooled = sqrt(sum(changeStream .^ 2, 2));
    noChangePooled = sqrt(sum(noChangeStream .^ 2, 2));
  case "Linf"
    changePooled = max(abs(changeStream), [], 2);
    noChangePooled = max(abs(noChangeStream), [], 2);
  case "Square"
    changePooled = sum(changeStream .^ 2, 2);
    noChangePooled = sum(noChangeStream .^ 2, 2);
  otherwise
    error('fitIDJModelRecoveryAcrossSessions:UnknownModel', ...
      'Unknown pooling model %s.', modelName);
end
decisionEvidence = changePooled - noChangePooled;
end

%% ------------------------------------------------------------------------
function sigma = fitDecisionSigma(decisionEvidence, correct)

scale = max(std(decisionEvidence), 1);
objective = @(logSigma) -bernoulliLogLikelihood( ...
  gaussianCDF(decisionEvidence ./ exp(logSigma)), correct);
sigmaBounds = log([scale * 1e-3 scale * 1e3]);
logSigma = fminbnd(objective, sigmaBounds(1), sigmaBounds(2), ...
  optimset('Display', 'off', 'TolX', 1e-8));
sigma = exp(logSigma);
end

%% ------------------------------------------------------------------------
function probability = gaussianCDF(z)

probability = 0.5 .* erfc(-double(z) ./ sqrt(2));
probability = min(max(probability, 1e-12), 1 - 1e-12);
end

%% ------------------------------------------------------------------------
function logLikelihood = bernoulliLogLikelihood(probability, correct)

probability = min(max(double(probability(:)), 1e-12), 1 - 1e-12);
correct = double(correct(:));
logLikelihood = sum(correct .* log(probability) + ...
  (1 - correct) .* log(1 - probability));
end

%% ------------------------------------------------------------------------
function G = makePsychometricTable(P)

keys = unique([P.unsignedAllocPC P.totalStepCohPC], 'rows');
n = size(keys, 1);
unsignedAllocPC = keys(:,1);
totalStepCohPC = keys(:,2);
nTrials = nan(n,1);
observedPCorrect = nan(n,1);
predictedPCorrectL1 = nan(n,1);
predictedPCorrectL2 = nan(n,1);
predictedPCorrectLinf = nan(n,1);
predictedPCorrectSquare = nan(n,1);
for k = 1:n
  idx = abs(P.unsignedAllocPC - unsignedAllocPC(k)) < 1e-6 & ...
    abs(P.totalStepCohPC - totalStepCohPC(k)) < 1e-6;
  nTrials(k) = sum(idx);
  observedPCorrect(k) = mean(P.correct(idx));
  predictedPCorrectL1(k) = mean(P.pL1(idx));
  predictedPCorrectL2(k) = mean(P.pL2(idx));
  predictedPCorrectLinf(k) = mean(P.pLinf(idx));
  predictedPCorrectSquare(k) = mean(P.pSquare(idx));
end
G = table(unsignedAllocPC, totalStepCohPC, nTrials, observedPCorrect, ...
  predictedPCorrectL1, predictedPCorrectL2, predictedPCorrectLinf, ...
  predictedPCorrectSquare);
G = sortrows(G, {'unsignedAllocPC','totalStepCohPC'}, {'descend','ascend'});
end

%% ------------------------------------------------------------------------
function S = subsetRecoveryResult(result, animalName)
% Preserve the result contract while restricting rows to one animal.
keep = result.sessionTable.animal == animalName;
S = result;
S = rmfield(S, 'animalResults');
S.animal = animalName;
S.sourceFiles = result.sourceFiles(keep);
S.sessionTable = result.sessionTable(keep,:);
S.kernelByFrameSession = result.kernelByFrameSession(:,keep);
S.kernelSourceBySession = result.kernelSourceBySession(keep);
S.kernelTrainingSessionCount = result.kernelTrainingSessionCount(keep);
S.meanLOOKernel = mean(S.kernelByFrameSession, 2);
S.sigmaBySessionModel = result.sigmaBySessionModel(keep,:);
S.fitLogLikelihood = result.fitLogLikelihood(keep,:);
S.calibrationLogLikelihood = result.calibrationLogLikelihood(keep,:);
S.referenceLogLikelihood = result.referenceLogLikelihood(keep,:);
S.dividedLogLikelihood = result.dividedLogLikelihood(keep,:);
S.scoredLogLikelihood = result.scoredLogLikelihood(keep,:);
S.totalScoredLogLikelihood = sum(S.scoredLogLikelihood, 1);
S.heldOutLogLikelihood = [];
S.totalHeldOutLogLikelihood = [];
if result.sigmaFitMode == "reference"
  S.heldOutLogLikelihood = result.heldOutLogLikelihood(keep,:);
  S.totalHeldOutLogLikelihood = sum(S.heldOutLogLikelihood, 1);
end
S.deltaLogLikelihoodVsLinf = S.totalScoredLogLikelihood - ...
  S.totalScoredLogLikelihood(3);
S.noiseAuditTable = result.noiseAuditTable(keep,:);
S.excludedNoiseFreeSessions = result.excludedNoiseFreeSessions( ...
  result.excludedNoiseFreeSessions.animal==animalName,:);
S.predictionTable = result.predictionTable( ...
  ismember(result.predictionTable.sourceSessionID, ...
  S.sessionTable.sourceSessionID),:);
S.psychometricTable = makePsychometricTable(S.predictionTable);
S.groupTable = result.groupTable(result.groupTable.group == animalName,:);
end

%% ------------------------------------------------------------------------
function groupTable = makeModelComparisonTable(result)
% Equal-weight session means; pooledPer1000 weights by scored trial count.
T = result.sessionTable;
LL = result.scoredLogLikelihood;
n = double(T.nScoredTrials);
comparisonColumns = [2 1 4];
comparisonLabels = ["Max - L2" "Max - L1" "Max - Square"];
groups = unique(T.animal, 'stable');
if numel(groups) > 1
  groups(end+1) = "Combined";
end
group = strings(0,1);
comparison = strings(0,1);
nSessions = zeros(0,1);
nDividedTrials = zeros(0,1);
nScoredTrials = zeros(0,1);
nPositive = zeros(0,1);
meanPer1000 = zeros(0,1);
sdPer1000 = zeros(0,1);
lower95 = zeros(0,1);
upper95 = zeros(0,1);
lowerOneSided95 = zeros(0,1);
oneSidedP = zeros(0,1);
totalDeltaLL = zeros(0,1);
pooledPer1000 = zeros(0,1);
for k = 1:numel(groups)
  keep = T.animal == groups(k);
  if groups(k) == "Combined"
    keep = true(height(T),1);
  end
  for model = 1:numel(comparisonColumns)
    d = LL(keep,3) - LL(keep,comparisonColumns(model));
    label = comparisonLabels(model);
    rate = 1000 * d ./ n(keep);
    [mu, sd, lo, hi, loOneSided, pValue] = sessionT(rate);
    group(end+1,1) = groups(k);
    comparison(end+1,1) = label;
    nSessions(end+1,1) = sum(keep);
    nDividedTrials(end+1,1) = sum(T.nDividedTrials(keep));
    nScoredTrials(end+1,1) = sum(n(keep));
    nPositive(end+1,1) = sum(d > 0);
    meanPer1000(end+1,1) = mu;
    sdPer1000(end+1,1) = sd;
    lower95(end+1,1) = lo;
    upper95(end+1,1) = hi;
    lowerOneSided95(end+1,1) = loOneSided;
    oneSidedP(end+1,1) = pValue;
    totalDeltaLL(end+1,1) = sum(d);
    pooledPer1000(end+1,1) = 1000 * sum(d) / sum(n(keep));
  end
end
groupTable = table(group, comparison, nSessions, nDividedTrials, ...
  nScoredTrials, ...
  nPositive, meanPer1000, sdPer1000, lower95, upper95, ...
  lowerOneSided95, oneSidedP, totalDeltaLL, pooledPer1000);
groupTable.scoreDescription = repmat(result.scoreDescription,height(groupTable),1);
end

%% ------------------------------------------------------------------------
function [mu, sd, lo, hi, loOneSided, pValue] = sessionT(rate)
% Provisional one-sided session t test; no inference with one session.
rate = double(rate(:));
n = numel(rate);
mu = mean(rate);
if n < 2
  sd = NaN; lo = NaN; hi = NaN;
  loOneSided = NaN; pValue = NaN;
  return;
end
sd = std(rate, 0);
se = sd / sqrt(n);
if se == 0
  lo = mu; hi = mu; loOneSided = mu;
  pValue = double(mu < 0) + 0.5 * double(mu == 0);
  return;
end
df = n - 1;
pValue = tUpperTail(mu / se, df);
twoSidedCritical = fzero(@(z) tUpperTail(z, df) - 0.025, [0 1000]);
oneSidedCritical = fzero(@(z) tUpperTail(z, df) - 0.05, [0 1000]);
lo = mu - twoSidedCritical * se;
hi = mu + twoSidedCritical * se;
loOneSided = mu - oneSidedCritical * se;
end

function p = tUpperTail(t, df)
halfTail = 0.5 * betainc(df / (df + t.^2), df/2, 0.5);
if t >= 0
  p = halfTail;
else
  p = 1 - halfTail;
end
end

%% ------------------------------------------------------------------------
function fig = plotRecoverySummary(result, visible)

colors = [0.10 0.45 0.80; 0.90 0.45 0.10; 0.20 0.60 0.35; ...
  0.75 0.20 0.55];
fig = figure('Color', 'w', 'Visible', visible, ...
  'Position', [100 100 1600 980]);
tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
kernelDescription = char(result.kernelMethod);
boxcar = startsWith(result.kernelSourceBySession,"uniform boxcar");
if all(boxcar)
  kernelDescription = 'uniform boxcar';
elseif any(boxcar)
  kernelDescription = [kernelDescription ' / boxcar'];
end
title(tl, sprintf('IDJ model recovery: %s (%s kernel)', ...
  char(result.animal), kernelDescription), 'FontWeight', 'bold');

ax = nexttile(tl, 1);
axis(ax, 'off');
info = sprintf(['%s\n\nSessions: %d\nAnalysis trials: %d\n' ...
  'Sigma fit: %s (%d trials)\nScored trials: %d\n%s\n\nPrepared: %s\n\n' ...
  'Kernel: noisy reference trials, same-animal LOO\nBoxcar sessions: %d\n' ...
  'Trials without step noise: %d\nExcluded noise-free sessions: %d'], ...
  char(result.animal), height(result.sessionTable), ...
  sum(result.sessionTable.nReferenceTrials) + ...
  sum(result.sessionTable.nDividedTrials), ...
  char(result.sigmaFitMode),sum(result.sessionTable.nFitTrials), ...
  sum(result.sessionTable.nScoredTrials),char(result.scoreDescription), ...
  datestr(result.createdAt, 'yyyy-mm-dd HH:MM'), ...
  sum(boxcar), ...
  sum(result.noiseAuditTable.nZeroNoiseTrials), ...
  height(result.excludedNoiseFreeSessions));
text(ax, 0.02, 0.98, info, 'Units', 'normalized', ...
  'VerticalAlignment', 'top', 'Interpreter', 'none', 'FontSize', 10);

ax = nexttile(tl, 2);
advantage = result.scoredLogLikelihood(:,3) - ...
  result.scoredLogLikelihood(:,[1 2 4]);
plotSessionAdvantages(ax, result, advantage, false);

ax = nexttile(tl, 3);
plotSessionAdvantages(ax, result, advantage, true);

ax = nexttile(tl, 6);
plot(ax, result.stepTimeMS, result.kernelByFrameSession, '-', ...
  'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
hold(ax, 'on');
plot(ax, result.stepTimeMS, result.meanLOOKernel, 'k-', 'LineWidth', 2);
xlabel(ax, 'Time from step onset (ms)');
ylabel(ax, 'Normalized weight');
title(ax, 'Temporal kernels (same-animal LOO / boxcar)');
grid(ax, 'on'); box(ax, 'off');

allocations = [100 70 50 100/3];
tiles = [4 5 7 8];
G = result.psychometricTable;
commonMaxStep = max(5, 5 * ceil(max(G.totalStepCohPC) / 5));
for allocationIndex = 1:numel(allocations)
  ax = nexttile(tl, tiles(allocationIndex));
  rows = abs(G.unsignedAllocPC - allocations(allocationIndex)) < 1e-6;
  x = G.totalStepCohPC(rows);
  hold(ax, 'on');
  plot(ax, x, G.observedPCorrect(rows), 'ko', 'MarkerFaceColor', 'k', ...
    'DisplayName', 'Observed');
  plot(ax, x, G.predictedPCorrectL1(rows), '-', 'Color', colors(1,:), ...
    'LineWidth', 1.5, 'DisplayName', 'L1');
  plot(ax, x, G.predictedPCorrectL2(rows), '-', 'Color', colors(2,:), ...
    'LineWidth', 1.5, 'DisplayName', 'L2');
  plot(ax, x, G.predictedPCorrectLinf(rows), '-', 'Color', colors(3,:), ...
    'LineWidth', 2, 'DisplayName', 'Linf');
  plot(ax, x, G.predictedPCorrectSquare(rows), '--', 'Color', colors(4,:), ...
    'LineWidth', 2, 'DisplayName', 'Square');
  yline(ax, 0.75, ':', 'HandleVisibility', 'off');
  xlabel(ax, 'Total coherence step (%)');
  ylabel(ax, 'P(correct)');
  xlim(ax, [0 commonMaxStep]);
  xticks(ax, 0:5:commonMaxStep);
  ylim(ax, [0.45 1.02]);
  if result.sigmaFitMode == "all"
    suffix = ' (fitted)';
  elseif allocations(allocationIndex) == result.referenceAllocationPC
    suffix = ' (sigma calibration)';
  else
    suffix = ' (held out)';
  end
  title(ax, sprintf('%s%s',idjAllocationLabel(allocations(allocationIndex)),suffix));
  grid(ax, 'on');
  box(ax, 'off');
  if allocationIndex == 1
    legend(ax, 'Location', 'southeast');
  end
end

ax = nexttile(tl, 9);
axis(ax, 'off');
G = result.groupTable;
G = G(G.group == result.animal,:);
lines = sprintf('Session-level model comparison\n%s\n\n', ...
  char(result.scoreDescription));
for k = 1:height(G)
  if isnan(G.oneSidedP(k))
    inference = 'CI and p: available after 2 sessions';
  else
    inference = sprintf('95%% CI [%.2f, %.2f], one-sided p = %.3g', ...
      G.lower95(k), G.upper95(k), G.oneSidedP(k));
  end
  lines = [lines sprintf(['%s: %d/%d sessions favor max\n' ...
    '  Mean %.2f LL/1000; total %+.1f\n  %s\n\n'], ...
    char(G.comparison(k)), G.nPositive(k), G.nSessions(k), ...
    G.meanPer1000(k), G.totalDeltaLL(k), inference)]; %#ok<AGROW>
end
if result.animal == "Combined" && numel(unique(result.sessionTable.animal)) > 1
  lines = [lines 'Combined p treats sessions, not animals, as units.'];
end
text(ax, 0.02, 0.98, lines, 'Units', 'normalized', ...
  'VerticalAlignment', 'top', 'Interpreter', 'none', 'FontSize', 9);
end

%% ------------------------------------------------------------------------
function plotSessionAdvantages(ax, result, advantage, cumulative)

n = height(result.sessionTable);
x = 1:n;
if cumulative
  values = cumsum(advantage, 1);
  if result.sigmaFitMode == "all"
    ylabel(ax, 'Cumulative fitted LL advantage');
  else
    ylabel(ax, 'Cumulative held-out LL advantage');
  end
  title(ax, sprintf('Cumulative: L1 %+.1f; L2 %+.1f; Square %+.1f', ...
    values(end,1), values(end,2), values(end,3)));
else
  values = advantage;
  if result.sigmaFitMode == "all"
    ylabel(ax, 'Session fitted LL advantage');
  else
    ylabel(ax, 'Session held-out LL advantage');
  end
  title(ax, 'Individual sessions');
end
isCombined = numel(unique(result.sessionTable.animal)) > 1;
plot(ax, x, values(:,1), 'o-', ...
  'LineWidth', 1.2, 'DisplayName', 'Max - L1');
hold(ax, 'on');
plot(ax, x, values(:,2), 's-', ...
  'LineWidth', 1.2, 'DisplayName', 'Max - L2');
plot(ax, x, values(:,3), 'd-', ...
  'LineWidth', 1.2, 'DisplayName', 'Max - Square');
yline(ax, 0, ':', 'HandleVisibility', 'off');
labels = string(result.sessionTable.sourceSessionID);
for k = 1:n
  date = regexp(char(labels(k)), '[0-9]{8}$', 'match', 'once');
  if ~isempty(date)
    labels(k) = string(date(5:6)) + "/" + string(date(7:8));
  end
  if isCombined
    animalName = char(result.sessionTable.animal(k));
    labels(k) = string(animalName(1)) + " " + labels(k);
  end
end
if isCombined
  boundaries = find(result.sessionTable.animal(2:end) ~= ...
    result.sessionTable.animal(1:end-1));
  for k = boundaries(:)'
    xline(ax, k + 0.5, ':', 'HandleVisibility', 'off');
  end
end
ticks = 1:max(1, ceil(n/12)):n;
set(ax, 'XTick', ticks, 'XTickLabel', labels(ticks));
xtickangle(ax, 45);
xlabel(ax, 'Session date');
legend(ax, 'Location', 'best');
grid(ax, 'on'); box(ax, 'off');
end
