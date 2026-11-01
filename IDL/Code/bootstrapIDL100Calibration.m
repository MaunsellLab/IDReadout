function [groupTable, correlationTable, bootstrap] = ...
    bootstrapIDL100Calibration(input, varargin)
% bootstrapIDL100Calibration  Effect of 100:0 sampling on divided-trial LL.
%
% [groupTable, correlationTable, bootstrap] = ...
%   bootstrapIDL100Calibration(result, 'NumBootstrap', 500, 'Seed', 1)
%
% Resamples WHOLE 100:0 trials within each coherence step and session.
% The same draw is used to refit all three model sigmas. Each session's
% saved temporal kernel, divided-trial predictors, and divided choices
% remain fixed. This isolates uncertainty from the session's own 100:0
% sigma calibration. It does not re-estimate the LOO kernels or include
% sampling variation in the divided trials.
%
% Source sessionAnalysis files in result.sourceFiles must be accessible.
% The function checks that recomputed probabilities reproduce the saved
% predictions before bootstrapping, preventing use of a stale result.

p = inputParser;
p.addParameter('NumBootstrap', 500, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2 && x == floor(x));
p.addParameter('Seed', 1, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x == floor(x));
p.parse(varargin{:});
B = double(p.Results.NumBootstrap);

if ischar(input) || (isstring(input) && isscalar(input))
  saved = load(char(input), 'result');
  result = saved.result;
else
  result = input;
end
assert(isequal(string(result.modelNames(:))', ["L1" "L2" "Linf"]), ...
  'bootstrapIDL100Calibration:ModelOrder', ...
  'Expected result.modelNames to be L1, L2, Linf.');

S = result.sessionTable;
P = result.predictionTable;
nSessions = height(S);
allocations = sort(unique(double(P.unsignedAllocPC( ...
  ~P.isCalibrationTrial))), 'descend');
nAllocations = numel(allocations);
assert(nAllocations >= 1, 'bootstrapIDL100Calibration:NoDividedTrials', ...
  'The result has no divided-allocation trials.');

% Dimensions: bootstrap replicate, session, allocation, comparison.
% Comparisons are Max-L2 and Max-L1, in that order.
rates = nan(B, nSessions, nAllocations, 2);
observedRates = nan(nSessions, nAllocations, 2);
nDivided = zeros(nSessions, nAllocations);
sigmas = nan(B, nSessions, 3);

oldRNG = rng;
restoreRNG = onCleanup(@() rng(oldRNG)); %#ok<NASGU>
rng(double(p.Results.Seed), 'twister');

for k = 1:nSessions
  source = char(result.sourceFiles(k));
  if ~isfile(source)
    [~, basename, extension] = fileparts(source);
    source = fullfile(char(result.sourceFolder), [basename extension]);
  end
  assert(isfile(source), 'bootstrapIDL100Calibration:MissingSessionFile', ...
    'Cannot find the saved sessionAnalysis source file: %s', source);
  loaded = load(source, 'sessionAnalysis');
  A = loaded.sessionAnalysis;
  T = A.trialTable;
  assert(string(A.sourceSessionID) == string(S.sourceSessionID(k)), ...
    'bootstrapIDL100Calibration:SessionOrder', ...
    'Source file and saved session table differ at session %d.', k);
  savedRows = string(P.sourceSessionID) == string(S.sourceSessionID(k));
  Q = P(savedRows,:);
  assert(height(Q) == height(T) && ...
    isequal(double(Q.sourceTrialNumber), double(T.sourceTrialNumber)) && ...
    isequal(logical(Q.correct), logical(T.correct)) && ...
    isequal(double(Q.unsignedAllocPC), double(T.unsignedAllocPC)) && ...
    isequal(double(Q.totalStepCohPC), double(T.totalStepCohPC)), ...
    'bootstrapIDL100Calibration:TrialOrder', ...
    'Saved prediction rows do not match source trials in session %d.', k);

  evidence = nan(height(T), 3);
  probabilityColumns = {'pL1','pL2','pLinf'};
  models = ["L1" "L2" "Linf"];
  for model = 1:3
    evidence(:,model) = modelEvidence(A, ...
      result.kernelByFrameSession(:,k), models(model));
    reconstructed = gaussianCDF(evidence(:,model) ./ ...
      result.sigmaBySessionModel(k,model));
    discrepancy = max(abs(reconstructed - double(Q.(probabilityColumns{model}))));
    assert(isfinite(discrepancy) && discrepancy < 1e-8, ...
      'bootstrapIDL100Calibration:PredictionMismatch', ...
      'Recomputed predictions differ in session %d, model %s (max error %.3g).', ...
      k, models(model), discrepancy);
  end

  calibration = logical(Q.isCalibrationTrial);
  divided = ~calibration;
  referenceRows = find(calibration);
  yReference = double(T.correct(referenceRows));
  XReference = evidence(referenceRows,:);
  step = double(T.totalStepCohPC(referenceRows));
  stepValues = unique(step);
  stepGroups = cell(numel(stepValues),1);
  for h = 1:numel(stepValues)
    stepGroups{h} = find(abs(step - stepValues(h)) < 1e-6);
  end

  dividedRows = find(divided);
  yDivided = double(T.correct(dividedRows));
  XDivided = evidence(dividedRows,:);
  allocDivided = double(T.unsignedAllocPC(dividedRows));
  allocationRows = cell(nAllocations,1);
  for j = 1:nAllocations
    allocationRows{j} = find(abs(allocDivided - allocations(j)) < 1e-6);
    nDivided(k,j) = numel(allocationRows{j});
    assert(nDivided(k,j) > 0, ...
      'bootstrapIDL100Calibration:EmptyAllocation', ...
      'Session %d has no trials at allocation %g.', k, allocations(j));
  end

  originalLL = trialLogLikelihood( ...
    XDivided, yDivided, result.sigmaBySessionModel(k,:));
  observedRates(k,:,:) = allocationRates( ...
    originalLL, allocationRows, nDivided(k,:));

  for b = 1:B
    drawn = zeros(numel(referenceRows),1);
    cursor = 0;
    for h = 1:numel(stepGroups)
      members = stepGroups{h};
      n = numel(members);
      drawn(cursor+(1:n)) = members(randi(n,n,1));
      cursor = cursor + n;
    end
    sigmaBoot = nan(1,3);
    for model = 1:3
      sigmaBoot(model) = fitDecisionSigma( ...
        XReference(drawn,model), yReference(drawn));
    end
    sigmas(b,k,:) = reshape(sigmaBoot, 1, 1, 3);
    ll = trialLogLikelihood(XDivided, yDivided, sigmaBoot);
    rates(b,k,:,:) = reshape( ...
      allocationRates(ll, allocationRows, nDivided(k,:)), ...
      1, 1, nAllocations, 2);
  end
end

comparisonNames = ["Max - L2" "Max - L1"];
animals = unique(string(S.animal), 'stable');
group = strings(0,1);
comparison = strings(0,1);
allocationPC = zeros(0,1);
nGroupSessions = zeros(0,1);
meanCalibrationTrials = zeros(0,1);
meanObservedPer1000 = zeros(0,1);
observedSessionSDPer1000 = zeros(0,1);
rmsCalibrationSEPer1000 = zeros(0,1);
calibrationVarianceOverObserved = zeros(0,1);

pairGroup = strings(0,1);
pairComparison = strings(0,1);
allocationA = zeros(0,1);
allocationB = zeros(0,1);
calibrationCorrelation = zeros(0,1);

for g = 1:numel(animals)
  sessions = find(string(S.animal) == animals(g));
  for c = 1:2
    covariance = nan(nAllocations);
    for a = 1:nAllocations
      xa = reshape(rates(:,sessions,a,c), B, []);
      xa = xa - mean(xa,1);
      for b = 1:nAllocations
        xb = reshape(rates(:,sessions,b,c), B, []);
        xb = xb - mean(xb,1);
        covariance(a,b) = sum(xa .* xb, 'all') / ...
          ((B - 1) * numel(sessions));
      end
      observed = observedRates(sessions,a,c);
      group(end+1,1) = animals(g);
      comparison(end+1,1) = comparisonNames(c);
      allocationPC(end+1,1) = allocations(a);
      nGroupSessions(end+1,1) = numel(sessions);
      meanCalibrationTrials(end+1,1) = ...
        mean(double(S.nCalibrationTrials(sessions)));
      meanObservedPer1000(end+1,1) = mean(observed);
      observedSessionSDPer1000(end+1,1) = std(observed);
      rmsCalibrationSEPer1000(end+1,1) = sqrt(covariance(a,a));
      calibrationVarianceOverObserved(end+1,1) = ...
        covariance(a,a) / observedSessionSDPer1000(end)^2;
    end
    for a = 1:nAllocations
      for b = a+1:nAllocations
        pairGroup(end+1,1) = animals(g);
        pairComparison(end+1,1) = comparisonNames(c);
        allocationA(end+1,1) = allocations(a);
        allocationB(end+1,1) = allocations(b);
        calibrationCorrelation(end+1,1) = covariance(a,b) / ...
          sqrt(covariance(a,a) * covariance(b,b));
      end
    end
  end
end

groupTable = table(group, comparison, allocationPC, nGroupSessions, ...
  meanCalibrationTrials, meanObservedPer1000, observedSessionSDPer1000, ...
  rmsCalibrationSEPer1000, calibrationVarianceOverObserved);
correlationTable = table(pairGroup, pairComparison, allocationA, ...
  allocationB, calibrationCorrelation, 'VariableNames', ...
  {'group','comparison','allocationA','allocationB', ...
  'calibrationCorrelation'});

bootstrap = struct();
bootstrap.numBootstrap = B;
bootstrap.seed = p.Results.Seed;
bootstrap.sessionIDs = string(S.sourceSessionID);
bootstrap.animals = string(S.animal);
bootstrap.allocations = allocations;
bootstrap.comparisons = comparisonNames;
bootstrap.observedRatePer1000 = observedRates;
bootstrap.ratePer1000 = rates;
bootstrap.sigmaByReplicateSessionModel = sigmas;
bootstrap.nDividedTrials = nDivided;
bootstrap.contract = ['Whole 100:0 trials resampled within step; ' ...
  'refit 3 sigmas; kernels and divided trials held fixed'];

if nargout == 0
  disp(groupTable);
  disp(correlationTable);
end
end

function rates = allocationRates(ll, allocationRows, nDivided)
% ll has columns L1, L2, Max.
rates = nan(1,numel(allocationRows),2);
for j = 1:numel(allocationRows)
  rows = allocationRows{j};
  rates(1,j,1) = 1000 * sum(ll(rows,3) - ll(rows,2)) / nDivided(j);
  rates(1,j,2) = 1000 * sum(ll(rows,3) - ll(rows,1)) / nDivided(j);
end
end

function ll = trialLogLikelihood(evidence, correct, sigma)
ll = nan(size(evidence));
for model = 1:3
  probability = gaussianCDF(evidence(:,model) ./ sigma(model));
  ll(:,model) = correct .* log(probability) + ...
    (1 - correct) .* log(1 - probability);
end
end

function sigma = fitDecisionSigma(evidence, correct)
% Match fitIDLModelRecoveryAcrossSessions exactly.
scale = max(std(evidence), 1);
objective = @(logSigma) -sumBernoulli( ...
  gaussianCDF(evidence ./ exp(logSigma)), correct);
sigmaBounds = log([scale * 1e-3 scale * 1e3]);
logSigma = fminbnd(objective, sigmaBounds(1), sigmaBounds(2), ...
  optimset('Display', 'off', 'TolX', 1e-8));
sigma = exp(logSigma);
end

function ll = sumBernoulli(probability, correct)
ll = sum(correct .* log(probability) + ...
  (1 - correct) .* log(1 - probability));
end

function probability = gaussianCDF(z)
probability = 0.5 .* erfc(-double(z) ./ sqrt(2));
probability = min(max(probability, 1e-12), 1 - 1e-12);
end

function evidence = modelEvidence(A, kernel, modelName)
% Match the analysis evidence calculation; a single kernel weights both
% streams and both patches before separate within-patch pooling.
T = A.trialTable;
n = height(T);
allocation = [double(T.alloc0) double(T.alloc1)];
stepByStream = allocation .* double(T.totalStepCohPC);
pedestal = double(T.prestepCohPC);
changeStream = nan(n,2);
noChangeStream = nan(n,2);
for stream = 1:2
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
    changePooled = sqrt(sum(changeStream.^2, 2));
    noChangePooled = sqrt(sum(noChangeStream.^2, 2));
  case "Linf"
    changePooled = max(abs(changeStream), [], 2);
    noChangePooled = max(abs(noChangeStream), [], 2);
  otherwise
    error('bootstrapIDL100Calibration:UnknownModel', ...
      'Unknown model %s.', modelName);
end
evidence = changePooled - noChangePooled;
end
