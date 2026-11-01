function [groupTable, sessionNoiseTable] = summarizeIDLTrialSamplingNoise(input)
% summarizeIDLTrialSamplingNoise  Trial sampling noise in IDL delta LL.
%
% [groupTable, sessionNoiseTable] = summarizeIDLTrialSamplingNoise(result)
% [groupTable, sessionNoiseTable] = ...
%   summarizeIDLTrialSamplingNoise('Data/AcrossSessions/ModelRecovery.mat')
%
% The stratified trial SE holds the number of trials at each sampled step
% fixed and uses the within-step variance of whole-trial LL contributions.
% It includes variation from stimulus-noise draws and choices. The separate
% choice-only SE holds the recorded predictors fixed. For each trial, the
% Max-minus-alternative LL contribution has two values:
% one if correct and one if incorrect. Given success probability q, its
% variance is q*(1-q) times the squared difference of those values.
%
% The choice SE uses a smoothed observed correct fraction for each
% session/allocation/step. Max and alternative model probabilities provide
% sensitivity checks. All SEs exclude uncertainty in the temporal kernel,
% 100:0 sigma, serial dependence, and changes in the trial/step design.
% The residual session SD is a descriptive variance subtraction. It can
% also reflect differences in step grids and fitted predictors, and is
% especially uncertain when an animal has few sessions.

if ischar(input) || (isstring(input) && isscalar(input))
  saved = load(char(input), 'result');
  result = saved.result;
else
  result = input;
end

P = result.predictionTable;
S = result.sessionTable;
[found, sessionIndex] = ismember( ...
  string(P.sourceSessionID), string(S.sourceSessionID));
assert(all(found), 'summarizeIDLTrialSamplingNoise:UnknownSession', ...
  'A prediction-table session is missing from the session table.');
calibration = logical(P.isCalibrationTrial);
allocation = double(P.unsignedAllocPC);
step = double(P.totalStepCohPC);
correct = double(P.correct);
allocations = sort(unique(allocation(~calibration)), 'descend');
assert(~isempty(allocations), ...
  'summarizeIDLTrialSamplingNoise:NoDividedTrials', ...
  'No divided-allocation trials were found.');

% Estimate q at each sampled step without using either pooling model.
% Jeffreys smoothing avoids q=0 or 1 for small or perfect step groups.
qObserved = nan(height(P),1);
for k = 1:height(S)
  for j = 1:numel(allocations)
    condition = ~calibration & sessionIndex == k & ...
      abs(allocation - allocations(j)) < 1e-6;
    stepValues = unique(step(condition));
    for h = 1:numel(stepValues)
      rows = condition & abs(step - stepValues(h)) < 1e-6;
      qObserved(rows) = (sum(correct(rows)) + 0.5) / (sum(rows) + 1);
    end
  end
end
assert(all(isfinite(qObserved(~calibration))), ...
  'summarizeIDLTrialSamplingNoise:MissingEmpiricalProbability', ...
  'Could not estimate a correct fraction for some divided trials.');

pMax = clippedProbability(P.pLinf);
comparisonNames = ["Max - L2", "Max - L1"];
alternatives = [clippedProbability(P.pL2), clippedProbability(P.pL1)];

sessionIDColumn = strings(0,1);
animalColumn = strings(0,1);
comparisonColumn = strings(0,1);
allocationColumn = zeros(0,1);
nTrialsColumn = zeros(0,1);
deltaLLColumn = zeros(0,1);
per1000Column = zeros(0,1);
trialSEColumn = zeros(0,1);
choiceSEColumn = zeros(0,1);
choiceSEMaxColumn = zeros(0,1);
choiceSEAltColumn = zeros(0,1);

for model = 1:2
  pAlt = alternatives(:,model);
  llCorrect = log(pMax ./ pAlt);
  llError = log((1 - pMax) ./ (1 - pAlt));
  swing = llCorrect - llError;
  trialDelta = correct .* llCorrect + (1 - correct) .* llError;
  varianceObserved = qObserved .* (1 - qObserved) .* swing.^2;
  varianceMax = pMax .* (1 - pMax) .* swing.^2;
  varianceAlt = pAlt .* (1 - pAlt) .* swing.^2;

  for k = 1:height(S)
    for j = 1:numel(allocations)
      rows = ~calibration & sessionIndex == k & ...
        abs(allocation - allocations(j)) < 1e-6;
      n = sum(rows);
      if n == 0
        continue
      end
      sessionIDColumn(end+1,1) = string(S.sourceSessionID(k));
      animalColumn(end+1,1) = string(S.animal(k));
      comparisonColumn(end+1,1) = comparisonNames(model);
      allocationColumn(end+1,1) = allocations(j);
      nTrialsColumn(end+1,1) = n;
      deltaLLColumn(end+1,1) = sum(trialDelta(rows));
      per1000Column(end+1,1) = 1000 * deltaLLColumn(end) / n;
      % Resample complete trials within each step while retaining step counts.
      trialVariance = 0;
      stepValues = unique(step(rows));
      for h = 1:numel(stepValues)
        withinStep = rows & abs(step - stepValues(h)) < 1e-6;
        nStep = sum(withinStep);
        if nStep < 2
          trialVariance = NaN;
          break
        end
        trialVariance = trialVariance + nStep * var(trialDelta(withinStep));
      end
      trialSEColumn(end+1,1) = 1000 * sqrt(trialVariance) / n;
      choiceSEColumn(end+1,1) = 1000 * sqrt(sum(varianceObserved(rows))) / n;
      choiceSEMaxColumn(end+1,1) = 1000 * sqrt(sum(varianceMax(rows))) / n;
      choiceSEAltColumn(end+1,1) = 1000 * sqrt(sum(varianceAlt(rows))) / n;
    end
  end
end

sessionNoiseTable = table(sessionIDColumn, animalColumn, ...
  comparisonColumn, allocationColumn, nTrialsColumn, deltaLLColumn, ...
  per1000Column, trialSEColumn, choiceSEColumn, ...
  choiceSEMaxColumn, choiceSEAltColumn, ...
  'VariableNames', {'sourceSessionID','animal','comparison', ...
  'allocationPC','nTrials','deltaLL','per1000', ...
  'trialSamplingSEPer1000','choiceSEPer1000', ...
  'choiceSEPer1000UnderMax', ...
  'choiceSEPer1000UnderAlternative'});

groups = unique(string(S.animal), 'stable');
groupColumn = strings(0,1);
comparisonSummary = strings(0,1);
allocationSummary = zeros(0,1);
nSessions = zeros(0,1);
meanTrialsPerSession = zeros(0,1);
meanPer1000 = zeros(0,1);
observedSDPer1000 = zeros(0,1);
rmsTrialSamplingSEPer1000 = zeros(0,1);
rmsChoiceSEPer1000 = zeros(0,1);
residualSessionSDPer1000 = zeros(0,1);
rmsChoiceSEUnderMax = zeros(0,1);
rmsChoiceSEUnderAlternative = zeros(0,1);

for g = 1:numel(groups)
  for model = 1:2
    for j = 1:numel(allocations)
      rows = sessionNoiseTable.animal == groups(g) & ...
        sessionNoiseTable.comparison == comparisonNames(model) & ...
        abs(sessionNoiseTable.allocationPC - allocations(j)) < 1e-6;
      Q = sessionNoiseTable(rows,:);
      if isempty(Q)
        continue
      end
      groupColumn(end+1,1) = groups(g);
      comparisonSummary(end+1,1) = comparisonNames(model);
      allocationSummary(end+1,1) = allocations(j);
      nSessions(end+1,1) = height(Q);
      meanTrialsPerSession(end+1,1) = mean(Q.nTrials);
      meanPer1000(end+1,1) = mean(Q.per1000);
      rmsTrialSamplingSEPer1000(end+1,1) = ...
        sqrt(mean(Q.trialSamplingSEPer1000.^2));
      rmsChoiceSEPer1000(end+1,1) = sqrt(mean(Q.choiceSEPer1000.^2));
      rmsChoiceSEUnderMax(end+1,1) = ...
        sqrt(mean(Q.choiceSEPer1000UnderMax.^2));
      rmsChoiceSEUnderAlternative(end+1,1) = ...
        sqrt(mean(Q.choiceSEPer1000UnderAlternative.^2));
      if height(Q) < 2
        observedSDPer1000(end+1,1) = NaN;
        residualSessionSDPer1000(end+1,1) = NaN;
      else
        observedSDPer1000(end+1,1) = std(Q.per1000);
        residualSessionSDPer1000(end+1,1) = sqrt(max(0, ...
          observedSDPer1000(end)^2 - rmsTrialSamplingSEPer1000(end)^2));
      end
    end
  end
end

groupTable = table(groupColumn, comparisonSummary, allocationSummary, ...
  nSessions, meanTrialsPerSession, meanPer1000, observedSDPer1000, ...
  rmsTrialSamplingSEPer1000, rmsChoiceSEPer1000, ...
  residualSessionSDPer1000, rmsChoiceSEUnderMax, ...
  rmsChoiceSEUnderAlternative, 'VariableNames', ...
  {'group','comparison','allocationPC','nSessions','meanTrialsPerSession', ...
  'meanPer1000', ...
  'observedSDPer1000','rmsTrialSamplingSEPer1000', ...
  'rmsChoiceSEPer1000', ...
  'residualSessionSDPer1000','rmsChoiceSEUnderMax', ...
  'rmsChoiceSEUnderAlternative'});

if nargout == 0
  disp(groupTable);
end
end

function p = clippedProbability(p)
p = min(max(double(p(:)), 1e-12), 1 - 1e-12);
end
