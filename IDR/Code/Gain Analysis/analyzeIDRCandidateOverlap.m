function O = analyzeIDRCandidateOverlap(D)
% analyzeIDRCandidateOverlap  Rectangular physical-overlap summaries for IDR.
%
% O = analyzeIDRCandidateOverlap(D)
%
% D is the enriched sideGainData struct produced by analyzeIDRSideGains.
% Probe predictors in this analysis are single-candidate physical streams,
% not the effective sum over the two yoked probe directions. All predictors
% are signed mean coherence over the step interval.
%
% The primary change-side margin is
%
%   changeMargin = changeProbe - (stepDeltaCoh + changePreferred)
%
% and a physical crossing is changeMargin > 0. The analogous no-change-side
% margin, noChangeProbe - noChangePreferred, is retained descriptively.
%
% OUTPUT TABLES
%   trialTable            One row per analyzed INC noise trial
%   offsetSummary         One row per probe offset
%   sessionOffsetSummary  One row per parent session x probe offset
%   conditionSummary      One row per offset x noise amplitudes x step


requiredFields = { ...
  'trialTable', 'stepFrames', ...
  'changePrefNoiseByFrameTrial', ...
  'noChangePrefNoiseByFrameTrial', ...
  'changeProbeCandidateNoiseByFrameTrial', ...
  'noChangeProbeCandidateNoiseByFrameTrial'};
for i = 1:numel(requiredFields)
  if ~isfield(D, requiredFields{i})
    error('analyzeIDRCandidateOverlap:MissingField', ...
      'Input sideGainData is missing %s.', requiredFields{i});
  end
end

T = D.trialTable;
stepFrames = D.stepFrames(:);

changePreferred = stepMean(D.changePrefNoiseByFrameTrial, stepFrames);
noChangePreferred = stepMean(D.noChangePrefNoiseByFrameTrial, stepFrames);
changeProbe = stepMean(D.changeProbeCandidateNoiseByFrameTrial, stepFrames);
noChangeProbe = stepMean(D.noChangeProbeCandidateNoiseByFrameTrial, stepFrames);
step = double(T.stepDeltaCohPC);

changeMargin = changeProbe - (step + changePreferred);
noChangeMargin = noChangeProbe - noChangePreferred;
changeCrossed = changeMargin > 0;
noChangeCrossed = noChangeMargin > 0;
changeEqual = changeMargin == 0;
noChangeEqual = noChangeMargin == 0;

trialTable = T(:, intersect({ ...
  'sessionIndex','sessionID','animal','probeSessionIndex','probeDirDeg', ...
  'nYokedProbeStreams','trialIdx','stepCohPC','stepDeltaCohPC', ...
  'prefCohNoisePC','probeCohNoisePC','correct','changeSide', ...
  'sessionThreshold','sourceFile'}, T.Properties.VariableNames, 'stable'));
trialTable.changePreferred = changePreferred;
trialTable.changeProbeCandidate = changeProbe;
trialTable.changeMargin = changeMargin;
trialTable.changeCrossed = changeCrossed;
trialTable.changeEqual = changeEqual;
trialTable.noChangePreferred = noChangePreferred;
trialTable.noChangeProbeCandidate = noChangeProbe;
trialTable.noChangeMargin = noChangeMargin;
trialTable.noChangeCrossed = noChangeCrossed;
trialTable.noChangeEqual = noChangeEqual;

valid = isfinite(trialTable.probeDirDeg) & ...
  isfinite(trialTable.stepDeltaCohPC) & ...
  isfinite(changePreferred) & isfinite(changeProbe) & ...
  isfinite(noChangePreferred) & isfinite(noChangeProbe);
trialTable = trialTable(valid,:);

offsetGroups = findgroups(trialTable.probeDirDeg);
offsetSummary = groupedSummary(trialTable, offsetGroups, ...
  {'probeDirDeg'});

[sessionOffsetGroups, sessionIndex, sessionID, animal, probeDirDeg] = ...
  findgroups(trialTable.sessionIndex, trialTable.sessionID, ...
  trialTable.animal, trialTable.probeDirDeg);
sessionOffsetSummary = groupedSummary(trialTable, sessionOffsetGroups, ...
  {'sessionIndex','sessionID','animal','probeDirDeg'}, ...
  sessionIndex, sessionID, animal, probeDirDeg);

conditionUse = isfinite(trialTable.prefCohNoisePC) & ...
  isfinite(trialTable.probeCohNoisePC);
conditionTrials = trialTable(conditionUse,:);
[conditionGroups, probeDirDeg, prefCohNoisePC, probeCohNoisePC, ...
  stepDeltaCohPC] = findgroups(conditionTrials.probeDirDeg, ...
  conditionTrials.prefCohNoisePC, conditionTrials.probeCohNoisePC, ...
  conditionTrials.stepDeltaCohPC);
conditionSummary = groupedSummary(conditionTrials, conditionGroups, ...
  {'probeDirDeg','prefCohNoisePC','probeCohNoisePC','stepDeltaCohPC'}, ...
  probeDirDeg, prefCohNoisePC, probeCohNoisePC, stepDeltaCohPC);

offsetSummary = sortrows(offsetSummary, 'probeDirDeg');
sessionOffsetSummary = sortrows(sessionOffsetSummary, ...
  {'probeDirDeg','sessionIndex'});
conditionSummary = sortrows(conditionSummary, ...
  {'probeDirDeg','prefCohNoisePC','probeCohNoisePC','stepDeltaCohPC'});

O = struct();
O.createdAt = datetime('now');
O.createdBy = mfilename;
O.predictorDefinition = 'signed mean physical coherence over step';
O.probeScaling = 'single-candidate probe predictor; effective yoked sum divided by nYokedProbeStreams';
O.changeMarginDefinition = 'changeProbeCandidate - (stepDeltaCohPC + changePreferred)';
O.noChangeMarginDefinition = 'noChangeProbeCandidate - noChangePreferred';
O.crossingDefinition = 'margin > 0; exact equality reported separately';
O.trialTable = trialTable;
O.offsetSummary = offsetSummary;
O.sessionOffsetSummary = sessionOffsetSummary;
O.conditionSummary = conditionSummary;
end

%% ------------------------------------------------------------------------
function x = stepMean(noiseByFrameTrial, stepFrames)

x = mean(double(noiseByFrameTrial(stepFrames,:)), 1, 'omitnan')';
end

%% ------------------------------------------------------------------------
function S = groupedSummary(T, G, keyNames, varargin)

nGroups = max(G);
if isempty(nGroups) || ~isfinite(nGroups)
  S = table();
  return
end

if isempty(varargin)
  keyValues = cell(1,numel(keyNames));
  for k = 1:numel(keyNames)
    keyValues{k} = splitapply(@(x) x(1), T.(keyNames{k}), G);
  end
else
  keyValues = varargin;
end

nTrials = zeros(nGroups,1);
nSessions = zeros(nGroups,1);
meanChangePreferred = nan(nGroups,1);
sdChangePreferred = nan(nGroups,1);
meanChangeProbe = nan(nGroups,1);
sdChangeProbe = nan(nGroups,1);
corrChangePrefProbe = nan(nGroups,1);
meanNoChangePreferred = nan(nGroups,1);
sdNoChangePreferred = nan(nGroups,1);
meanNoChangeProbe = nan(nGroups,1);
sdNoChangeProbe = nan(nGroups,1);
corrNoChangePrefProbe = nan(nGroups,1);
meanStep = nan(nGroups,1);
sdStep = nan(nGroups,1);
medianStep = nan(nGroups,1);
q25Step = nan(nGroups,1);
q75Step = nan(nGroups,1);
q95Step = nan(nGroups,1);
minStep = nan(nGroups,1);
maxStep = nan(nGroups,1);
meanStepOverSDChangePreferred = nan(nGroups,1);
meanStepOverSDChangeProbe = nan(nGroups,1);
meanChangeMargin = nan(nGroups,1);
sdChangeMargin = nan(nGroups,1);
medianChangeMargin = nan(nGroups,1);
q90ChangeMargin = nan(nGroups,1);
q95ChangeMargin = nan(nGroups,1);
q975ChangeMargin = nan(nGroups,1);
q99ChangeMargin = nan(nGroups,1);
maxChangeMargin = nan(nGroups,1);
nChangeCrossed = zeros(nGroups,1);
pChangeCrossed = nan(nGroups,1);
pChangeCrossedCI95Low = nan(nGroups,1);
pChangeCrossedCI95High = nan(nGroups,1);
nChangeEqual = zeros(nGroups,1);
meanNoChangeMargin = nan(nGroups,1);
sdNoChangeMargin = nan(nGroups,1);
medianNoChangeMargin = nan(nGroups,1);
q95NoChangeMargin = nan(nGroups,1);
q99NoChangeMargin = nan(nGroups,1);
maxNoChangeMargin = nan(nGroups,1);
nNoChangeCrossed = zeros(nGroups,1);
pNoChangeCrossed = nan(nGroups,1);
pNoChangeCrossedCI95Low = nan(nGroups,1);
pNoChangeCrossedCI95High = nan(nGroups,1);
nNoChangeEqual = zeros(nGroups,1);

for i = 1:nGroups
  X = T(G == i,:);
  nTrials(i) = height(X);
  nSessions(i) = numel(unique(X.sessionIndex));

  [meanChangePreferred(i), sdChangePreferred(i)] = meanSD(X.changePreferred);
  [meanChangeProbe(i), sdChangeProbe(i)] = meanSD(X.changeProbeCandidate);
  corrChangePrefProbe(i) = safeCorrelation( ...
    X.changePreferred, X.changeProbeCandidate);
  [meanNoChangePreferred(i), sdNoChangePreferred(i)] = meanSD(X.noChangePreferred);
  [meanNoChangeProbe(i), sdNoChangeProbe(i)] = meanSD(X.noChangeProbeCandidate);
  corrNoChangePrefProbe(i) = safeCorrelation( ...
    X.noChangePreferred, X.noChangeProbeCandidate);

  step = double(X.stepDeltaCohPC);
  [meanStep(i), sdStep(i)] = meanSD(step);
  qStep = percentileValues(step, [25 50 75 95]);
  q25Step(i) = qStep(1);
  medianStep(i) = qStep(2);
  q75Step(i) = qStep(3);
  q95Step(i) = qStep(4);
  minStep(i) = min(step, [], 'omitnan');
  maxStep(i) = max(step, [], 'omitnan');
  meanStepOverSDChangePreferred(i) = safeRatio( ...
    meanStep(i), sdChangePreferred(i));
  meanStepOverSDChangeProbe(i) = safeRatio( ...
    meanStep(i), sdChangeProbe(i));

  changeMargin = double(X.changeMargin);
  [meanChangeMargin(i), sdChangeMargin(i)] = meanSD(changeMargin);
  q = percentileValues(changeMargin, [50 90 95 97.5 99 100]);
  medianChangeMargin(i) = q(1);
  q90ChangeMargin(i) = q(2);
  q95ChangeMargin(i) = q(3);
  q975ChangeMargin(i) = q(4);
  q99ChangeMargin(i) = q(5);
  maxChangeMargin(i) = q(6);
  nChangeCrossed(i) = sum(X.changeCrossed);
  nChangeEqual(i) = sum(X.changeEqual);
  pChangeCrossed(i) = nChangeCrossed(i) / nTrials(i);
  ci = wilsonInterval(nChangeCrossed(i), nTrials(i));
  pChangeCrossedCI95Low(i) = ci(1);
  pChangeCrossedCI95High(i) = ci(2);

  noChangeMargin = double(X.noChangeMargin);
  [meanNoChangeMargin(i), sdNoChangeMargin(i)] = meanSD(noChangeMargin);
  q = percentileValues(noChangeMargin, [50 95 99 100]);
  medianNoChangeMargin(i) = q(1);
  q95NoChangeMargin(i) = q(2);
  q99NoChangeMargin(i) = q(3);
  maxNoChangeMargin(i) = q(4);
  nNoChangeCrossed(i) = sum(X.noChangeCrossed);
  nNoChangeEqual(i) = sum(X.noChangeEqual);
  pNoChangeCrossed(i) = nNoChangeCrossed(i) / nTrials(i);
  ci = wilsonInterval(nNoChangeCrossed(i), nTrials(i));
  pNoChangeCrossedCI95Low(i) = ci(1);
  pNoChangeCrossedCI95High(i) = ci(2);
end

S = table(keyValues{:}, nTrials, nSessions, ...
  meanChangePreferred, sdChangePreferred, meanChangeProbe, sdChangeProbe, ...
  corrChangePrefProbe, meanNoChangePreferred, sdNoChangePreferred, ...
  meanNoChangeProbe, sdNoChangeProbe, corrNoChangePrefProbe, ...
  meanStep, sdStep, medianStep, q25Step, q75Step, q95Step, minStep, maxStep, ...
  meanStepOverSDChangePreferred, meanStepOverSDChangeProbe, ...
  meanChangeMargin, sdChangeMargin, medianChangeMargin, ...
  q90ChangeMargin, q95ChangeMargin, q975ChangeMargin, q99ChangeMargin, ...
  maxChangeMargin, nChangeCrossed, pChangeCrossed, ...
  pChangeCrossedCI95Low, pChangeCrossedCI95High, nChangeEqual, ...
  meanNoChangeMargin, sdNoChangeMargin, medianNoChangeMargin, ...
  q95NoChangeMargin, q99NoChangeMargin, maxNoChangeMargin, ...
  nNoChangeCrossed, pNoChangeCrossed, ...
  pNoChangeCrossedCI95Low, pNoChangeCrossedCI95High, nNoChangeEqual, ...
  'VariableNames', [keyNames, { ...
  'nTrials','nSessions', ...
  'meanChangePreferred','sdChangePreferred','meanChangeProbe','sdChangeProbe', ...
  'corrChangePrefProbe','meanNoChangePreferred','sdNoChangePreferred', ...
  'meanNoChangeProbe','sdNoChangeProbe','corrNoChangePrefProbe', ...
  'meanStep','sdStep','medianStep','q25Step','q75Step','q95Step','minStep','maxStep', ...
  'meanStepOverSDChangePreferred','meanStepOverSDChangeProbe', ...
  'meanChangeMargin','sdChangeMargin','medianChangeMargin', ...
  'q90ChangeMargin','q95ChangeMargin','q975ChangeMargin','q99ChangeMargin', ...
  'maxChangeMargin','nChangeCrossed','pChangeCrossed', ...
  'pChangeCrossedCI95Low','pChangeCrossedCI95High','nChangeEqual', ...
  'meanNoChangeMargin','sdNoChangeMargin','medianNoChangeMargin', ...
  'q95NoChangeMargin','q99NoChangeMargin','maxNoChangeMargin', ...
  'nNoChangeCrossed','pNoChangeCrossed', ...
  'pNoChangeCrossedCI95Low','pNoChangeCrossedCI95High','nNoChangeEqual'}]);
end

%% ------------------------------------------------------------------------
function [m, s] = meanSD(x)

x = double(x(:));
x = x(isfinite(x));
if isempty(x)
  m = nan;
  s = nan;
else
  m = mean(x);
  s = std(x);
end
end

%% ------------------------------------------------------------------------
function r = safeCorrelation(x, y)

x = double(x(:));
y = double(y(:));
use = isfinite(x) & isfinite(y);
x = x(use);
y = y(use);
if numel(x) < 2 || std(x) == 0 || std(y) == 0
  r = nan;
else
  C = corrcoef(x,y);
  r = C(1,2);
end
end

%% ------------------------------------------------------------------------
function r = safeRatio(a, b)

if isfinite(a) && isfinite(b) && b > 0
  r = a / b;
else
  r = nan;
end
end

%% ------------------------------------------------------------------------
function q = percentileValues(x, percentages)

x = sort(double(x(:)));
x = x(isfinite(x));
if isempty(x)
  q = nan(size(percentages));
  return
end

% Linear interpolation matching the usual empirical quantile convention.
n = numel(x);
p = double(percentages(:)') / 100;
h = 1 + (n - 1) .* p;
lo = floor(h);
hi = ceil(h);
w = h - lo;
q = (1-w) .* x(lo)' + w .* x(hi)';
end

%% ------------------------------------------------------------------------
function ci = wilsonInterval(k, n)

if n <= 0
  ci = [nan nan];
  return
end
z = 1.95996398454005;
p = k / n;
denom = 1 + z^2/n;
center = (p + z^2/(2*n)) / denom;
halfWidth = z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / denom;
ci = [max(0,center-halfWidth), min(1,center+halfWidth)];
end
