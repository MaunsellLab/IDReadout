function [allocationTable, sessionAllocationTable] = summarizeIDLAllocationSupport(input)
% summarizeIDLAllocationSupport  Held-out model support by signal allocation.
%
% [allocationTable, sessionAllocationTable] = ...
%   summarizeIDLAllocationSupport(result)
% [allocationTable, sessionAllocationTable] = ...
%   summarizeIDLAllocationSupport('Data/AcrossSessions/ModelRecovery.mat')
%
% Uses the per-trial predictions already stored by
% fitIDLModelRecoveryAcrossSessions. The 100:0 calibration trials are
% excluded. Positive observed delta LL favors Max. Expected KL columns
% describe discrimination at the ACTUAL sampled steps and fitted
% predictions; they are not a full-pipeline power calculation.

if ischar(input) || (isstring(input) && isscalar(input))
  saved = load(char(input), 'result');
  result = saved.result;
else
  result = input;
end

P = result.predictionTable;
S = result.sessionTable;
sessionIDs = string(S.sourceSessionID);
[found, sessionIndex] = ismember(string(P.sourceSessionID), sessionIDs);
assert(all(found), 'summarizeIDLAllocationSupport:UnknownSession', ...
  'The prediction table contains a session absent from the session table.');
allocation = double(P.unsignedAllocPC);
calibration = logical(P.isCalibrationTrial);
allocations = sort(unique(allocation(~calibration)), 'descend');
assert(~isempty(allocations), ...
  'summarizeIDLAllocationSupport:NoDividedTrials', ...
  'There are no divided-allocation prediction trials.');

observed = double(P.correct);
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
klMaxColumn = zeros(0,1);
klAltColumn = zeros(0,1);

for model = 1:2
  pAlt = alternatives(:,model);
  trialDelta = observed .* log(pMax ./ pAlt) + ...
    (1 - observed) .* log((1 - pMax) ./ (1 - pAlt));
  klMax = pMax .* log(pMax ./ pAlt) + ...
    (1 - pMax) .* log((1 - pMax) ./ (1 - pAlt));
  klAlt = pAlt .* log(pAlt ./ pMax) + ...
    (1 - pAlt) .* log((1 - pAlt) ./ (1 - pMax));

  for k = 1:height(S)
    for j = 1:numel(allocations)
      rows = ~calibration & sessionIndex == k & ...
        abs(allocation - allocations(j)) < 1e-6;
      n = sum(rows);
      if n == 0
        continue
      end
      sessionIDColumn(end+1,1) = sessionIDs(k);
      animalColumn(end+1,1) = string(S.animal(k));
      comparisonColumn(end+1,1) = comparisonNames(model);
      allocationColumn(end+1,1) = allocations(j);
      nTrialsColumn(end+1,1) = n;
      deltaLLColumn(end+1,1) = sum(trialDelta(rows));
      per1000Column(end+1,1) = 1000 * deltaLLColumn(end) / n;
      klMaxColumn(end+1,1) = 1000 * mean(klMax(rows));
      klAltColumn(end+1,1) = 1000 * mean(klAlt(rows));
    end
  end
end

sessionAllocationTable = table(sessionIDColumn, animalColumn, ...
  comparisonColumn, allocationColumn, nTrialsColumn, deltaLLColumn, ...
  per1000Column, klMaxColumn, klAltColumn, 'VariableNames', ...
  {'sourceSessionID','animal','comparison','allocationPC','nTrials', ...
  'deltaLL','per1000','expectedKLPer1000UnderMax', ...
  'expectedKLPer1000UnderAlternative'});

groups = unique(string(S.animal), 'stable');
if numel(groups) > 1
  groups(end+1,1) = "Combined";
end
groupColumn = strings(0,1);
comparisonSummary = strings(0,1);
allocationSummary = zeros(0,1);
nSessions = zeros(0,1);
nTrials = zeros(0,1);
nPositive = zeros(0,1);
totalDeltaLL = zeros(0,1);
pooledPer1000 = zeros(0,1);
meanSessionPer1000 = zeros(0,1);
sdSessionPer1000 = zeros(0,1);
expectedKLPer1000UnderMax = zeros(0,1);
expectedKLPer1000UnderAlternative = zeros(0,1);

for g = 1:numel(groups)
  for model = 1:2
    for j = 1:numel(allocations)
      rows = sessionAllocationTable.comparison == comparisonNames(model) & ...
        abs(sessionAllocationTable.allocationPC - allocations(j)) < 1e-6;
      if groups(g) ~= "Combined"
        rows = rows & sessionAllocationTable.animal == groups(g);
      end
      A = sessionAllocationTable(rows,:);
      if isempty(A)
        continue
      end
      groupColumn(end+1,1) = groups(g);
      comparisonSummary(end+1,1) = comparisonNames(model);
      allocationSummary(end+1,1) = allocations(j);
      nSessions(end+1,1) = height(A);
      nTrials(end+1,1) = sum(A.nTrials);
      nPositive(end+1,1) = sum(A.deltaLL > 0);
      totalDeltaLL(end+1,1) = sum(A.deltaLL);
      pooledPer1000(end+1,1) = 1000 * totalDeltaLL(end) / nTrials(end);
      meanSessionPer1000(end+1,1) = mean(A.per1000);
      sdSessionPer1000(end+1,1) = std(A.per1000);
      expectedKLPer1000UnderMax(end+1,1) = ...
        sum(A.nTrials .* A.expectedKLPer1000UnderMax) / nTrials(end);
      expectedKLPer1000UnderAlternative(end+1,1) = ...
        sum(A.nTrials .* A.expectedKLPer1000UnderAlternative) / nTrials(end);
    end
  end
end

allocationTable = table(groupColumn, comparisonSummary, ...
  allocationSummary, nSessions, nTrials, nPositive, totalDeltaLL, ...
  pooledPer1000, meanSessionPer1000, sdSessionPer1000, ...
  expectedKLPer1000UnderMax, expectedKLPer1000UnderAlternative, ...
  'VariableNames', {'group','comparison','allocationPC','nSessions', ...
  'nTrials','nPositive','totalDeltaLL','pooledPer1000', ...
  'meanSessionPer1000','sdSessionPer1000', ...
  'expectedKLPer1000UnderMax','expectedKLPer1000UnderAlternative'});

if nargout == 0
  disp(allocationTable);
end
end

function p = clippedProbability(p)
p = min(max(double(p(:)), 1e-12), 1 - 1e-12);
end
