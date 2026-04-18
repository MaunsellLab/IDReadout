function scaleSummary = computeScaleEstimators(summary)
% computeScaleEstimators
%
% Compute multiple scale estimators from per-session summary table.
%
% Inputs:
%   summary  (table from collectKernelSessionSummary)
%
% Outputs:
%   scaleSummary (table with one row per sideType x stepType)
%
% Estimators:
%   1) mean(sessionScale)
%   2) pref-energy-weighted mean
%   3) aggregate numerator/denominator (sum crossDot / sum prefEnergy)
%   4) trial-weighted mean (optional but useful)
%
% Also reports:
%   - number of sessions
%   - total trials
%   - mean cosine similarity

sideTypes = unique(summary.sideType);
stepTypes = unique(summary.stepType);

rows = struct( ...
  'sideType', {}, ...
  'stepType', {}, ...
  'sideLabel', {}, ...
  'stepLabel', {}, ...
  'nSessions', {}, ...
  'totalTrials', {}, ...
  'meanScale', {}, ...
  'prefE_weighted', {}, ...
  'numden', {}, ...
  'trialWeighted', {}, ...
  'meanCosSim', {} );

rowCounter = 0;

for sideType = sideTypes'
  for stepType = stepTypes'

    idx = summary.sideType == sideType & summary.stepType == stepType;

    if ~any(idx)
      continue;
    end

    s = summary(idx, :);

    sessScale = s.sessionScale;
    prefE     = s.prefEnergy;
    crossDot  = s.crossDot;
    nTrials   = s.nTrialsStep;
    cosSim    = s.cosSim;

    % Remove NaNs consistently
    valid = isfinite(sessScale) & isfinite(prefE) & prefE > 0;

    sessScale = sessScale(valid);
    prefE     = prefE(valid);
    crossDot  = crossDot(valid);
    nTrials   = nTrials(valid);
    cosSim    = cosSim(valid);

    if isempty(sessScale)
      continue;
    end

    % --- Estimators ---

    % 1. Mean of session scales
    meanScale = mean(sessScale);

    % 2. Pref-energy-weighted mean
    prefE_weighted = sum(prefE .* sessScale) / sum(prefE);

    % 3. Aggregate numerator/denominator
    numden = sum(crossDot) / sum(prefE);

    % 4. Trial-weighted mean (optional)
    trialWeighted = sum(nTrials .* sessScale) / sum(nTrials);

    % --- Diagnostics ---
    nSessions = numel(sessScale);
    totalTrials = sum(nTrials);
    meanCosSim = mean(cosSim, 'omitnan');

    rowCounter = rowCounter + 1;

    rows(rowCounter).sideType = sideType;
    rows(rowCounter).stepType = stepType;

    rows(rowCounter).sideLabel = s.sideLabel(1);
    rows(rowCounter).stepLabel = s.stepLabel(1);

    rows(rowCounter).nSessions = nSessions;
    rows(rowCounter).totalTrials = totalTrials;

    rows(rowCounter).meanScale = meanScale;
    rows(rowCounter).prefE_weighted = prefE_weighted;
    rows(rowCounter).numden = numden;
    rows(rowCounter).trialWeighted = trialWeighted;

    rows(rowCounter).meanCosSim = meanCosSim;
  end
end

scaleSummary = struct2table(rows);

% Sort for readability
scaleSummary = sortrows(scaleSummary, {'sideType', 'stepType'});

% --- Print compact summary ---
fprintf('\n=== Scale Estimator Summary ===\n');

for i = 1:height(scaleSummary)
  r = scaleSummary(i,:);
  fprintf('%-8s %-3s | n=%2d | mean=% .3f | prefE=% .3f | numden=% .3f | trial=% .3f | cos=%.2f\n', ...
    r.sideLabel{1}, r.stepLabel{1}, ...
    r.nSessions, ...
    r.meanScale, ...
    r.prefE_weighted, ...
    r.numden, ...
    r.trialWeighted, ...
    r.meanCosSim);
end

fprintf('================================\n\n');

end