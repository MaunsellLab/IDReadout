function psych = fitIDQInitialSessionThresholds(trialTable, betaWeibull, lapse, targetPerformance)
% fitIDQInitialSessionThresholds
%
% Fit one fixed-beta, fixed-lapse Weibull threshold per session.
%
% Inputs:
%   trialTable          across-session trial table with sessionIndex
%   betaWeibull         fixed Weibull shape
%   lapse               fixed lapse
%   targetPerformance   threshold-defining performance, e.g. 0.75
%
% Output:
%   psych struct with per-session fit table and aligned coherence values

sessionIDs = unique(trialTable.sessionIndex);
nSessions = numel(sessionIDs);

sessionIndex = nan(nSessions, 1);
sessionName = strings(nSessions, 1);
nTrials = nan(nSessions, 1);
nCorrect = nan(nSessions, 1);
meanCorrect = nan(nSessions, 1);
noisyStepCoh = nan(nSessions, 1);

threshold = nan(nSessions, 1);
alpha = nan(nSessions, 1);
negLogLikelihood = nan(nSessions, 1);
exitflag = nan(nSessions, 1);

alignedCoh = nan(height(trialTable), 1);
noisyStepAlignedCoh = nan(height(trialTable), 1);

for i = 1:nSessions

  s = sessionIDs(i);
  idx = trialTable.sessionIndex == s;

  fit = fitIDQSessionWeibullFixedBeta( ...
    trialTable.stepCoh(idx), ...
    trialTable.correct(idx), ...
    betaWeibull, ...
    lapse, ...
    targetPerformance);

  sessionIndex(i) = s;
  sessionName(i) = trialTable.sessionName(find(idx, 1, 'first'));
  nTrials(i) = fit.nTrials;
  nCorrect(i) = fit.nCorrect;
  meanCorrect(i) = fit.nCorrect / fit.nTrials;

  noisyCohThis = unique(trialTable.noisyStepCoh(idx));
  if numel(noisyCohThis) ~= 1
    error('fitIDQInitialSessionThresholds:NoisyCohMismatch', ...
      'Expected one noisyStepCoh for session %d.', s);
  end
  noisyStepCoh(i) = noisyCohThis;

  threshold(i) = fit.threshold;
  alpha(i) = fit.alpha;
  negLogLikelihood(i) = fit.negLogLikelihood;
  exitflag(i) = fit.exitflag;

  alignedCoh(idx) = trialTable.stepCoh(idx) ./ fit.threshold;
  noisyStepAlignedCoh(idx) = trialTable.noisyStepCoh(idx) ./ fit.threshold;

end

sessionFits = table( ...
  sessionIndex, ...
  sessionName, ...
  nTrials, ...
  nCorrect, ...
  meanCorrect, ...
  noisyStepCoh, ...
  threshold, ...
  alpha, ...
  negLogLikelihood, ...
  exitflag);

psych = struct();

psych.betaWeibull = betaWeibull;
psych.lapse = lapse;
psych.targetPerformance = targetPerformance;
psych.sessionFits = sessionFits;
psych.alignedCoh = alignedCoh;
psych.noisyStepAlignedCoh = noisyStepAlignedCoh;

end