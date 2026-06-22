function gainFit = fitIDQNoiseGain(trialTable, predictorName, sessionFits, alignedWeibull, targetPerformance)
% fitIDQNoiseGain
%
% Fit one shared gain for rectangular-weighted signed noise.
%
% Model:
%   effectiveCoh = stepCoh + gain * rectNoisePredictor
%   P(correct)  = Weibull(effectiveCoh; threshold_session, betaWeibull, lapse)
%
% Flat-readout prediction:
%   gain = 1
%
% Inputs:
%   trialTable        across-session trial table
%   sessionFits       pass-2 session fit table, with sessionIndex, threshold
%   alignedWeibull    pass-2 across aligned Weibull fit
%   targetPerformance threshold-defining performance, e.g. 0.75

idx = trialTable.hasStepNoise;

sessionIndex = trialTable.sessionIndex(idx);
stepCoh = double(trialTable.stepCoh(idx));
xNoise = double(trialTable.(predictorName)(idx));
correct = double(trialTable.correct(idx));

valid = isfinite(sessionIndex) & isfinite(stepCoh) & ...
        isfinite(xNoise) & isfinite(correct);

sessionIndex = sessionIndex(valid);
stepCoh = stepCoh(valid);
xNoise = xNoise(valid);
correct = correct(valid);

betaWeibull = alignedWeibull.betaWeibull;
lapse = alignedWeibull.lapse;

sessionThreshold = nan(size(sessionIndex));
for i = 1:height(sessionFits)
    idxSession = sessionIndex == sessionFits.sessionIndex(i);
    sessionThreshold(idxSession) = sessionFits.threshold(i);
end

if any(~isfinite(sessionThreshold))
    error('fitIDQNoiseGainRect:MissingSessionThreshold', ...
        'Could not assign threshold to all noisy trials.');
end

sessionAlpha = idqWeibullAlphaForThreshold( ...
    targetPerformance, sessionThreshold, betaWeibull, lapse);

gainFit = struct();
gainFit.predictor = predictorName;
gainFit.model = sprintf('Weibull(stepCoh + gain * %s)', predictorName);
gainFit.flatPrediction = 1;

gainFit.gain = NaN;
gainFit.SE = NaN;
gainFit.CI95 = [NaN NaN];
gainFit.negLogLikelihood = NaN;
gainFit.exitflag = NaN;
gainFit.nTrials = numel(correct);
gainFit.nCorrect = sum(correct == 1);
gainFit.betaWeibull = betaWeibull;
gainFit.lapse = lapse;
gainFit.targetPerformance = targetPerformance;

objective = @(gain) negLogLikelihoodGain(gain, stepCoh, xNoise, correct, sessionAlpha, betaWeibull, lapse);

% Use a generous but finite bracket. Noise and signal are in the same units,
% and gain near 1 is the main prediction.
lb = -5;
ub = 5;

opts = optimset( ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'MaxFunEvals', 2000, ...
    'MaxIter', 2000);

[gainHat, nll, exitflag] = fminbnd(objective, lb, ub, opts);

gainFit.gain = gainHat;
gainFit.negLogLikelihood = nll;
gainFit.exitflag = exitflag;

% Approximate SE from finite-difference second derivative.
h = 1e-4;
f0 = objective(gainHat);
fp = objective(gainHat + h);
fm = objective(gainHat - h);
secondDeriv = (fp - 2 * f0 + fm) / h^2;

gainFit.secondDerivative = secondDeriv;

if isfinite(secondDeriv) && secondDeriv > 0
    gainFit.SE = sqrt(1 / secondDeriv);
    gainFit.CI95 = gainHat + [-1 1] * 1.96 * gainFit.SE;
end

% Useful descriptive diagnostic.
gainFit.meanXNoise = mean(xNoise, 'omitnan');
gainFit.stdXNoise = std(xNoise, 'omitnan');
gainFit.meanXNoiseCorrect = mean(xNoise(correct == 1), 'omitnan');
gainFit.meanXNoiseError = mean(xNoise(correct == 0), 'omitnan');
gainFit.meanDiffCorrectMinusError = ...
    gainFit.meanXNoiseCorrect - gainFit.meanXNoiseError;
gainFit.nEffectiveCohClipped = sum(stepCoh + gainHat .* xNoise < 0);
gainFit.zVsFlat = (gainFit.gain - gainFit.flatPrediction) / gainFit.SE;
gainFit.zVsZero = gainFit.gain / gainFit.SE;
gainFit.pVsFlat = 2 * normcdf(-abs(gainFit.zVsFlat));
gainFit.pVsZero = 2 * normcdf(-abs(gainFit.zVsZero));

end

%% ------------------------------------------------------------------------
function nll = negLogLikelihoodGain(gain, stepCoh, xNoise, correct, ...
    sessionAlpha, betaWeibull, lapse)

effectiveCoh = stepCoh + gain .* xNoise;

% Pedestal/noise implementation should keep relevant effective coherences
% nonnegative, but guard against numerical excursions during wide gain search.
effectiveCoh = max(effectiveCoh, 0);

p = idqWeibullP(effectiveCoh, sessionAlpha, betaWeibull, lapse);
p = min(max(p, eps), 1 - eps);

nll = -sum(correct .* log(p) + (1 - correct) .* log(1 - p));

end