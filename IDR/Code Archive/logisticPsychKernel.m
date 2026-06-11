function [kernelLogOdds, kernelSE, P0] = logisticPsychKernel(cohMat, trialOutcome)
% logisticPsychKernel  Logistic-regression-based psychophysical kernel.
%
%   [kernelLogOdds, kernelSE, P0] = logisticPsychKernel(cohMat, trialOutcome)
%
%   INPUTS:
%     cohMat       : m x nTrials matrix of coherence noise values
%                    (rows = time bins, columns = trials).
%                    Units are whatever you used for coherence (e.g., %).
%
%     trialOutcome : nTrials x 1 vector of trial outcomes
%                    0 = correct, 1 = wrong.
%
%   OUTPUTS:
%     kernelLogOdds : m x 1 vector of regression weights beta(t).
%                     Units: change in log-odds(error) per unit coherence
%                     at time t.
%
%     kernelSE      : m x 1 vector of standard errors of beta(t).
%
%     P0            : scalar, baseline error probability in this dataset.
%
%   Near the baseline error rate P0, a 1-unit increase in coherence at
%   time t changes error probability by approximately:
%
%       dP_error/ds  ≈  P0 * (1 - P0) * kernelLogOdds(t)
%
%   so you can convert to an approximate probability change per unit
%   coherence if you wish.

    % Basic checks
    [m, nTrials] = size(cohMat);
    if nTrials ~= numel(trialOutcome)
        error('cohMat has %d columns but trialOutcome has %d elements.', ...
              nTrials, numel(trialOutcome));
    end

    % Design matrix: trials x predictors
    X = cohMat.';             % nTrials x m
    y = trialOutcome(:);      % nTrials x 1 (0/1)

    % Baseline error rate for interpretation
    P0 = mean(y);
    if P0 == 0 || P0 == 1
        error(['logisticPsychKernel:AllSameOutcome', ...
               ' — All outcomes are %d. Need both correct and wrong trials.'], P0);
    end

    % Fit logistic regression:
    %   logit( P(error) ) = beta0 + sum_t beta_t * s_t
    %
    % 'binomial' + default 'logit' link gives standard logistic regression.
    glm = fitglm(X, y, 'Distribution', 'binomial', 'Link', 'logit');

    % Coefficients: first is intercept, rest are the kernel over time
    beta    = glm.Coefficients.Estimate;  % (m+1) x 1
    betaSE  = glm.Coefficients.SE;        % (m+1) x 1

    kernelLogOdds = beta(2:end);          % m x 1
    kernelSE      = betaSE(2:end);        % m x 1

    % (Optional sanity check)
    if numel(kernelLogOdds) ~= m
        error('Unexpected number of coefficients (%d), expected %d.', ...
              numel(kernelLogOdds), m);
    end
end
