function [kernelLogOdds, kernelSE, P0] = logisticPsychKernel1D(cohMat, trialOutcome)
% logisticPsychKernel1D  Robust logistic-regression kernel, 1D per time bin.
%
%   [kernelLogOdds, kernelSE, P0] = logisticPsychKernel1D(cohMat, trialOutcome)
%
%   INPUTS:
%     cohMat       : m x nTrials matrix of coherence noise values
%                    (rows = time bins, columns = trials).
%     trialOutcome : nTrials x 1 vector, 0 = correct, 1 = wrong.
%
%   OUTPUTS:
%     kernelLogOdds : m x 1 vector of regression weights beta(t).
%                     Units: change in log-odds(error) per unit coherence
%                     at time t.
%
%     kernelSE      : m x 1 vector of standard errors of beta(t).
%
%     P0            : scalar baseline error probability in this dataset.
%
%   Near P0, a 1-unit increase in coherence at time t changes error
%   probability by approximately:
%
%       dP_error/ds ≈ P0 * (1 - P0) * kernelLogOdds(t)

    [m, nTrials] = size(cohMat);
    if nTrials ~= numel(trialOutcome)
        error('cohMat has %d columns but trialOutcome has %d elements.', ...
              nTrials, numel(trialOutcome));
    end
    y = 1 - trialOutcome(:);  % (now 1 = correct, 0 = error)
    P0 = mean(y);
    if P0 == 0 || P0 == 1
        error('logisticPsychKernel1D:AllSameOutcome', ...
              'All outcomes are %d; need both correct and wrong trials.', P0);
    end

    kernelLogOdds = nan(m, 1);
    kernelSE      = nan(m, 1);

    for t = 1:m
        x = cohMat(t, :).';  % nTrials x 1 predictor for time bin t

        % If this time bin has no variance across trials, kernel is undefined
        if all(x == x(1))
            kernelLogOdds(t) = 0;
            kernelSE(t)      = NaN;
            continue;
        end

        u = unique(x);

        % ---- Special handling for binary predictor (common in your case) ----
        if numel(u) == 2
            % Count correct/error at each coherence level
            counts = zeros(2, 2); % rows: u(1), u(2); cols: correct(1), error(2)
            for j = 1:2
                mask = (x == u(j));
                n_j   = sum(mask);
                err_j = sum(y(mask) == 1);
                cor_j = n_j - err_j;
                counts(j, 1) = cor_j;
                counts(j, 2) = err_j;
            end

            % If any cell is zero, we have (quasi-)separation → avoid glmfit
            if any(counts(:) == 0)
                % Jeffreys prior (0.5 pseudo-count) to avoid 0 or 1 probabilities
                logitVals = zeros(2,1);
                varLogit  = zeros(2,1);
                for j = 1:2
                    n_j   = sum(x == u(j));
                    err_j = sum(y(x == u(j)) == 1);
                    p_err_j = (err_j + 0.5) / (n_j + 1); % Beta(0.5,0.5) posterior mean
                    logitVals(j) = log(p_err_j / (1 - p_err_j));
                    % Approximate var(logit(p)) ≈ 1 / (n * p * (1-p))
                    varLogit(j) = 1 / max(n_j * p_err_j * (1 - p_err_j), eps);
                end

                dx = u(2) - u(1);
                beta_t = (logitVals(2) - logitVals(1)) / dx;
                se_t   = sqrt(sum(varLogit)) / abs(dx);

                kernelLogOdds(t) = beta_t;
                kernelSE(t)      = se_t;
                continue;
            end
        end

        % ---- Default: regular 1D logistic regression ----
        try
            [b, ~, stats] = glmfit(x, y, 'binomial', 'link', 'logit');
            kernelLogOdds(t) = b(2);       % slope
            kernelSE(t)      = stats.se(2);
        catch ME
            warning('Time bin %d: glmfit failed (%s). Setting kernel to 0.', ...
                    t, ME.message);
            kernelLogOdds(t) = 0;
            kernelSE(t)      = NaN;
        end
    end
end
