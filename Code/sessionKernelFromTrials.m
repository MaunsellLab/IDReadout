function [kernels, kVars, kStats] = sessionKernelFromTrials(sessionData, trialIdx)
%SESSIONKERNELFROMTRIALS Recompute session kernels from trial-level noise matrices.
%
%   [kernels, kVars, kStats] = sessionKernelFromTrials(sessionData, trialIdx)
%
% Inputs
%   sessionData struct with fields:
%       header          : session header struct
%       prefNoise       : [nFrames x nTrials] pref-direction noise, in % coherence
%       probeNoise      : [nFrames x nTrials] probe-direction noise, in % coherence
%       trialOutcomes   : [1 x nTrials] or [nTrials x 1], 0=correct, 1=wrong
%       changeSides     : [1 x nTrials] or [nTrials x 1], 0=RF changed, 1=Opp changed
%       changeIndices   : [1 x nTrials] or [nTrials x 1], 1=DEC, 2=INC
%
%   trialIdx            : vector of trial indices to use. May contain repeats
%                         for bootstrap resampling. If omitted or empty, all trials
%                         are used.
%
% Outputs
%   kernels             : [5 x 2 x nFrames] mean(correct)-mean(wrong) kernels
%                         dim1 = [diff change noChange RF Opp]
%                         dim2 = [DEC INC]
%   kVars               : [5 x 2 x nFrames] approximate variance of kernel estimate
%   kStats              : struct with counts and bookkeeping
%
% Notes
%   - This function assumes invalid trials have already been removed upstream.
%   - Probe noise is stored raw and projected onto the drift/pref axis here.
%   - If previous kernel code used a different sign convention for DEC,
%     change only the helper 'signedNoiseMatrices' below.

    if nargin < 2 || isempty(trialIdx)
        trialIdx = 1:size(sessionData.prefNoise, 2);
    end

    % Pull selected trials. trialIdx may contain repeats (bootstrap).
    prefNoise    = sessionData.prefNoise(:, trialIdx);
    probeNoise   = sessionData.probeNoise(:, trialIdx);
    trialOutcomes = sessionData.trialOutcomes(trialIdx);
    changeSides   = sessionData.changeSides(trialIdx);
    changeIndices = sessionData.changeIndices(trialIdx);

    % Standardize row vectors for logical indexing.
    trialOutcomes = reshape(trialOutcomes, 1, []);
    changeSides   = reshape(changeSides,   1, []);
    changeIndices = reshape(changeIndices, 1, []);

    nFrames = size(prefNoise, 1);

    % ---------------------------------------------------------------------
    % Project probe noise onto drift/pref axis.
    %
    % Assumption:
    %   header contains the probe offset relative to pref/drift in degrees.
    %   Replace 'probeDirOffsetDeg' below with the actual field path if needed.
    % ---------------------------------------------------------------------
    probeOffsetDeg = sessionData.header.probeDirDeg.data;
    probeProj = cosd(probeOffsetDeg) .* probeNoise;

    % Build signed evidence matrices for change and no-change patches.
    %
    % Convention here:
    %   - For INC trials, evidence favoring the correct choice is positive as-is.
    %   - For DEC trials, multiply by -1 so positive still means evidence
    %     favoring the correct choice.
    %
    % If your previous kernels were not sign-aligned this way, edit helper only.
    [changeMat, noChangeMat, rfMat, oppMat] = ...
        signedNoiseMatrices(prefNoise, probeProj, changeSides, changeIndices);

    % Preallocate.
    kernels = nan(5, 2, nFrames);
    kVars   = nan(5, 2, nFrames);

    % Names just for bookkeeping.
    kernelNames = {'diff', 'change', 'noChange', 'RF', 'Opp'};

    % Compute per-step-type kernels.
    % dim2: 1 = DEC, 2 = INC
    for iStep = 1:2
        useStep = (changeIndices == iStep);

        % Outcome masks within this step type.
        hitMask  = useStep & (trialOutcomes == 0);
        errMask  = useStep & (trialOutcomes == 1);

        % Number of trials.
        nHit = nnz(hitMask);
        nErr = nnz(errMask);

        % Fill stats.
        kStats.nHit(iStep) = nHit;
        kStats.nErr(iStep) = nErr;
        kStats.stepName{iStep} = ternary(iStep == 1, 'DEC', 'INC');

        if nHit == 0 || nErr == 0
            warning('sessionKernelFromTrials:EmptyCondition', ...
                'No hit or error trials for %s in selected trial set.', kStats.stepName{iStep});
            continue;
        end

        % Pull matrices for hits/errors.
        chHit = changeMat(:, hitMask);
        chErr = changeMat(:, errMask);

        ncHit = noChangeMat(:, hitMask);
        ncErr = noChangeMat(:, errMask);

        rfHit = rfMat(:, hitMask);
        rfErr = rfMat(:, errMask);

        opHit = oppMat(:, hitMask);
        opErr = oppMat(:, errMask);

        % Mean(correct) - mean(wrong)
        kChange   = mean(chHit, 2) - mean(chErr, 2);
        kNoChange = mean(ncHit, 2) - mean(ncErr, 2);
        kRF       = mean(rfHit, 2) - mean(rfErr, 2);
        kOpp      = mean(opHit, 2) - mean(opErr, 2);
        kDiff     = kChange - kNoChange;

        kernels(1, iStep, :) = kDiff;
        kernels(2, iStep, :) = kChange;
        kernels(3, iStep, :) = kNoChange;
        kernels(4, iStep, :) = kRF;
        kernels(5, iStep, :) = kOpp;

        % Variance of difference of means by frame.
        vChange   = var(chHit, 0, 2) ./ nHit + var(chErr, 0, 2) ./ nErr;
        vNoChange = var(ncHit, 0, 2) ./ nHit + var(ncErr, 0, 2) ./ nErr;
        vRF       = var(rfHit, 0, 2) ./ nHit + var(rfErr, 0, 2) ./ nErr;
        vOpp      = var(opHit, 0, 2) ./ nHit + var(opErr, 0, 2) ./ nErr;

        % Approximate variance for diff = change - noChange.
        % Ignores covariance between change and noChange components.
        % If you want exact framewise variance, compute it from trialwise
        % diff matrices directly.
        vDiff = vChange + vNoChange;

        kVars(1, iStep, :) = vDiff;
        kVars(2, iStep, :) = vChange;
        kVars(3, iStep, :) = vNoChange;
        kVars(4, iStep, :) = vRF;
        kVars(5, iStep, :) = vOpp;
    end

    kStats.kernelNames = kernelNames;
    kStats.probeOffsetDeg = probeOffsetDeg;
    kStats.nFrames = nFrames;
    kStats.nTrialsUsed = numel(trialIdx);
end


function [changeMat, noChangeMat, rfMat, oppMat] = ...
    signedNoiseMatrices(prefNoise, probeProj, changeSides, changeIndices)
% Build trialwise evidence matrices in a common "correct-evidence-positive" sign.
%
% Inputs
%   prefNoise    : [nFrames x nTrials]
%   probeProj    : [nFrames x nTrials] projected probe noise along pref/drift axis
%   changeSides  : [1 x nTrials], 0=RF changed, 1=Opp changed
%   changeIndices: [1 x nTrials], 1=DEC, 2=INC
%
% Outputs
%   changeMat    : evidence on changed patch
%   noChangeMat  : evidence on unchanged patch
%   rfMat        : evidence on RF patch
%   oppMat       : evidence on Opp patch
%
% Assumption:
%   - RF patch corresponds to pref noise stream.
%   - Opp patch corresponds to projected probe stream.
%   - DEC trials are sign-flipped so positive always favors the correct response.

    nTrials = size(prefNoise, 2);

    % Trialwise sign alignment: DEC -> -1, INC -> +1
    stepSign = ones(1, nTrials);
    stepSign(changeIndices == 1) = -1;   % DEC
    stepSign(changeIndices == 2) = +1;   % INC

    % RF / Opp evidence matrices, aligned so + means evidence for correct answer.
    rfMat  = prefNoise  .* stepSign;
    oppMat = probeProj  .* stepSign;

    % Change / noChange determined by change side.
    changeMat   = nan(size(prefNoise));
    noChangeMat = nan(size(prefNoise));

    rfChanged = (changeSides == 0);
    oppChanged = (changeSides == 1);

    changeMat(:, rfChanged)    = rfMat(:, rfChanged);
    noChangeMat(:, rfChanged)  = oppMat(:, rfChanged);

    changeMat(:, oppChanged)   = oppMat(:, oppChanged);
    noChangeMat(:, oppChanged) = rfMat(:, oppChanged);
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end