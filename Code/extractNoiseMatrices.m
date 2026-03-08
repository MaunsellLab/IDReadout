function [prefMat, probeMat, trialOutcomes] = extractNoiseMatrices(header, trials, stepTypes, sideType)
% extractNoiseMatrices  Compile coherence-noise evidence into matrices for kernel estimation.
%
%   [prefMat, probeMat, trialOutcomes] = extractNoiseMatrices(header, trials, stepTypes, sideType)
%
%   INPUT:
%     header    : header struct from an IDKReadout data file
%     trials    : cell array of trial structs from an IDKReadout data file
%     stepTypes : vector of valid step types (1 = decrement, 2 = increment)
%     sideType  : which evidence stream to return
%                0 = DIFFERENCE evidence  (change - noChange)   [recommended for behavior]
%                1 = RF patch only
%                2 = opposite patch only
%                3 = change patch only     (legacy behavior)
%                4 = noChange patch only
%
%   OUTPUT:
%     prefMat       : m x nTrials matrix for preferred-direction noise evidence
%     probeMat      : m x nTrials matrix for probe-direction noise evidence
%     trialOutcomes : 1 x nTrials vector of trialEnd values (0=correct, 1=wrong)
%
%   Notes:
%   - This function aligns evidence by the per-trial timestamp vectors (e.g., changeTimesMS,
%     noChangeTimesMS). It does NOT assume synchronous noise frames across patches or streams.
%   - For sideType==0, the returned matrices are difference evidence:
%         dm_pref(t)  = pref_change(t)  - pref_noChange(t)
%         dm_probe(t) = probe_change(t) - probe_noChange(t)

nTrials  = numel(trials);
if nTrials == 0
    error('extractNoiseMatrices:EmptyInput', 'Input "trials" is empty.');
end

% ---- Frame/grid parameters ----
msPerVFrame = 1000.0 / header.frameRateHz.data(1);
m = round((header.preStepMS.data(1) + header.stepMS.data(1)) / msPerVFrame);

% ---- Find valid trials (certified, correct/wrong, requested step type, has noise) ----
validIdx = [];
trialOutcomes = [];
for k = 1:nTrials
    tr = trials{k};
    if ~isfield(tr, 'trialEnd') || ~isfield(tr, 'trialCertify') || ~isfield(tr, 'trial')
        continue;
    end
    tCert = tr.trialCertify.data;
    tEnd  = tr.trialEnd.data;
    tStep = tr.trial.data.changeIndex + 1; % 1=dec, 2=inc

    if ~(tCert == 0 && ismember(tEnd, [0 1]) && ismember(tStep, stepTypes))
        continue;
    end

    % Require at least one non-zero noise sample on the streams we need.
    hasChangeNoise   = ~(isNoNoise(tr.changePrefCohsPC.data(:))   && isNoNoise(tr.changeProbeCohsPC.data(:)));
    hasNoChangeNoise = ~(isNoNoise(tr.noChangePrefCohsPC.data(:)) && isNoNoise(tr.noChangeProbeCohsPC.data(:)));

    switch sideType
        case 0  % difference needs both
            hasNoise = hasChangeNoise || hasNoChangeNoise;
        case {1,2,3} % these rely on change patch if it is the selected patch on this trial, else noChange patch
            hasNoise = hasChangeNoise || hasNoChangeNoise;
        case 4
            hasNoise = hasNoChangeNoise;
        otherwise
            error('extractNoiseMatrices:BadSideType', 'Unknown sideType=%d.', sideType);
    end

    if ~hasNoise
        continue;
    end

    validIdx(end+1) = k;          %#ok<AGROW>
    trialOutcomes(end+1) = tEnd;  %#ok<AGROW>
end

nValid = numel(validIdx);
if nValid == 0
    error('extractNoiseMatrices:NoValidTrials', 'No valid matching trials were found.');
end

% ---- Preallocate output ----
prefMat  = nan(m, nValid);
probeMat = nan(m, nValid);

% ---- Extract evidence ----
for kk = 1:nValid
    tr = trials{validIdx(kk)};

    % Build matrices for each patch, aligned by that patch's own timestamps.
    prefChange    = fillFromTimes(tr.changePrefCohsPC.data(:),    tr.changeTimesMS.data(:),    m, msPerVFrame);
    probeChange   = fillFromTimes(tr.changeProbeCohsPC.data(:),   tr.changeTimesMS.data(:),    m, msPerVFrame);
    prefNoChange  = fillFromTimes(tr.noChangePrefCohsPC.data(:),  tr.noChangeTimesMS.data(:),  m, msPerVFrame);
    probeNoChange = fillFromTimes(tr.noChangeProbeCohsPC.data(:), tr.noChangeTimesMS.data(:),  m, msPerVFrame);

    switch sideType
        case 0  % difference evidence
            prefMat(:, kk)  = prefChange    - prefNoChange;
            probeMat(:, kk) = probeChange   - probeNoChange;

        case 3  % change patch only (legacy behavior)
            prefMat(:, kk)  = prefChange;
            probeMat(:, kk) = probeChange;

        case 4  % noChange patch only
            prefMat(:, kk)  = prefNoChange;
            probeMat(:, kk) = probeNoChange;

        case 1  % RF patch only
            if tr.trial.data.changeSide == 0
                % RF is the change patch on this trial
                prefMat(:, kk)  = prefChange;
                probeMat(:, kk) = probeChange;
            else
                % RF is the noChange patch on this trial
                prefMat(:, kk)  = prefNoChange;
                probeMat(:, kk) = probeNoChange;
            end

        case 2  % opposite patch only
            if tr.trial.data.changeSide == 1
                % opposite patch is the change patch on this trial
                prefMat(:, kk)  = prefChange;
                probeMat(:, kk) = probeChange;
            else
                % opposite patch is the noChange patch on this trial
                prefMat(:, kk)  = prefNoChange;
                probeMat(:, kk) = probeNoChange;
            end
    end
end
end

% ===== Helper: reconstruct sample-and-hold time series on the video-frame grid =====
function v = fillFromTimes(cohsPC, timesMS, m, msPerVFrame)
    v = nan(m,1);
    if isempty(timesMS) || isempty(cohsPC)
        return;
    end
    nTimes = min(numel(timesMS), numel(cohsPC));
    timesMS = timesMS(1:nTimes);
    cohsPC  = cohsPC(1:nTimes);

    for tIndex = 1:nTimes
        % Map timestamp to 1-indexed video-frame bins (sample-and-hold).
        % Use floor(+1) rather than round to avoid occasional overlaps/skips.
        t0 = timesMS(tIndex);
        theVFrame = floor(t0 / msPerVFrame) + 1;
        if theVFrame < 1
            theVFrame = 1;
        elseif theVFrame > m
            continue;
        end

        if tIndex < nTimes
            t1 = timesMS(tIndex + 1);
            nextVFrame = floor(t1 / msPerVFrame) + 1;
        else
            nextVFrame = m + 1;
        end

        if nextVFrame <= theVFrame
            continue;
        end
        if nextVFrame > m + 1
            nextVFrame = m + 1;
        end
        v(theVFrame:nextVFrame-1) = cohsPC(tIndex);
    end
end

function tf = isNoNoise(x)
    tf = isscalar(x) && isequal(x, 0);
end
