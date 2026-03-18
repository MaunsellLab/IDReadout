function R = strategyForensics(sessionOrFolder, stepTypes, opts)
% strategyForensics  Single-session OR batch forensics for drift+probe noise strategy.
%
% USAGE (single session):
%   R = strategyForensics('IDReadout_Meetz_20260226.mat', [1 2], opts)
%
% USAGE (batch folder):
%   R = strategyForensics('../Data/', [1 2], opts)
%   R = strategyForensics([], [1 2], opts)   % uses defaultDataFolder()
%
% Batch mode returns:
%   R.summary : per-session table of key metrics and GLM coefficients
%   R.sessions: struct array with per-session outputs (optional, can be large)
%
% Requires fields per trial:
%   tr.targetChosen.data in {0,1}
%   tr.trial.data.changeSide in {0,1}
%   base/step coherence somewhere retrievable by getBaseStep()
%   tr.changePrefCohsPC.data, tr.noChangePrefCohsPC.data (deltas)
%   tr.changeProbeCohsPC.data, tr.noChangeProbeCohsPC.data (coupled probe delta)
%   tr.changeTimesMS.data, tr.noChangeTimesMS.data
%
% Splits analyses by stepType: 2=INC, 1=DEC.

if nargin < 1, sessionOrFolder = defaultDataFolder(); end
if nargin < 2 || isempty(stepTypes), stepTypes = [1 2]; end
if nargin < 3, opts = struct(); end

% Defaults
if ~isfield(opts,'decisionMS'),        opts.decisionMS = 250; end
if ~isfield(opts,'cosProbe'),          opts.cosProbe = cosd(45); end
if ~isfield(opts,'printModels'),       opts.printModels = true; end
if ~isfield(opts,'thrFrac'),           opts.thrFrac = 0.70; end
if ~isfield(opts,'hiCoh'),             opts.hiCoh = 65; end
if ~isfield(opts,'loCoh'),             opts.loCoh = 15; end
if ~isfield(opts,'useRT'),             opts.useRT = true; end
if ~isfield(opts,'excludedFiles'),     opts.excludedFiles = defaultExcluded(); end
if ~isfield(opts,'storePerSession'),   opts.storePerSession = false; end
if ~isfield(opts,'batchPlot'),         opts.batchPlot = false; end

% Heuristic: if input is empty or a folder -> batch mode
% if isempty(sessionOrFolder)
%     sessionOrFolder = defaultDataFolder();
% end

if isfolder(sessionOrFolder)
    % -------- BATCH MODE --------
    dataFolder = char(sessionOrFolder);
    matFiles = dir(fullfile(dataFolder, '*.mat'));

    if isempty(matFiles)
        error('strategyForensics:NoMatFiles', 'No .mat files found in %s', dataFolder);
    end

    % In batch, default to no plotting unless explicitly requested
    if ~isfield(opts,'plot'), opts.plot = false; end

    summaryRows = [];
    sessions = struct([]);

    fprintf('\nBatch forensics in folder: %s\n', dataFolder);

    for f = 1:numel(matFiles)
        fileName = matFiles(f).name;
        [~, baseName, ~] = fileparts(fileName);
        if endsWith(fileName, '_fileInfo.mat')
            continue;
        end
        if ~isempty(opts.excludedFiles) && ismember(baseName, opts.excludedFiles)
            continue;
        end
        fullPath = fullfile(dataFolder, fileName);
        fprintf('\nProcessing %s (%d/%d)\n', fileName, f, numel(matFiles));
        try
            Ri = strategyForensics_one(fullPath, stepTypes, opts);
            row = Ri.summaryRow;
            row.file = string(fileName);
            summaryRows = [summaryRows; row]; %#ok<AGROW>

            if opts.storePerSession
                sessions(end+1).file = fileName; %#ok<AGROW>
                sessions(end).R = Ri; 
            end
        catch ME
            warning('strategyForensics:SessionFailed', 'FAILED %s: %s', fileName, ME.message);
            row = emptySummaryRow();
            row.file = string(fileName);
            row.failed = true;
            row.failMsg = string(ME.message);
            summaryRows = [summaryRows; row]; %#ok<AGROW>
        end
    end

    % Convert struct rows to table
    summary = struct2table(summaryRows);

    % Optional batch plots
    if opts.batchPlot
        batchPlots(summary);
    end

    R = struct();
    R.summary = summary;
    if opts.storePerSession
        R.sessions = sessions;
    end
    R.meta = struct('folder', dataFolder, 'opts', opts);
    return;
else
    % -------- SINGLE SESSION MODE --------
    if ~isfield(opts, 'plot'), opts.plot = true; end

    sessionOrFolder = char(sessionOrFolder);

    if isfile(sessionOrFolder)
        % already a valid full or relative path
        sessionFile = sessionOrFolder;
    else
        % treat as bare filename and look in default data folder
        candidate = fullfile(defaultDataFolder(), sessionOrFolder);
        if isfile(candidate)
            sessionFile = candidate;
        else
            error('strategyForensics:SessionFileNotFound', ...
                'Could not find session file: %s', sessionOrFolder);
        end
    end

    R = strategyForensics_one(sessionFile, stepTypes, opts);
end
end

% ===================== SINGLE SESSION CORE =====================

function R = strategyForensics_one(sessionFile, stepTypes, opts)

S = load(sessionFile);
header = S.header;
trials = S.trials;

msPerVFrame = 1000.0 / header.frameRateHz.data;
preStepMS   = header.preStepMS.data(1);
totalVFrames = round((header.preStepMS.data(1) + header.stepMS.data(1)) / msPerVFrame);

postStartFrame = max(1, round(preStepMS / msPerVFrame));
postEndFrame   = min(totalVFrames, postStartFrame + round(opts.decisionMS / msPerVFrame) - 1);
postFrames     = postStartFrame:postEndFrame;
nStepVFrames   = numel(postFrames);
proj = 2 * opts.cosProbe; % coupled probes at ±45°

% ---- Select noise trials ----
validIdx = [];
EOT      = [];    % 0 correct, 1 wrong
stepType = [];    % 1 dec, 2 inc
baseCoh  = [];    % per-trial base coherence (PC)
stepCoh  = [];    % per-trial step coherence (PC)
RT       = [];

for k = 1:numel(trials)
    tr = trials{k};
    tCert = tr.trialCertify.data;
    tEnd  = tr.trialEnd.data;
    tStep = tr.trial.data.changeIndex + 1; % 1 dec, 2 inc

    hasNoise = ~(isNoNoise(tr.changePrefCohsPC.data(:)) && isNoNoise(tr.changeProbeCohsPC.data(:)));
    if ~(tCert == 0 && ismember(tEnd,[0 1]) && ismember(tStep, stepTypes) && hasNoise)
        continue;
    end
    b = getScalarField(tr.trial.data, 'baseCohPC');
    s = getScalarField(tr.trial.data, 'stepCohPC');
    if isnan(b) || isnan(s)
        continue;
    end
    validIdx(end+1) = k; %#ok<AGROW>
    EOT(end+1)      = tEnd; %#ok<AGROW>
    stepType(end+1) = tStep; %#ok<AGROW>
    baseCoh(end+1)  = b; %#ok<AGROW>
    stepCoh(end+1)  = s; %#ok<AGROW>

    if opts.useRT && isfield(tr,'reactTimeMS') && isfield(tr.reactTimeMS,'data') && ~isempty(tr.reactTimeMS.data)
        RT(end+1) = tr.reactTimeMS.data(1); %#ok<AGROW>
    else
        RT(end+1) = NaN; %#ok<AGROW>
    end
end

nValid = numel(validIdx);
if nValid == 0
    error('strategyForensics:NoValidTrials','No valid noise trials found in %s', sessionFile);
end

% isCorrect is [nValid x 1] logical (or 0/1)
isCorrect = (EOT(:) == 0);

% ---- Extract noise deltas per patch per frame ----
prefCh  = nan(nStepVFrames, nValid);
probeCh = nan(nStepVFrames, nValid);
prefNc  = nan(nStepVFrames, nValid);
probeNc = nan(nStepVFrames, nValid);
chgSide  = [];    % 0=RF has change, 1=opposite has change
typeStep  = [];

for i = 1:nValid
    tr = trials{validIdx(i)};

    prefC  = tr.changePrefCohsPC.data(:);
    probeC = tr.changeProbeCohsPC.data(:);
    tC     = tr.changeTimesMS.data(:);

    prefN  = tr.noChangePrefCohsPC.data(:);
    probeN = tr.noChangeProbeCohsPC.data(:);
    tN     = tr.noChangeTimesMS.data(:);

    xC = fillFrames(prefC,  tC, totalVFrames, msPerVFrame);
    yC = fillFrames(probeC, tC, totalVFrames, msPerVFrame);
    xN = fillFrames(prefN,  tN, totalVFrames, msPerVFrame);
    yN = fillFrames(probeN, tN, totalVFrames, msPerVFrame);

    prefCh(:, i)  = xC(postFrames);
    probeCh(:, i) = yC(postFrames);
    prefNc(:, i)  = xN(postFrames);
    probeNc(:, i) = yN(postFrames);
    chgSide(end+1) = tr.trial.data.changeSide; %#ok<AGROW> 0 = RF is change patch, 1 = RF is no-change patch
    typeStep(end+1) = tr.trial.data.changeIndex; %#ok<AGROW>  1 dec, 2 inc 
end

% ---- Absolute drift-axis noise coherence per patch ----
nCh = prefCh + proj * probeCh;
nNc = prefNc + proj * probeNc;
dn  = nCh - nNc;

isRFChange = (chgSide(:) == 0);     % RF contains the change patch on these trials
isDec = typeStep(:) == 1;
% isInc = ~isDec;
nRF  = nan(nStepVFrames, nValid);
nOpp = nan(nStepVFrames, nValid);

nRF(:, isRFChange)   = nCh(:, isRFChange);
nOpp(:, isRFChange)  = nNc(:, isRFChange);

nRF(:, ~isRFChange)  = nNc(:, ~isRFChange);
nOpp(:, ~isRFChange) = nCh(:, ~isRFChange);

% Per-trial window averages (in % coherence units)
% +1 INC, -1 DEC
stepSign = ones(nValid,1);
stepSign(isDec) = -1;

nCh_signed = nCh .* stepSign';
nNc_signed = nNc .* stepSign';
dn_signed  = dn  .* stepSign';

nRF_signed  = nRF  .* stepSign';
nOpp_signed = nOpp .* stepSign';

nRF_bar  = mean(nRF_signed,  1)';
nOpp_bar = mean(nOpp_signed, 1)';
dn_signed_bar   = mean(dn_signed,   1)';
sumNcNeg = -sum(nNc_signed, 1)';

% ============================================================
% DIAGNOSTIC: compare candidate scalar evidence variables
% ============================================================

% sumDM_chNc  = sum(dn_signed, 1)';                  % current: change - nochange
% sumDM_rfOpp = sum(nRF_signed - nOpp_signed, 1)';  % RF - Opp
% sumCh       = sum(nCh_signed, 1)';                % change patch only
% sumNcNeg    = -sum(nNc_signed, 1)';               % no-change patch, sign flipped
% 
% fprintf('\n[Decision-variable diagnostics]\n');
% localDiagScalar('sumDM_chNc',  sumDM_chNc,  isCorrect);
% localDiagScalar('sumDM_rfOpp', sumDM_rfOpp, isCorrect);
% localDiagScalar('sumCh',       sumCh,       isCorrect);
% localDiagScalar('sumNcNeg',    sumNcNeg,    isCorrect);
% 
% % optional quick plots
% localPvPlot(201, sumDM_chNc,  isCorrect, 'sumDM ch-nc');
% localPvPlot(202, sumDM_rfOpp, isCorrect, 'sumDM RF-Opp');
% localPvPlot(203, sumCh,       isCorrect, 'sum change');
% localPvPlot(204, sumNcNeg,    isCorrect, 'neg sum no-change');
% 
% sumDM_fixed = sum((nRF - nOpp).*stepSign',1)';
% 
% localDiagScalar('sumDM_fixed', sumDM_fixed, isCorrect);
% localPvPlot(205, sumDM_fixed, isCorrect, 'sumDM fixed');
% 
% 
% 
% PvIN_Plot(100, dn(:, isInc),  isCorrect(isInc), [], 'INC Steps');
% PvIN_Plot(101, -dn(:, isDec), isCorrect(isDec), [], 'DEC Steps');
% PvIN_Plot(102, dn_signed,     isCorrect, [], 'INC & DEC Steps');
% 
% fprintf('corr(sumDM_chNc,  sumDM_rfOpp) = %.4f\n', ...
%     corr(sumDM_chNc, sumDM_rfOpp, 'rows', 'complete'));

% --- Conditioned shifts (Correct - Error) ---
% nRF_decAligned  = nRF_bar  .* stepSign;
% nOpp_decAligned = nOpp_bar .* stepSign;
dn_decAligned   = dn_signed_bar   .* stepSign;

% d_nRF  = mean(nRF_decAligned(isCorrect))  - mean(nRF_decAligned(~isCorrect));
% d_nOpp = mean(nOpp_decAligned(isCorrect)) - mean(nOpp_decAligned(~isCorrect));
d_dn   = mean(dn_decAligned(isCorrect))   - mean(dn_decAligned(~isCorrect));
% Optional SEMs (difference of means)
semDiff = @(x,y) sqrt(var(x)/numel(x) + var(y)/numel(y));
% sem_nRF  = semDiff(nRF_bar(isCorrect),  nRF_bar(~isCorrect));
% sem_nOpp = semDiff(nOpp_bar(isCorrect), nOpp_bar(~isCorrect));
sem_dn_signed   = semDiff(dn_signed_bar(isCorrect),   dn_signed_bar(~isCorrect));

% fprintf('Δmean(RF)   = %.3f (SEM %.3f) %%coh\n', d_nRF, sem_nRF);
% fprintf('Δmean(Opp)  = %.3f (SEM %.3f) %%coh\n', d_nOpp, sem_nOpp);
fprintf('Δmean(change - noChange)   = %.3f%% coh (SEM %.3f)\n', d_dn, sem_dn_signed);
% fprintf('Check: Δdm ≈ ΔRF - ΔOpp = %.3f %%coh\n', d_mRF - d_mOpp);

obsMax = max(abs(dn_signed(:)));
thr = opts.thrFrac * obsMax;

% ============================================================
% PURE-NOISE PREDICTORS
% dn_signed must be signed so that + means evidence for correct
% ============================================================

[nFrames, nTrials] = size(dn_signed);

% ---- basic net evidence ----
sumDM  = sum(dn_signed, 1)';          % integrated signed evidence
meanDM = mean(dn_signed, 1)';         % average signed evidence per frame

% ---- occupancy ----
isPos  = dn_signed > 0;
isNeg  = dn_signed < 0;
isZero = dn_signed == 0;

countPos  = sum(isPos,  1)';
countNeg  = sum(isNeg,  1)';
countZero = sum(isZero, 1)';

% ---- run structure ----
longestPosRun = longestRunPerTrial(isPos);
longestNegRun = longestRunPerTrial(isNeg);

% number of sign switches, ignoring zeros by working from sign(dn_signed)
sgn = sign(dn_signed);   % -1, 0, +1
nSwitch = zeros(nTrials,1);
for t = 1:nTrials
    s = sgn(:,t);
    s = s(s ~= 0);       % ignore zero-evidence frames
    if numel(s) <= 1
        nSwitch(t) = 0;
    else
        nSwitch(t) = sum(diff(s) ~= 0);
    end
end

% ---- temporal weighting ----
midFrame = floor(nFrames/2);

sumEarly = sum(dn_signed(1:midFrame, :), 1)';
sumLate  = sum(dn_signed(midFrame+1:end, :), 1)';
dEarlyLate = sumEarly - sumLate;

% optional linear ramp projection:
% positive means evidence was weighted earlier more than later
ramp = linspace(1, -1, nFrames)';     % early positive, late negative
rampDM = (ramp' * dn_signed)';        % nTrials x 1

% ---- cumulative evidence features ----
cumDM = cumsum(dn_signed, 1);
maxCumDM = max(cumDM, [], 1)';        % strongest favorable excursion
minCumDM = min(cumDM, [], 1)';        % strongest unfavorable excursion
% finalCumDM = cumDM(end, :)';          % same as sumDM, included for clarity

% ---- outcome / trial info ----
y = isCorrect(:);

% If you still want bookkeeping variables in the table, keep them;
% but they are not needed for the pure-noise GLMs.
tbl = table(y, stepType(:), RT(:), ...
    sumDM, meanDM, ...
    countPos, countNeg, countZero, ...
    longestPosRun, longestNegRun, nSwitch, ...
    sumEarly, sumLate, dEarlyLate, rampDM, ...
    maxCumDM, minCumDM, ...
    'VariableNames', { ...
    'isCorrect','stepType','RT', ...
    'sumDM','meanDM', ...
    'countPos','countNeg','countZero', ...
    'longestPosRun','longestNegRun','nSwitch', ...
    'sumEarly','sumLate','dEarlyLate','rampDM', ...
    'maxCumDM','minCumDM'});

% Per-trial window means (% coh)
nRF_bar  = mean(nRF_signed,  1)';
nOpp_bar = mean(nOpp_signed, 1)';
dn_bar   = mean(dn_signed,   1)';          % change - nochange (global sign)
dnRF_bar = mean(nRF_signed - nOpp_signed, 1)';     % RF - Opp (flips with changeSide)

isInc = (tbl.stepType==2);
isDec = (tbl.stepType==1);
y     = tbl.isCorrect;

% Helper printer
prt = @(tag, dRF, dOpp, dDM, dDMRF) fprintf( ...
  '%s: ΔRF=%.3f  ΔOpp=%.3f  Δdm=%.3f  Δ(dmRF)=%.3f  (all in %%coh)\n', ...
  tag, dRF, dOpp, dDM, dDMRF);

% INC (correct - error)
dRF  = mean(nRF_bar(isInc & y))  - mean(nRF_bar(isInc & ~y));
dOpp = mean(nOpp_bar(isInc & y)) - mean(nOpp_bar(isInc & ~y));
dDM  = mean(dn_bar(isInc & y))   - mean(dn_bar(isInc & ~y));
dDMRF= mean(dnRF_bar(isInc & y)) - mean(dnRF_bar(isInc & ~y));
prt('INC', dRF, dOpp, dDM, dDMRF);

% DEC: flip sign of dn_signed so "correct-direction evidence" is positive
% DEC
dRF  = mean(nRF_bar(isDec & y))  - mean(nRF_bar(isDec & ~y));
dOpp = mean(nOpp_bar(isDec & y)) - mean(nOpp_bar(isDec & ~y));
dm_bar_dec = -dn_bar;   % sign flip for DEC trials
dDM  = mean(dm_bar_dec(isDec & y)) - mean(dm_bar_dec(isDec & ~y));
dDMRF = mean(dnRF_bar(isDec & y)) - mean(dnRF_bar(isDec & ~y));
% dRF  = mean(mRF_bar(isDec & y))  - mean(mRF_bar(isDec & ~y));
% dOpp = mean(mOpp_bar(isDec & y)) - mean(mOpp_bar(isDec & ~y));
% dDM  = mean((-dm_bar)(isDec & y)) - mean((-dm_bar)(isDec & ~y));   % <-- sign flip
% dDMRF= mean(dmRF_bar(isDec & y)) - mean(dmRF_bar(isDec & ~y));
prt('DEC', dRF, dOpp, dDM, dDMRF);

% And the identity check (per condition) should now be:
% For INC with change-nochange sign:      Δdm ≈ ΔmCh - ΔmNc  (not ΔRF-ΔOpp)
% For RF-minus-Opp sign (dmRF):           Δ(dmRF) ≈ ΔRF - ΔOpp
fprintf('Check (both conds): Δ(dmRF) ≈ ΔRF - ΔOpp\n');

% if opts.printModels
%     fprintf('\n=== %s ===\n', sessionFile);
%     fprintf('Selected noise trials: %d (INC %d, DEC %d)\n', height(tbl), height(tblInc), height(tblDec));
%     fprintf('P(correct): all %.3f | INC %.3f | DEC %.3f\n', mean(tbl.isCorrect), mean(tblInc.isCorrect), mean(tblDec.isCorrect));
%     fprintf('Post-step window: frames %d..%d (%d frames, %.1f ms)\n', postStartFrame, postEndFrame, nStepVFrames, nStepVFrames*msPerVFrame);
%     fprintf('proj (2cos45): %.3f | thr: %.3f (thrFrac %.2f of observed max |dn_signed| %.3f)\n', proj, thr, opts.thrFrac, obsMax);
% end

if opts.printModels
    fprintf('\n=== %s ===\n', sessionFile);
    fprintf('Selected noise trials: %d\n', height(tbl));
    fprintf('P(correct): %.3f\n', mean(tbl.isCorrect));
    fprintf('Post-step window: frames %d..%d (%d frames, %.1f ms)\n', ...
        postStartFrame, postEndFrame, nFrames, nFrames*msPerVFrame);

    dMean = mean(tbl.sumDM(tbl.isCorrect==1)) - mean(tbl.sumDM(tbl.isCorrect==0));
    fprintf('Δmean(sumDM): %.3f %%coh·frames\n', dMean);

    dMeanFrame = mean(tbl.meanDM(tbl.isCorrect==1)) - mean(tbl.meanDM(tbl.isCorrect==0));
    fprintf('Δmean(meanDM): %.3f %%coh/frame\n', dMeanFrame);

    fprintf('Mean longestPosRun: correct %.2f | error %.2f frames\n', ...
        mean(tbl.longestPosRun(tbl.isCorrect==1)), ...
        mean(tbl.longestPosRun(tbl.isCorrect==0)));

    fprintf('Mean longestNegRun: correct %.2f | error %.2f frames\n', ...
        mean(tbl.longestNegRun(tbl.isCorrect==1)), ...
        mean(tbl.longestNegRun(tbl.isCorrect==0)));

    fprintf('Mean nSwitch:       correct %.2f | error %.2f\n', ...
        mean(tbl.nSwitch(tbl.isCorrect==1)), ...
        mean(tbl.nSwitch(tbl.isCorrect==0)));

    fprintf('Mean dEarlyLate:    correct %.3f | error %.3f\n', ...
        mean(tbl.dEarlyLate(tbl.isCorrect==1)), ...
        mean(tbl.dEarlyLate(tbl.isCorrect==0)));
end

models = struct();
if height(tbl) >= 40
    % ---- baseline ----
    m0 = safeFitglm(tbl, ...
        'isCorrect ~ 1', ...
        'Distribution', 'binomial');

    % ---- net evidence only ----
    m1 = safeFitglm(tbl, ...
        'isCorrect ~ sumDM', ...
        'Distribution', 'binomial');

    % ---- net evidence + run structure ----
    m2 = safeFitglm(tbl, ...
        'isCorrect ~ sumDM + longestPosRun + longestNegRun + nSwitch', ...
        'Distribution', 'binomial', ...
        'LikelihoodPenalty', 'jeffreys-prior');

    % ---- net evidence + temporal weighting ----
    % Using sumEarly + sumLate lets the model discover weighting asymmetry.
    % If collinearity is annoying, replace with: sumDM + dEarlyLate
    m3 = safeFitglm(tbl, ...
        'isCorrect ~ sumDM + dEarlyLate', ...
        'Distribution', 'binomial', ...
        'LikelihoodPenalty', 'jeffreys-prior');

    % ---- optional: cumulative evidence structure ----
    m4 = safeFitglm(tbl, ...
        'isCorrect ~ sumDM + maxCumDM + minCumDM', ...
        'Distribution', 'binomial', ...
        'LikelihoodPenalty', 'jeffreys-prior');

    models.pooled = struct('m0',m0,'m1',m1,'m2',m2,'m3',m3,'m4',m4);

    if opts.printModels
        fprintf('\n[Pooled signed-evidence models]\n');
        fprintf('Deviance m0 %.2f\n', m0.Deviance);
        fprintf('Deviance m1 %.2f   [sumDM]\n', m1.Deviance);
        fprintf('Deviance m2 %.2f   [sumDM + runs]\n', m2.Deviance);
        fprintf('Deviance m3 %.2f   [sumDM + early/late]\n', m3.Deviance);
        fprintf('Deviance m4 %.2f   [sumDM + cumulative excursions]\n', m4.Deviance);

        fprintf('Dev explained m1: %.5f\n', (m0.Deviance - m1.Deviance)/m0.Deviance);
        fprintf('Dev explained m2: %.5f\n', (m0.Deviance - m2.Deviance)/m0.Deviance);
        fprintf('Dev explained m3: %.5f\n', (m0.Deviance - m3.Deviance)/m0.Deviance);
        fprintf('Dev explained m4: %.5f\n', (m0.Deviance - m4.Deviance)/m0.Deviance);

        fprintf('Increment over m1 (runs):       %.5f\n', (m1.Deviance - m2.Deviance)/m0.Deviance);
        fprintf('Increment over m1 (early/late): %.5f\n', (m1.Deviance - m3.Deviance)/m0.Deviance);
        fprintf('Increment over m1 (cum):        %.5f\n', (m1.Deviance - m4.Deviance)/m0.Deviance);
    end
end

prefVals  = cellfun(@(t) t.trial.data.prefCohNoisePC,  trials);
probeVals = cellfun(@(t) t.trial.data.probeCohNoisePC, trials);  
prefCohNoisePC  = max(prefVals(prefVals  > 0));
probeCohNoisePC = max(probeVals(probeVals > 0));

trialForecast = struct();
if isfield(models,'pooled') && isfield(models.pooled,'m1')
    m1 = models.pooled.m1;
    beta1 = coefPvOnly(m1, 'sumDM');

    trialForecast = predictKernelTrials(beta1, mean(tbl.isCorrect), ...
        std(tbl.sumDM, 0, 'omitnan'), ...
        prefCohNoisePC, probeCohNoisePC, proj, msPerVFrame, ...
        'refWindowMS', numel(postFrames) * msPerVFrame, ...
        'windowMS', [50 125 250], ...
        'zTarget', 2);
end

summaryRow = makeSummaryRow(sessionFile, tbl, models, opts, trialForecast);

R = struct();
R.tblAll = tbl;
R.models = models;
R.meta = struct( ...
    'postFrames', postFrames, ...
    'msPerVFrame', msPerVFrame, ...
    'proj', proj, ...
    'thr', thr, ...
    'obsMaxAbsDm', obsMax);
R.summaryRow = summaryRow;
R.trialForecast = trialForecast;

% ----- pooled slope interpretation -----
m1 = R.models.pooled.m1;
beta1 = m1.Coefficients.Estimate('sumDM');
se1   = m1.Coefficients.SE('sumDM');

T = numel(R.meta.postFrames);
dSum_50to75   = 1.0986 / beta1;
dMeanDm_50to75 = dSum_50to75 / T;

% ----- pooled difference-of-means summary -----
x = R.tblAll.sumDM;
y = R.tblAll.isCorrect;

dMu = mean(x(y==1)) - mean(x(y==0));
semMu = sqrt(var(x(y==1))/sum(y==1) + var(x(y==0))/sum(y==0));

% Optional:
fprintf('Pooled: beta1 = %.4g (SE %.4g)\n', beta1, se1);
fprintf('Pooled: ΔsumDM for 50%%→75%% = %.3f %%coh·frames\n', dSum_50to75);
fprintf('Pooled: Equivalent Δmean dn_signed = %.3f %%coh\n', dMeanDm_50to75);
fprintf('Pooled: Δmean(sumDM) = %.2f %%coh·frames => %.3f %%coh\n', dMu, dMu/T);
fprintf('Pooled: SEM(Δmean dn_signed) = %.3f %%coh\n', semMu/T);
end

% ===================== SUMMARY ROW =====================

function row = makeSummaryRow(sessionFile, tbl, models, ~, trialForecast)

if nargin < 5 || isempty(trialForecast)
    trialForecast = struct();
end
row = emptySummaryRow();
[~, baseName, ext] = fileparts(sessionFile);
row.file = string(sessionFile);
row.session = string([baseName ext]);
row.failed = false;

% ---------------- basic trial counts ----------------
row.nTrials = height(tbl);
if row.nTrials > 0
    row.pCorrect = mean(tbl.isCorrect);
end

if ismember('stepType', tbl.Properties.VariableNames)
    isInc = tbl.stepType == 2;
    isDec = tbl.stepType == 1;
    row.nInc = sum(isInc);
    row.nDec = sum(isDec);
    if row.nInc > 0, row.pInc = mean(tbl.isCorrect(isInc)); end
    if row.nDec > 0, row.pDec = mean(tbl.isCorrect(isDec)); end
end

% ---------------- descriptive predictor means ----------------
if row.nTrials > 0
    yc = tbl.isCorrect == 1;
    ye = tbl.isCorrect == 0;

    if any(yc) && any(ye)
        if ismember('sumDM', tbl.Properties.VariableNames)
            row.dMean_sumDM = mean(tbl.sumDM(yc)) - mean(tbl.sumDM(ye));
        end
        if ismember('meanDM', tbl.Properties.VariableNames)
            row.dMean_meanDM = mean(tbl.meanDM(yc)) - mean(tbl.meanDM(ye));
        end
        if ismember('longestPosRun', tbl.Properties.VariableNames)
            row.dMean_longestPosRun = mean(tbl.longestPosRun(yc)) - mean(tbl.longestPosRun(ye));
        end
        if ismember('longestNegRun', tbl.Properties.VariableNames)
            row.dMean_longestNegRun = mean(tbl.longestNegRun(yc)) - mean(tbl.longestNegRun(ye));
        end
        if ismember('nSwitch', tbl.Properties.VariableNames)
            row.dMean_nSwitch = mean(tbl.nSwitch(yc)) - mean(tbl.nSwitch(ye));
        end
        if ismember('dEarlyLate', tbl.Properties.VariableNames)
            row.dMean_dEarlyLate = mean(tbl.dEarlyLate(yc)) - mean(tbl.dEarlyLate(ye));
        end
        if ismember('maxCumDM', tbl.Properties.VariableNames)
            row.dMean_maxCumDM = mean(tbl.maxCumDM(yc)) - mean(tbl.maxCumDM(ye));
        end
        if ismember('minCumDM', tbl.Properties.VariableNames)
            row.dMean_minCumDM = mean(tbl.minCumDM(yc)) - mean(tbl.minCumDM(ye));
        end
    end
end

% ---------------- pooled model summaries ----------------
if isfield(models,'pooled') && isfield(models.pooled,'m0')
    mp = models.pooled;

    if isfield(mp,'m0'), row.dev_m0 = mp.m0.Deviance; end
    if isfield(mp,'m1'), row.dev_m1 = mp.m1.Deviance; end
    if isfield(mp,'m2'), row.dev_m2 = mp.m2.Deviance; end
    if isfield(mp,'m3'), row.dev_m3 = mp.m3.Deviance; end
    if isfield(mp,'m4'), row.dev_m4 = mp.m4.Deviance; end

    if ~isnan(row.dev_m0) && ~isnan(row.dev_m1)
        row.dDev_m1 = row.dev_m0 - row.dev_m1;
        row.devExpl_m1 = row.dDev_m1 / row.dev_m0;
    end
    if ~isnan(row.dev_m0) && ~isnan(row.dev_m2)
        row.dDev_m2 = row.dev_m0 - row.dev_m2;
        row.devExpl_m2 = row.dDev_m2 / row.dev_m0;
    end
    if ~isnan(row.dev_m0) && ~isnan(row.dev_m3)
        row.dDev_m3 = row.dev_m0 - row.dev_m3;
        row.devExpl_m3 = row.dDev_m3 / row.dev_m0;
    end
    if ~isnan(row.dev_m0) && ~isnan(row.dev_m4)
        row.dDev_m4 = row.dev_m0 - row.dev_m4;
        row.devExpl_m4 = row.dDev_m4 / row.dev_m0;
    end

    if ~isnan(row.dev_m1) && ~isnan(row.dev_m2)
        row.incOver_m1_m2 = row.dev_m1 - row.dev_m2;
    end
    if ~isnan(row.dev_m1) && ~isnan(row.dev_m3)
        row.incOver_m1_m3 = row.dev_m1 - row.dev_m3;
    end
    if ~isnan(row.dev_m1) && ~isnan(row.dev_m4)
        row.incOver_m1_m4 = row.dev_m1 - row.dev_m4;
    end

    % ---- coefficients from m1
    if isfield(mp,'m1')
        [row.b_m1_sumDM, row.p_m1_sumDM] = coefPv(mp.m1, 'sumDM');
    end

    % ---- coefficients from m2 (run structure)
    if isfield(mp,'m2')
        [row.b_m2_sumDM, row.p_m2_sumDM] = coefPv(mp.m2, 'sumDM');
        [row.b_m2_longestPosRun, row.p_m2_longestPosRun] = coefPv(mp.m2, 'longestPosRun');
        [row.b_m2_longestNegRun, row.p_m2_longestNegRun] = coefPv(mp.m2, 'longestNegRun');
        [row.b_m2_nSwitch, row.p_m2_nSwitch] = coefPv(mp.m2, 'nSwitch');
    end

    % ---- coefficients from m3 (temporal weighting)
    if isfield(mp,'m3')
        [row.b_m3_sumDM, row.p_m3_sumDM] = coefPv(mp.m3, 'sumDM');
        [row.b_m3_dEarlyLate, row.p_m3_dEarlyLate] = coefPv(mp.m3, 'dEarlyLate');
    end

    % ---- coefficients from m4 (cumulative excursions)
    if isfield(mp,'m4')
        [row.b_m4_sumDM, row.p_m4_sumDM] = coefPv(mp.m4, 'sumDM');
        [row.b_m4_maxCumDM, row.p_m4_maxCumDM] = coefPv(mp.m4, 'maxCumDM');
        [row.b_m4_minCumDM, row.p_m4_minCumDM] = coefPv(mp.m4, 'minCumDM');
    end

    % ---------------- trial-count forecasts ----------------
% ---------------- trial-count forecasts ----------------
if isstruct(trialForecast) && ~isempty(trialForecast)

    % main forecast values
    if isfield(trialForecast, 'windowMS')
        w = trialForecast.windowMS;

        if isfield(trialForecast, 'prefN')
            for i = 1:numel(w)
                switch w(i)
                    case 50
                        row.Npred_pref_50 = trialForecast.prefN(i);
                    case 125
                        row.Npred_pref_125 = trialForecast.prefN(i);
                    case 250
                        row.Npred_pref_250 = trialForecast.prefN(i);
                end
            end
        end

        if isfield(trialForecast, 'probeN')
            for i = 1:numel(w)
                switch w(i)
                    case 50
                        row.Npred_probe_50 = trialForecast.probeN(i);
                    case 125
                        row.Npred_probe_125 = trialForecast.probeN(i);
                    case 250
                        row.Npred_probe_250 = trialForecast.probeN(i);
                end
            end
        end
    end

    % useful debug / bookkeeping values
    if isfield(trialForecast, 'beta1')
        row.Npred_beta1 = trialForecast.beta1;
    end
    if isfield(trialForecast, 'pCorrect')
        row.Npred_pCorrect = trialForecast.pCorrect;
    end
    if isfield(trialForecast, 'stdSumDM')
        row.Npred_stdSumDM = trialForecast.stdSumDM;
    end
    if isfield(trialForecast, 'gamma')
        row.Npred_gamma = trialForecast.gamma;
    end
    if isfield(trialForecast, 'rPref')
        row.Npred_rPref = trialForecast.rPref;
    end
    if isfield(trialForecast, 'rProbe')
        row.Npred_rProbe = trialForecast.rProbe;
    end
    if isfield(trialForecast, 'refWindowMS')
        row.Npred_refWindowMS = trialForecast.refWindowMS;
    end
    if isfield(trialForecast, 'nRefFrames')
        row.Npred_nRefFrames = trialForecast.nRefFrames;
    end
    if isfield(trialForecast, 'zTarget')
        row.Npred_zTarget = trialForecast.zTarget;
    end
    if isfield(trialForecast, 'valid')
        row.Npred_valid = double(trialForecast.valid);
    end
    % optional bookkeeping
    if isfield(trialForecast, 'beta1')
        row.Npred_beta1 = trialForecast.beta1;
    end
    if isfield(trialForecast, 'pCorrect')
        row.Npred_pCorrect = trialForecast.pCorrect;
    end
end
end

end


function [b,p] = coefPv(mdl, name)
b = NaN; p = NaN;
if isempty(mdl) || ~isprop(mdl,'Coefficients')
    return;
end
T = mdl.Coefficients;
ix = strcmp(T.Properties.RowNames, name);
if any(ix)
    b = T.Estimate(ix);
    p = T.pValue(ix);
end
end

function b = coefPvOnly(mdl, name)
b = NaN;
if isempty(mdl) || ~isprop(mdl,'Coefficients')
    return;
end
T = mdl.Coefficients;
ix = strcmp(T.Properties.RowNames, name);
if any(ix)
    b = T.Estimate(ix);
end
end

function row = emptySummaryRow()
row = struct();

row.file = "";
row.session = "";
row.failed = false;
row.failMsg = "";

% ---- counts / performance ----
row.nTrials = NaN;
row.nInc = NaN;
row.nDec = NaN;
row.pCorrect = NaN;
row.pInc = NaN;
row.pDec = NaN;

% ---- descriptive correct-error differences ----
row.dMean_sumDM = NaN;
row.dMean_meanDM = NaN;
row.dMean_longestPosRun = NaN;
row.dMean_longestNegRun = NaN;
row.dMean_nSwitch = NaN;
row.dMean_dEarlyLate = NaN;
row.dMean_maxCumDM = NaN;
row.dMean_minCumDM = NaN;

% ---- model deviances ----
row.dev_m0 = NaN;
row.dev_m1 = NaN;
row.dev_m2 = NaN;
row.dev_m3 = NaN;
row.dev_m4 = NaN;

row.dDev_m1 = NaN;
row.dDev_m2 = NaN;
row.dDev_m3 = NaN;
row.dDev_m4 = NaN;

row.devExpl_m1 = NaN;
row.devExpl_m2 = NaN;
row.devExpl_m3 = NaN;
row.devExpl_m4 = NaN;

row.incOver_m1_m2 = NaN;
row.incOver_m1_m3 = NaN;
row.incOver_m1_m4 = NaN;

% ---- coefficients: m1 ----
row.b_m1_sumDM = NaN;
row.p_m1_sumDM = NaN;

% ---- coefficients: m2 ----
row.b_m2_sumDM = NaN;
row.p_m2_sumDM = NaN;
row.b_m2_longestPosRun = NaN;
row.p_m2_longestPosRun = NaN;
row.b_m2_longestNegRun = NaN;
row.p_m2_longestNegRun = NaN;
row.b_m2_nSwitch = NaN;
row.p_m2_nSwitch = NaN;

% ---- coefficients: m3 ----
row.b_m3_sumDM = NaN;
row.p_m3_sumDM = NaN;
row.b_m3_dEarlyLate = NaN;
row.p_m3_dEarlyLate = NaN;

% ---- coefficients: m4 ----
row.b_m4_sumDM = NaN;
row.p_m4_sumDM = NaN;
row.b_m4_maxCumDM = NaN;
row.p_m4_maxCumDM = NaN;
row.b_m4_minCumDM = NaN;
row.p_m4_minCumDM = NaN;

% ---- trial forecasts: preferred ----
row.Npred_pref_50 = NaN;
row.Npred_pref_125 = NaN;
row.Npred_pref_250 = NaN;

% ---- trial forecasts: probe ----
row.Npred_probe_50 = NaN;
row.Npred_probe_125 = NaN;
row.Npred_probe_250 = NaN;

% ---- forecast debug / bookkeeping ----
row.Npred_beta1 = NaN;
row.Npred_pCorrect = NaN;
row.Npred_stdSumDM = NaN;
row.Npred_gamma = NaN;
row.Npred_rPref = NaN;
row.Npred_rProbe = NaN;
row.Npred_refWindowMS = NaN;
row.Npred_nRefFrames = NaN;
row.Npred_zTarget = NaN;
row.Npred_valid = NaN;


end

% ===================== SAFE GLM =====================

function mdl = safeFitglm(tbl, formula, varargin)
% Remove predictors that are constant within this subset to avoid rank-deficient warnings.

parts = split(strrep(formula,' ',''), '~');
yName = parts{1};
rhs   = parts{2};

if strcmp(rhs,'1')
    mdl = fitglm(tbl, formula, varargin{:});
    return;
end

preds = split(rhs, '+');

keep = strings(0);
for i = 1:numel(preds)
    p = preds{i};
    if strcmp(p,'1') || strcmp(p,'0'), continue; end
    if ~ismember(p, tbl.Properties.VariableNames), continue; end

    v = tbl.(p);
    if islogical(v), v = double(v); end
    if all(isnan(v)), continue; end
    if nanstd(v) < eps, continue; end

    keep(end+1) = p; %#ok<AGROW>
end

if isempty(keep)
    newFormula = yName + "~1";
else
    newFormula = yName + "~" + strjoin(keep, "+");
end

mdl = fitglm(tbl, char(newFormula), varargin{:});
end

% ===================== HELPERS =====================

function x = fillFrames(valuesPC, timesMS, m, msPerVFrame)
x = nan(m,1);
nTimes = length(timesMS);
for tIndex = 1:nTimes
    if tIndex == 1
        theVFrame = 1;
    else
        theVFrame = max(1, round(timesMS(tIndex) / msPerVFrame));
    end
    if tIndex < nTimes
        nextVFrame = round(timesMS(tIndex + 1) / msPerVFrame);
    else
        nextVFrame = m + 1;
    end
    v0 = max(1, theVFrame);
    v1 = min(m, nextVFrame-1);
    if v1 >= v0
        x(v0:v1) = valuesPC(tIndex);
    end
end
end

function tf = isNoNoise(x)
tf = isscalar(x) && isequal(x, 0);
end

function v = getScalarField(tr, name)
v = NaN;
if isfield(tr, name)
    f = tr.(name);
    if isstruct(f) && isfield(f,'data') && ~isempty(f.data)
        v = f.data(1);
    elseif isnumeric(f) && ~isempty(f)
        v = f(1);
    end
end
end

% function [b, s] = getBaseStep(tr)
% % Try multiple likely field locations for base and step coherence (PC).
% b = NaN; s = NaN;
% 
% b = firstNonNan(b, getScalarField(tr, 'baseCohPC'));
% s = firstNonNan(s, getScalarField(tr, 'stepCohPC'));
% 
% if isfield(tr, 'trial')
%     b = firstNonNan(b, getScalarField(tr.trial, 'baseCohPC'));
%     s = firstNonNan(s, getScalarField(tr.trial, 'stepCohPC'));
%     if isfield(tr.trial, 'data')
%         b = firstNonNan(b, getScalarField(tr.trial.data, 'baseCohPC'));
%         s = firstNonNan(s, getScalarField(tr.trial.data, 'stepCohPC'));
%     end
% end
% 
% if isfield(tr,'blockStatus') && isfield(tr.blockStatus,'data')
%     bs = tr.blockStatus.data;
%     b = firstNonNan(b, getScalarField(bs,'baseCohPC'));
%     s = firstNonNan(s, getScalarField(bs,'stepCohPC'));
% end
% end

% function v = firstNonNan(v, cand)
% if isnan(v) && ~isnan(cand)
%     v = cand;
% end
% end

% function stratPlot(y, x, xlab)
% ok = ~isnan(x) & ~isnan(y);
% x = x(ok); y = y(ok);
% 
% if numel(x) < 10
%     plot(nan, nan); title([xlab ' (too few points)']); return;
% end
% 
% q = quantile(x, [0 .25 .5 .75 1]);
% qU = unique(q);
% 
% if numel(qU) < 3
%     u = unique(x);
%     if numel(u) > 6
%         idx = round(linspace(1,numel(u),6));
%         u = u(idx);
%     end
%     edges = [-inf; u(:); inf];
% else
%     edges = qU(:); edges(1) = -inf; edges(end) = inf;
% end
% 
% bin = discretize(x, edges);
% K = max(bin);
% 
% p = accumarray(bin, y, [K 1], @mean, NaN);
% n = accumarray(bin, y, [K 1], @numel, NaN);
% 
% plot(1:K, p, '-o'); ylim([0 1]); xlim([0.5 K+0.5]);
% set(gca,'XTick',1:K);
% xlabel('Bin'); ylabel('P(correct)');
% title(xlab);
% for k = 1:K
%     if ~isnan(p(k))
%         text(k, p(k), sprintf(' n=%d', n(k)), 'VerticalAlignment','bottom');
%     end
% end
% end

function excl = defaultExcluded()
excl = {'IDReadout_Meetz_20260113', 'IDReadout_Meetz_20260114', 'IDReadout_Meetz_20260114_2', ...
    'IDReadout_Meetz_20260114_3', 'IDReadout_Meetz_20260114_4', 'IDReadout_Meetz_20260115', ...
    'IDReadout_Meetz_20260116', 'IDReadout_Meetz_20260209', 'IDReadout_Meetz_20260210'};
end

function folder = defaultDataFolder()
% Try to mimic your kernelAverage behavior if dataFolderPath() exists.
folder = pwd;
try
    if exist('dataFolderPath','file') == 2
        folder = char(dataFolderPath());
    end
catch
end
folder = fullfile(folder, 'Data');
end

function batchPlots(summary)
% Simple batch plots if you want them.
if isempty(summary) || ~istable(summary), return; end
if any(summary.failed)
    fprintf('\nSome sessions failed; see summary.failMsg.\n');
end

% Example: deviance improvements
figure('Name','Batch: deviance improvements','Color','w');
plot(summary.dDevDec_m2, '-o'); hold on;
plot(summary.dDevInc_m2, '-o');
legend('DEC dDev (m0-m2)','INC dDev (m0-m2)');
xlabel('Session index'); ylabel('Deviance improvement');

figure('Name','Batch: boundary-hit rates','Color','w');
plot(summary.pLoCh, '-o'); hold on;
plot(summary.pHiCh, '-o');
legend('DEC P(countLoCh>0)','INC P(countHiCh>0)');
xlabel('Session index'); ylabel('Probability');
end

function L = longestRunPerTrial(mask)
% mask: nFrames x nTrials logical
% L   : nTrials x 1, longest consecutive run of true values in each trial

[~, nTrials] = size(mask);
L = zeros(nTrials,1);

for t = 1:nTrials
    x = mask(:,t);
    if ~any(x)
        L(t) = 0;
        continue;
    end
    dx = diff([false; x; false]);
    runStarts = find(dx == 1);
    runEnds   = find(dx == -1) - 1;
    L(t) = max(runEnds - runStarts + 1);
end
end

% function localDiagScalar(name, x, y)
% x = x(:);
% y = double(y(:));
% 
% ok = ~isnan(x) & ~isnan(y);
% x = x(ok);
% y = y(ok);
% 
% if numel(x) < 10
%     fprintf('%s: too few valid points\n', name);
%     return;
% end
% 
% pPos = NaN;
% pNeg = NaN;
% if any(x > 0), pPos = mean(y(x > 0)); end
% if any(x < 0), pNeg = mean(y(x < 0)); end
% 
% r = corr(x, y, 'rows', 'complete');
% 
% tbl = table(y, x, 'VariableNames', {'isCorrect','x'});
% m0 = fitglm(tbl, 'isCorrect ~ 1', 'Distribution', 'binomial');
% m1 = fitglm(tbl, 'isCorrect ~ x', 'Distribution', 'binomial');
% devExpl = (m0.Deviance - m1.Deviance) / m0.Deviance;
% 
% fprintf('%-12s  P(corr|x>0)=%.4f  P(corr|x<0)=%.4f  corr=%.4f  devExpl=%.6f\n', ...
%     name, pPos, pNeg, r, devExpl);
% end

% function localPvPlot(figNum, x, y, plotTitle)
% x = x(:);
% y = double(y(:));
% 
% ok = ~isnan(x) & ~isnan(y);
% x = x(ok);
% y = y(ok);
% 
% if numel(x) < 20
%     return;
% end
% 
% edges = quantile(x, 0:0.1:1);
% edges = unique(edges);
% 
% if numel(edges) < 3
%     return;
% end
% 
% nBins = numel(edges) - 1;
% xb = nan(nBins,1);
% pb = nan(nBins,1);
% nb = nan(nBins,1);
% 
% for b = 1:nBins
%     if b < nBins
%         idx = x >= edges(b) & x < edges(b+1);
%     else
%         idx = x >= edges(b) & x <= edges(b+1);
%     end
%     if any(idx)
%         xb(b) = mean(x(idx));
%         pb(b) = mean(y(idx));
%         nb(b) = sum(idx);
%     end
% end
% 
% figure(figNum); clf;
% plot(xb, pb, 'o-');
% ylim([0 1]);
% xlabel('Scalar evidence');
% ylabel('P(correct)');
% title(plotTitle);
% grid on;
% 
% for b = 1:nBins
%     if ~isnan(xb(b))
%         text(xb(b), pb(b), sprintf(' n=%d', nb(b)), ...
%             'VerticalAlignment', 'bottom');
%     end
% end
% end