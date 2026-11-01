function results = refitIDRProbeBiasGain(fitPath, varargin)
%REFITIDRPROBEBIASGAIN Refit IDR probe sessions with bias and gain terms.
%
% results = refitIDRProbeBiasGain(fitPath)
%
% Starting from IDRGainProbeOffsetFits_<Animal>.mat, this function refits each
% included probeSession using four nested models:
%
%   null:       deltaC_eff = max(deltaC + gPref*nPref, 0)
%   gain-only:  deltaC_eff = max(deltaC + gPref*nPref + gProbe*nProbe, 0)
%   bias-only:  deltaC_eff = max(deltaC + bProbe + gPref*nPref, 0)
%   joint:      deltaC_eff = max(deltaC + bProbe + gPref*nPref ...
%                               + gProbe*nProbe, 0)
%
% P(correct) is the same Weibull used by the source fit. Parent alpha, parent
% preferred gain, shared beta, and shared lapse remain fixed. Positive bProbe
% denotes added coherence-equivalent evidence and therefore a lower nominal
% threshold. The gain-only model is recomputed as an internal check against the
% stored fit.
%
% Name-value options
%   'OutputFolder'      output folder (default: source MAT folder)
%   'OutputStem'        output filename stem (default derived from source)
%   'GainLimit'         absolute bound on gProbe (default 5)
%   'BiasLimit'         absolute bProbe bound, percent coherence (default 30)
%   'Combine179And180'  pool 180 deg with 179 deg (default true)
%   'NumBootstrap'      parent-session bootstrap samples (default 5000)
%   'RandomSeed'        bootstrap seed (default 1)
%   'MakePlot'          export summary PDF (default true)
%   'SaveOutputs'       save MAT and CSV outputs (default true)
%   'Verbose'           print progress (default true)

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'fitPath', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'OutputFolder', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'OutputStem', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'GainLimit', 5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'BiasLimit', 30, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'Combine179And180', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'NumBootstrap', 5000, @(x) isnumeric(x) && isscalar(x) && x >= 0 && mod(x,1) == 0);
addParameter(p, 'RandomSeed', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'MakePlot', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'SaveOutputs', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) && isscalar(x));
parse(p, fitPath, varargin{:});
opts = p.Results;

fitPath = char(fitPath);
if ~isfile(fitPath)
  error('refitIDRProbeBiasGain:MissingInput', 'Input MAT file not found: %s', fitPath);
end
S = load(fitPath, 'fitData');
if ~isfield(S, 'fitData') || ~isfield(S.fitData, 'probeTable') || ...
    ~isfield(S.fitData, 'trialTable')
  error('refitIDRProbeBiasGain:MissingData', ...
    'Input must contain fitData.probeTable and fitData.trialTable.');
end
P = S.fitData.probeTable;
T = S.fitData.trialTable;

probeRequired = {'probeSessionIndex','probeSessionName','parentSessionName', ...
  'probeDirDeg','alphaParent','gPrefParent','betaShared','lapseShared', ...
  'nIncTrials','fitIncluded','gProbe','nllProbe'};
trialRequired = {'probeSessionIndex','deltaC','nPref','nProbe','correct'};
missingProbe = probeRequired(~ismember(probeRequired, P.Properties.VariableNames));
missingTrial = trialRequired(~ismember(trialRequired, T.Properties.VariableNames));
if ~isempty(missingProbe) || ~isempty(missingTrial)
  error('refitIDRProbeBiasGain:MissingVariables', ...
    'Missing probe variables [%s]; missing trial variables [%s].', ...
    strjoin(missingProbe, ', '), strjoin(missingTrial, ', '));
end

fitMask = logical(P.fitIncluded) & isfinite(P.gProbe) & isfinite(P.alphaParent) & ...
  P.alphaParent > 0 & isfinite(P.gPrefParent) & isfinite(P.betaShared) & ...
  P.betaShared > 0 & isfinite(P.lapseShared) & P.lapseShared >= 0 & P.lapseShared < 0.25;
sourceRows = find(fitMask);
if isempty(sourceRows)
  error('refitIDRProbeBiasGain:NoFits', 'No valid fitted probe sessions were available.');
end

nFits = numel(sourceRows);
rows = repmat(emptyResultRow(), nFits, 1);
optimOptions = optimset('Display', 'off', 'MaxIter', 2000, ...
  'MaxFunEvals', 8000, 'TolX', 1e-9, 'TolFun', 1e-9);

for iFit = 1:nFits
  iSource = sourceRows(iFit);
  probeIndex = double(P.probeSessionIndex(iSource));
  trialMask = double(T.probeSessionIndex) == probeIndex;
  if ~any(trialMask)
    error('refitIDRProbeBiasGain:MissingTrials', ...
      'No trial rows matched probeSessionIndex %g.', probeIndex);
  end

  D = struct();
  D.deltaC = double(T.deltaC(trialMask));
  D.nPref = double(T.nPref(trialMask));
  D.nProbe = double(T.nProbe(trialMask));
  D.correct = double(T.correct(trialMask));
  D.alpha = double(P.alphaParent(iSource));
  D.gPref = double(P.gPrefParent(iSource));
  D.beta = double(P.betaShared(iSource));
  D.lapse = double(P.lapseShared(iSource));

  oldGain = double(P.gProbe(iSource));
  nullNLL = modelNLLPhysical(0, 0, D);
  gainFit = fitGainOnly(D, oldGain, opts.GainLimit, optimOptions);
  biasInitial = initialBiasFromMean(D, opts.BiasLimit);
  biasFit = fitBiasOnly(D, biasInitial, opts.BiasLimit, optimOptions);
  jointFit = fitJoint(D, biasFit.bProbe, gainFit.gProbe, ...
    opts.BiasLimit, opts.GainLimit, optimOptions);
  localCovariance = localBiasGainCovariance(D, jointFit.bProbe, jointFit.gProbe, ...
    opts.BiasLimit, opts.GainLimit);

  r = emptyResultRow();
  r.sourceProbeTableRow = iSource;
  r.probeSessionIndex = probeIndex;
  r.probeSessionName = string(P.probeSessionName(iSource));
  r.parentSessionName = string(P.parentSessionName(iSource));
  r.probeDirDeg = double(P.probeDirDeg(iSource));
  if opts.Combine179And180 && r.probeDirDeg == 180
    r.probeDirDeg = 179;
  end
  r.nTrials = sum(trialMask);
  r.alphaParent = D.alpha;
  r.gPrefParent = D.gPref;
  r.betaShared = D.beta;
  r.lapseShared = D.lapse;
  r.parentC75 = parentC75Local(D.alpha, D.beta, D.lapse);
  r.gProbeStored = oldGain;
  r.gProbeStoredOverGPref = oldGain ./ D.gPref;
  r.nllGainStored = double(P.nllProbe(iSource));
  r.nllNull = nullNLL;
  r.bProbeBiasOnly = biasFit.bProbe;
  r.nllBiasOnly = biasFit.nll;
  r.exitflagBiasOnly = biasFit.exitflag;
  r.gProbeGainOnly = gainFit.gProbe;
  r.gProbeGainOnlyOverGPref = gainFit.gProbe ./ D.gPref;
  r.nllGainOnly = gainFit.nll;
  r.exitflagGainOnly = gainFit.exitflag;
  r.bProbeJoint = jointFit.bProbe;
  r.gProbeJoint = jointFit.gProbe;
  r.gProbeJointOverGPref = jointFit.gProbe ./ D.gPref;
  r.nllJoint = jointFit.nll;
  r.exitflagJoint = jointFit.exitflag;
  r.biasEquivalentThresholdRatio = (r.parentC75 - jointFit.bProbe) ./ r.parentC75;
  r.deltaGainJointMinusStored = jointFit.gProbe - oldGain;
  r.deltaGainRatioJointMinusStored = r.gProbeJointOverGPref - r.gProbeStoredOverGPref;
  r.deltaNLLBiasAfterGain = max(gainFit.nll - jointFit.nll, 0);
  r.pBiasAfterGain = chiSquareUpperLocal(2 .* r.deltaNLLBiasAfterGain, 1);
  r.deltaNLLGainAfterBias = max(biasFit.nll - jointFit.nll, 0);
  r.pGainAfterBias = chiSquareUpperLocal(2 .* r.deltaNLLGainAfterBias, 1);
  r.deltaNLLJointVsNull = max(nullNLL - jointFit.nll, 0);
  r.pJointVsNull = chiSquareUpperLocal(2 .* r.deltaNLLJointVsNull, 2);
  r.localSEBProbe = localCovariance.seB;
  r.localSEGProbe = localCovariance.seG;
  r.localCorrBiasGain = localCovariance.corrBG;
  r.biasNearBound = abs(jointFit.bProbe) > 0.95 .* opts.BiasLimit;
  r.gainNearBound = abs(jointFit.gProbe) > 0.95 .* opts.GainLimit;
  r.storedGainNLLDifference = gainFit.nll - r.nllGainStored;
  rows(iFit) = r;

  if opts.Verbose && (mod(iFit, 20) == 0 || iFit == nFits)
    fprintf('Refit %d/%d probe sessions\n', iFit, nFits);
  end
end

refitTable = struct2table(rows);
[offsetSummary, bootstrap] = summarizeOffsetsWithParentBootstrap(refitTable, ...
  opts.NumBootstrap, opts.RandomSeed);

population = struct();
population.nProbeSessions = height(refitTable);
population.nParentSessions = numel(unique(refitTable.parentSessionName));
population.nOffsets = height(offsetSummary);
population.totalDeltaNLLBiasAfterGain = sum(refitTable.deltaNLLBiasAfterGain);
population.totalDeltaNLLGainAfterBias = sum(refitTable.deltaNLLGainAfterBias);
population.pAggregateBiasAfterGain = chiSquareUpperLocal( ...
  2 .* population.totalDeltaNLLBiasAfterGain, population.nProbeSessions);
population.pAggregateGainAfterBias = chiSquareUpperLocal( ...
  2 .* population.totalDeltaNLLGainAfterBias, population.nProbeSessions);
population.nBiasAfterGainP005 = sum(refitTable.pBiasAfterGain < 0.05);
population.nGainAfterBiasP005 = sum(refitTable.pGainAfterBias < 0.05);
population.medianLocalCorrBiasGain = median(refitTable.localCorrBiasGain, 'omitnan');
population.corrBiasGainAcrossSessions = correlationLocal( ...
  refitTable.bProbeJoint, refitTable.gProbeJointOverGPref);
population.maxAbsStoredGainNLLDifference = max(abs(refitTable.storedGainNLLDifference));
population.nBiasNearBound = sum(refitTable.biasNearBound);
population.nGainNearBound = sum(refitTable.gainNearBound);

results = struct();
results.version = 1;
results.refitTable = refitTable;
results.offsetSummary = offsetSummary;
results.bootstrap = bootstrap;
results.population = population;
results.options = opts;
results.sourceFitPath = fitPath;
results.method = ['Per-probeSession nested null, gain-only, bias-only, and joint bias+gain fits. ' ...
  'Parent alpha/gPref and shared beta/lapse fixed from source fit. Offset summaries use ' ...
  'parent-session cluster bootstrap without refitting.'];
results.model = ['deltaC_eff = max(deltaC + bProbe + gPrefParent*nPref ' ...
  '+ gProbe*nProbe, 0)'];

[sourceFolder, sourceStem] = fileparts(fitPath);
if isempty(opts.OutputFolder)
  outputFolder = sourceFolder;
else
  outputFolder = char(opts.OutputFolder);
end
if isempty(outputFolder)
  outputFolder = pwd;
end
if isempty(opts.OutputStem)
  outputStem = [sourceStem '_BiasGainRefit'];
else
  outputStem = char(opts.OutputStem);
end
if (opts.MakePlot || opts.SaveOutputs) && ~isfolder(outputFolder)
  mkdir(outputFolder);
end

results.outputPaths = struct('pdf', '', 'mat', '', 'refitCSV', '', 'offsetSummaryCSV', '');
if opts.MakePlot
  pdfPath = fullfile(outputFolder, [outputStem '.pdf']);
  makeSummaryPlot(refitTable, offsetSummary, population, opts, pdfPath);
  results.outputPaths.pdf = pdfPath;
end
if opts.SaveOutputs
  matPath = fullfile(outputFolder, [outputStem '.mat']);
  refitCSVPath = fullfile(outputFolder, [outputStem '_ProbeSessions.csv']);
  offsetCSVPath = fullfile(outputFolder, [outputStem '_OffsetSummary.csv']);
  results.outputPaths.mat = matPath;
  results.outputPaths.refitCSV = refitCSVPath;
  results.outputPaths.offsetSummaryCSV = offsetCSVPath;
  writetable(refitTable, refitCSVPath);
  writetable(offsetSummary, offsetCSVPath);
  save(matPath, 'results', '-v7.3');
end

if opts.Verbose
  fprintf('\nIDR probe bias+gain refit:\n');
  fprintf('  probe sessions:                    %d\n', population.nProbeSessions);
  fprintf('  total Delta NLL, bias after gain:  %.4g\n', population.totalDeltaNLLBiasAfterGain);
  fprintf('  aggregate bias-after-gain p:       %.4g\n', population.pAggregateBiasAfterGain);
  fprintf('  total Delta NLL, gain after bias:  %.4g\n', population.totalDeltaNLLGainAfterBias);
  fprintf('  aggregate gain-after-bias p:       %.4g\n', population.pAggregateGainAfterBias);
  fprintf('  bias after gain p<0.05:            %d/%d\n', ...
    population.nBiasAfterGainP005, population.nProbeSessions);
  fprintf('  gain after bias p<0.05:            %d/%d\n', ...
    population.nGainAfterBiasP005, population.nProbeSessions);
  fprintf('  median local corr(b,g):            %.4g\n', population.medianLocalCorrBiasGain);
  if opts.MakePlot, fprintf('Saved PDF: %s\n', results.outputPaths.pdf); end
  if opts.SaveOutputs, fprintf('Saved MAT: %s\n', results.outputPaths.mat); end
end
end

% -------------------------------------------------------------------------
function row = emptyResultRow()
names = {'sourceProbeTableRow','probeSessionIndex','probeDirDeg','nTrials', ...
  'alphaParent','gPrefParent','betaShared','lapseShared','parentC75', ...
  'gProbeStored','gProbeStoredOverGPref','nllGainStored','nllNull', ...
  'bProbeBiasOnly','nllBiasOnly','exitflagBiasOnly','gProbeGainOnly', ...
  'gProbeGainOnlyOverGPref','nllGainOnly','exitflagGainOnly','bProbeJoint', ...
  'gProbeJoint','gProbeJointOverGPref','nllJoint','exitflagJoint', ...
  'biasEquivalentThresholdRatio','deltaGainJointMinusStored', ...
  'deltaGainRatioJointMinusStored','deltaNLLBiasAfterGain','pBiasAfterGain', ...
  'deltaNLLGainAfterBias','pGainAfterBias','deltaNLLJointVsNull','pJointVsNull', ...
  'localSEBProbe','localSEGProbe','localCorrBiasGain','storedGainNLLDifference'};
row = struct();
for i = 1:numel(names), row.(names{i}) = NaN; end
row.probeSessionName = "";
row.parentSessionName = "";
row.biasNearBound = false;
row.gainNearBound = false;
end

% -------------------------------------------------------------------------
function fit = fitGainOnly(D, initialGain, gainLimit, options)
starts = unique([initialGain; 0; 0.5; -0.5]);
best = struct('gProbe', NaN, 'nll', Inf, 'exitflag', NaN);
for i = 1:numel(starts)
  theta0 = boundedToTheta(starts(i), gainLimit);
  objective = @(theta) modelNLLPhysical(0, gainLimit .* tanh(theta(1)), D);
  [theta, nll, exitflag] = fminsearch(objective, theta0, options);
  if nll < best.nll
    best.gProbe = gainLimit .* tanh(theta(1));
    best.nll = nll;
    best.exitflag = exitflag;
  end
end
fit = best;
end

% -------------------------------------------------------------------------
function fit = fitBiasOnly(D, initialBias, biasLimit, options)
starts = unique([initialBias; 0; 2; -2]);
best = struct('bProbe', NaN, 'nll', Inf, 'exitflag', NaN);
for i = 1:numel(starts)
  theta0 = boundedToTheta(starts(i), biasLimit);
  objective = @(theta) modelNLLPhysical(biasLimit .* tanh(theta(1)), 0, D);
  [theta, nll, exitflag] = fminsearch(objective, theta0, options);
  if nll < best.nll
    best.bProbe = biasLimit .* tanh(theta(1));
    best.nll = nll;
    best.exitflag = exitflag;
  end
end
fit = best;
end

% -------------------------------------------------------------------------
function fit = fitJoint(D, bInitial, gInitial, biasLimit, gainLimit, options)
starts = [bInitial gInitial; bInitial 0; 0 gInitial; 0 0; 2 0.5; -2 0.5];
best = struct('bProbe', NaN, 'gProbe', NaN, 'nll', Inf, 'exitflag', NaN);
for i = 1:size(starts,1)
  theta0 = [boundedToTheta(starts(i,1), biasLimit), ...
    boundedToTheta(starts(i,2), gainLimit)];
  objective = @(theta) modelNLLPhysical(biasLimit .* tanh(theta(1)), ...
    gainLimit .* tanh(theta(2)), D);
  [theta, nll, exitflag] = fminsearch(objective, theta0, options);
  if nll < best.nll
    best.bProbe = biasLimit .* tanh(theta(1));
    best.gProbe = gainLimit .* tanh(theta(2));
    best.nll = nll;
    best.exitflag = exitflag;
  end
end
fit = best;
end

% -------------------------------------------------------------------------
function theta = boundedToTheta(value, limit)
z = max(min(value ./ limit, 0.95), -0.95);
theta = atanh(z);
end

% -------------------------------------------------------------------------
function b0 = initialBiasFromMean(D, biasLimit)
pObserved = mean(D.correct);
margin = 0.5 ./ (numel(D.correct) + 1);
pObserved = min(max(pObserved, 0.5 + margin), 1 - D.lapse - margin);
z = -log(1 - (pObserved - 0.5) ./ (0.5 - D.lapse));
equivalentC = D.alpha .* z .^ (1 ./ D.beta);
b0 = equivalentC - mean(D.deltaC);
b0 = max(min(b0, 0.8 .* biasLimit), -0.8 .* biasLimit);
end

% -------------------------------------------------------------------------
function nll = modelNLLPhysical(bProbe, gProbe, D)
deltaCEff = max(D.deltaC + bProbe + D.gPref .* D.nPref + gProbe .* D.nProbe, 0);
p = 0.5 + (0.5 - D.lapse) .* (1 - exp(-((deltaCEff ./ D.alpha) .^ D.beta)));
p = min(max(p, 1e-12), 1 - 1e-12);
nll = -sum(D.correct .* log(p) + (1 - D.correct) .* log(1 - p));
if ~isfinite(nll), nll = realmax; end
end

% -------------------------------------------------------------------------
function C = localBiasGainCovariance(D, b, g, biasLimit, gainLimit)
C = struct('seB', NaN, 'seG', NaN, 'corrBG', NaN);
if abs(b) > 0.9 .* biasLimit || abs(g) > 0.9 .* gainLimit
  return;
end
hB = max(0.02, 0.002 .* max(1, abs(b)));
hG = max(0.002, 0.002 .* max(1, abs(g)));
f00 = modelNLLPhysical(b, g, D);
fBp = modelNLLPhysical(b+hB, g, D);
fBm = modelNLLPhysical(b-hB, g, D);
fGp = modelNLLPhysical(b, g+hG, D);
fGm = modelNLLPhysical(b, g-hG, D);
fPP = modelNLLPhysical(b+hB, g+hG, D);
fPM = modelNLLPhysical(b+hB, g-hG, D);
fMP = modelNLLPhysical(b-hB, g+hG, D);
fMM = modelNLLPhysical(b-hB, g-hG, D);
H = [(fBp - 2*f00 + fBm)./hB.^2, (fPP-fPM-fMP+fMM)./(4*hB*hG); ...
     (fPP-fPM-fMP+fMM)./(4*hB*hG), (fGp - 2*f00 + fGm)./hG.^2];
if all(isfinite(H), 'all') && all(eig((H+H')./2) > 0)
  covariance = inv((H+H')./2);
  if all(diag(covariance) > 0)
    C.seB = sqrt(covariance(1,1));
    C.seG = sqrt(covariance(2,2));
    C.corrBG = covariance(1,2) ./ (C.seB .* C.seG);
  end
end
end

% -------------------------------------------------------------------------
function c75 = parentC75Local(alpha, beta, lapse)
z75 = -log(1 - 0.25 ./ (0.5 - lapse));
c75 = alpha .* z75 .^ (1 ./ beta);
end

% -------------------------------------------------------------------------
function p = chiSquareUpperLocal(x, df)
p = gammainc(max(x,0)./2, df./2, 'upper');
end

% -------------------------------------------------------------------------
function r = correlationLocal(x, y)
m = isfinite(x) & isfinite(y);
if sum(m) < 3
  r = NaN;
else
  c = corrcoef(x(m), y(m));
  r = c(1,2);
end
end

% -------------------------------------------------------------------------
function [summary, bootstrap] = summarizeOffsetsWithParentBootstrap(T, nBoot, seed)
offsets = sort(unique(T.probeDirDeg));
nOffsets = numel(offsets);
variables = {'probeDirDeg','nProbeSessions','nParents','medianBProbe', ...
  'bProbeCILow','bProbeCIHigh','medianBiasThresholdRatio', ...
  'biasThresholdRatioCILow','biasThresholdRatioCIHigh','medianStoredGainRatio', ...
  'storedGainRatioCILow','storedGainRatioCIHigh','medianJointGainRatio', ...
  'jointGainRatioCILow','jointGainRatioCIHigh','medianGainRatioChange', ...
  'gainRatioChangeCILow','gainRatioChangeCIHigh','medianDeltaNLLBiasAfterGain', ...
  'fractionBiasAfterGainP005','medianDeltaNLLGainAfterBias','fractionGainAfterBiasP005'};
summary = array2table(nan(nOffsets, numel(variables)), 'VariableNames', variables);

parents = unique(T.parentSessionName, 'stable');
nParents = numel(parents);
bootB = nan(nBoot,nOffsets);
bootTR = nan(nBoot,nOffsets);
bootOldG = nan(nBoot,nOffsets);
bootNewG = nan(nBoot,nOffsets);
bootDG = nan(nBoot,nOffsets);
oldRng = rng;
cleanupRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(seed, 'twister');
for iBoot = 1:nBoot
  sampled = parents(randi(nParents,nParents,1));
  sets = cell(nParents,1);
  for iParent = 1:nParents
    sets{iParent} = find(T.parentSessionName == sampled(iParent));
  end
  idx = vertcat(sets{:});
  for iOffset = 1:nOffsets
    m = T.probeDirDeg(idx) == offsets(iOffset);
    if any(m)
      ii = idx(m);
      bootB(iBoot,iOffset) = median(T.bProbeJoint(ii), 'omitnan');
      bootTR(iBoot,iOffset) = median(T.biasEquivalentThresholdRatio(ii), 'omitnan');
      bootOldG(iBoot,iOffset) = median(T.gProbeStoredOverGPref(ii), 'omitnan');
      bootNewG(iBoot,iOffset) = median(T.gProbeJointOverGPref(ii), 'omitnan');
      bootDG(iBoot,iOffset) = median(T.deltaGainRatioJointMinusStored(ii), 'omitnan');
    end
  end
end

for iOffset = 1:nOffsets
  m = T.probeDirDeg == offsets(iOffset);
  summary.probeDirDeg(iOffset) = offsets(iOffset);
  summary.nProbeSessions(iOffset) = sum(m);
  summary.nParents(iOffset) = numel(unique(T.parentSessionName(m)));
  summary.medianBProbe(iOffset) = median(T.bProbeJoint(m), 'omitnan');
  summary.medianBiasThresholdRatio(iOffset) = median(T.biasEquivalentThresholdRatio(m), 'omitnan');
  summary.medianStoredGainRatio(iOffset) = median(T.gProbeStoredOverGPref(m), 'omitnan');
  summary.medianJointGainRatio(iOffset) = median(T.gProbeJointOverGPref(m), 'omitnan');
  summary.medianGainRatioChange(iOffset) = median(T.deltaGainRatioJointMinusStored(m), 'omitnan');
  summary.medianDeltaNLLBiasAfterGain(iOffset) = median(T.deltaNLLBiasAfterGain(m), 'omitnan');
  summary.fractionBiasAfterGainP005(iOffset) = mean(T.pBiasAfterGain(m) < 0.05);
  summary.medianDeltaNLLGainAfterBias(iOffset) = median(T.deltaNLLGainAfterBias(m), 'omitnan');
  summary.fractionGainAfterBiasP005(iOffset) = mean(T.pGainAfterBias(m) < 0.05);
  summary{iOffset, {'bProbeCILow','bProbeCIHigh'}} = percentileLocal(bootB(:,iOffset), [2.5 97.5]);
  summary{iOffset, {'biasThresholdRatioCILow','biasThresholdRatioCIHigh'}} = percentileLocal(bootTR(:,iOffset), [2.5 97.5]);
  summary{iOffset, {'storedGainRatioCILow','storedGainRatioCIHigh'}} = percentileLocal(bootOldG(:,iOffset), [2.5 97.5]);
  summary{iOffset, {'jointGainRatioCILow','jointGainRatioCIHigh'}} = percentileLocal(bootNewG(:,iOffset), [2.5 97.5]);
  summary{iOffset, {'gainRatioChangeCILow','gainRatioChangeCIHigh'}} = percentileLocal(bootDG(:,iOffset), [2.5 97.5]);
end
bootstrap = struct('bProbe', bootB, 'biasThresholdRatio', bootTR, ...
  'storedGainRatio', bootOldG, 'jointGainRatio', bootNewG, ...
  'gainRatioChange', bootDG, 'offsets', offsets, 'nBootstrap', nBoot, ...
  'resamplingUnit', 'parentSession');
end

% -------------------------------------------------------------------------
function q = percentileLocal(x, pct)
x = sort(x(isfinite(x)));
if isempty(x), q = nan(size(pct)); return; end
if numel(x) == 1, q = repmat(x,size(pct)); return; end
position = 1 + (numel(x)-1) .* pct ./ 100;
lo = floor(position); hi = ceil(position); weight = position-lo;
% q = x(lo).*(1-weight) + x(hi).*weight;
q = reshape(x(lo),1,[]).*(1-weight) + ...
  reshape(x(hi),1,[]).*weight;
end

% -------------------------------------------------------------------------
function makeSummaryPlot(T, S, population, opts, pdfPath)
offsets = S.probeDirDeg(:)';
nOffsets = numel(offsets);
fig = figure('Color','w','Units','inches','Position',[1 1 14 8]);
layout = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
title(layout,'IDR probe bias + gain refit','FontWeight','bold');

ax = nexttile(layout,1); hold(ax,'on');
errorbar(ax,1:nOffsets,S.medianBProbe,S.medianBProbe-S.bProbeCILow, ...
  S.bProbeCIHigh-S.medianBProbe,'ko-','MarkerFaceColor','k','LineWidth',1.4);
yline(ax,0,':'); formatOffsetAxis(ax,offsets);
ylabel(ax,'b_{probe} (% coherence)'); title(ax,'Mean evidence term'); grid(ax,'on');

ax = nexttile(layout,2); hold(ax,'on');
errorbar(ax,1:nOffsets,S.medianStoredGainRatio, ...
  S.medianStoredGainRatio-S.storedGainRatioCILow, ...
  S.storedGainRatioCIHigh-S.medianStoredGainRatio,'o-','Color',[.55 .55 .55], ...
  'MarkerFaceColor',[.75 .75 .75],'LineWidth',1.2);
errorbar(ax,1:nOffsets,S.medianJointGainRatio, ...
  S.medianJointGainRatio-S.jointGainRatioCILow, ...
  S.jointGainRatioCIHigh-S.medianJointGainRatio,'ko-','MarkerFaceColor','k','LineWidth',1.4);
formatOffsetAxis(ax,offsets); ylabel(ax,'g_{probe} / g_{pref}');
legend(ax,{'Stored gain-only','Joint bias+gain'},'Location','best');
title(ax,'Gain before and after bias'); grid(ax,'on');

ax = nexttile(layout,3); hold(ax,'on');
errorbar(ax,1:nOffsets,S.medianGainRatioChange, ...
  S.medianGainRatioChange-S.gainRatioChangeCILow, ...
  S.gainRatioChangeCIHigh-S.medianGainRatioChange,'ko-','MarkerFaceColor','k','LineWidth',1.4);
yline(ax,0,':'); formatOffsetAxis(ax,offsets);
ylabel(ax,'Joint - stored gain ratio'); title(ax,'Gain change after adding bias'); grid(ax,'on');

ax = nexttile(layout,4); hold(ax,'on');
errorbar(ax,1:nOffsets,S.medianBiasThresholdRatio, ...
  S.medianBiasThresholdRatio-S.biasThresholdRatioCILow, ...
  S.biasThresholdRatioCIHigh-S.medianBiasThresholdRatio,'ko-','MarkerFaceColor','k','LineWidth',1.4);
yline(ax,1,':'); formatOffsetAxis(ax,offsets);
ylabel(ax,'(parent c_{75} - b_{probe}) / parent c_{75}');
title(ax,'Threshold shift represented by bias'); grid(ax,'on');

ax = nexttile(layout,5); hold(ax,'on');
colors = lines(max(nOffsets,7));
for iOffset = 1:nOffsets
  m = T.probeDirDeg == offsets(iOffset);
  scatter(ax,T.bProbeJoint(m),T.gProbeJointOverGPref(m),25,colors(iOffset,:), ...
    'filled','MarkerFaceAlpha',0.55);
end
xline(ax,0,':'); yline(ax,0,':');
xlabel(ax,'b_{probe} (% coherence)'); ylabel(ax,'Joint g_{probe} / g_{pref}');
title(ax,sprintf('Across-session r = %.3f',population.corrBiasGainAcrossSessions));
grid(ax,'on'); box(ax,'off');

ax = nexttile(layout,6); axis(ax,'off');
textLines = {
  sprintf('ProbeSessions: %d',population.nProbeSessions)
  sprintf('Parent sessions: %d',population.nParentSessions)
  sprintf('Bias limit: +/- %.1f%%',opts.BiasLimit)
  sprintf('Gain limit: +/- %.1f',opts.GainLimit)
  sprintf('Total DeltaNLL bias after gain: %.3f',population.totalDeltaNLLBiasAfterGain)
  sprintf('Aggregate bias-after-gain p: %.3g',population.pAggregateBiasAfterGain)
  sprintf('Total DeltaNLL gain after bias: %.3f',population.totalDeltaNLLGainAfterBias)
  sprintf('Aggregate gain-after-bias p: %.3g',population.pAggregateGainAfterBias)
  sprintf('Bias after gain p<.05: %d/%d',population.nBiasAfterGainP005,population.nProbeSessions)
  sprintf('Gain after bias p<.05: %d/%d',population.nGainAfterBiasP005,population.nProbeSessions)
  sprintf('Median local corr(b,g): %.3f',population.medianLocalCorrBiasGain)
  sprintf('Near bias bound: %d',population.nBiasNearBound)
  sprintf('Near gain bound: %d',population.nGainNearBound)
  sprintf('Max gain-only NLL discrepancy: %.3g',population.maxAbsStoredGainNLLDifference)
  };
text(ax,0.02,0.98,strjoin(textLines,newline),'Units','normalized', ...
  'VerticalAlignment','top','FontName','Menlo','FontSize',9.5,'Interpreter','none');

exportgraphics(fig,pdfPath,'ContentType','vector');
close(fig);
end

% -------------------------------------------------------------------------
function formatOffsetAxis(ax, offsets)
xlim(ax,[0.5 numel(offsets)+0.5]); xticks(ax,1:numel(offsets));
labels = arrayfun(@(x) sprintf('%g%c',x,char(176)),offsets,'UniformOutput',false);
xticklabels(ax,labels); xlabel(ax,'Probe offset'); box(ax,'off');
end
