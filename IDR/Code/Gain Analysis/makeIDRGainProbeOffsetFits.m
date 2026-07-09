function fitData = makeIDRGainProbeOffsetFits(varargin)
% makeIDRGainProbeOffsetFits
% First-pass probe-noise gain fits for IDR Inc/ChangeSide probe sessions.
%
% This function starts the probe stage of the Weibull/gain framework.
% It uses the stabilized preferred-noise fit as calibration and then fits a
% separate probe gain for each probe-session file. 
%
% Required upstream products:
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainAcrossSessionPrefFit_<Animal>.mat
%   Data/AcrossOffsetSummaries/BetaWeights_<Animal>.mat
%   Data/FullSessions/BetaAnalysis/<parentFileName>.mat
%   Data/Probe*/ProbeSessions/*.mat
%
% Model for each probe-session file, restricted to INC trials:
%
%   deltaC_eff = max(deltaC + gPref_parent*nPref + gProbe*nProbe, 0)
%   P(correct) = Weibull(deltaC_eff; alpha_parent, betaShared, lapseShared)
%
% alpha_parent and gPref_parent come from the parent full-session row of the
% stabilized preferred fit. betaShared/lapseShared are shared preferred-fit
% parameters. The only fitted parameter for each probe-session file is
% gProbe, a dimensionless coherence-equivalent gain.
%
% This version is intentionally conservative. It estimates per-file probe
% gains and writes a table suitable for inspection before any pooling or
% hierarchical offset-level model is introduced.
%
% Outputs:
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainProbeOffsetFits_<Animal>.mat
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainProbeOffsetFits_<Animal>_probeTable.csv
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainProbeOffsetFits_<Animal>_offsetSummary.csv

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'Bin180Into179', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'GainLimit', 5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'MinIncTrials', 20, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'ProbeDirDeg', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x)));
addParameter(p, 'Replace', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opts = p.Results;

dataRoot = fullfile(domainFolder(mfilename('fullpath')), 'Data');
outputFolder = validFolder(fullfile(dataRoot, 'AcrossOffsetSummaries', 'GainAnalysis'));
animalTag = char(string(opts.Animal));

matPath = fullfile(outputFolder, sprintf('IDRGainProbeOffsetFits_%s.mat', animalTag));
% probeCsvPath = fullfile(outputFolder, sprintf('IDRGainProbeOffsetFits_%s_probeTable.csv', animalTag));
% offsetCsvPath = fullfile(outputFolder, sprintf('IDRGainProbeOffsetFits_%s_offsetSummary.csv', animalTag));

if isfile(matPath) && ~opts.Replace
  if opts.Verbose
    fprintf('IDR probe gain fits already exist; loading %s\n', matPath);
  end
  S = load(matPath, 'fitData');
  fitData = S.fitData;
  return;
end

prefFitPath = fullfile(outputFolder, sprintf('IDRGainAcrossSessionPrefFit_%s.mat', animalTag));
PF = load(prefFitPath, 'fitData');
prefFit = PF.fitData;

weightPath = fullfile(dataRoot, 'AcrossOffsetSummaries', sprintf('BetaWeights_%s.mat', animalTag));
W = load(weightPath, 'weightData');
weightData = W.weightData;

probeDirs = dir(fullfile(dataRoot, 'Probe*'));
probeDirs = probeDirs([probeDirs.isdir]);
[~, order] = sort({probeDirs.name});
probeDirs = probeDirs(order);

allRows = repmat(emptyProbeRow(), 0, 1);
allTrialTables = cell(0,1);

for iProbeDir = 1:numel(probeDirs)
  probeTag = probeDirs(iProbeDir).name;
  sessionFolder = fullfile(probeDirs(iProbeDir).folder, probeTag, 'ProbeSessions');
  if ~isfolder(sessionFolder)
    continue;
  end

  try
    selectorArgs = {'Animal', opts.Animal};
    if ~isempty(opts.ProbeDirDeg)
      selectorArgs = [selectorArgs, {'ProbeDirDeg', opts.ProbeDirDeg}]; %#ok<AGROW>
    end
    [selectedFiles, fileInfo] = selectAnalysisFiles(sessionFolder, selectorArgs{:});
  catch ME
    warning('makeIDRGainProbeOffsetFits:SelectionFailed', ...
      'Skipping %s because selection failed: %s', sessionFolder, ME.message);
    continue;
  end

  for iFile = 1:numel(selectedFiles)
    sourcePath = selectedFiles{iFile};
    row = emptyProbeRow();
    row.animal = string(fileInfo.animal{iFile});
    row.probeTag = string(probeTag);
    row.probeSessionName = string(stripExtension(fileInfo.fileName{iFile}));
    row.fileName = string(fileInfo.fileName{iFile});
    row.filePath = string(sourcePath);
    row.probeDirDeg = fileInfo.probeDirDeg(iFile);
    row.parentFileName = string(fileInfo.parentFileName{iFile});

    if opts.Verbose
      fprintf('Fitting probe gain: %s\n', fileInfo.fileName{iFile});
    end

    try
      D = loadProbeGainData(sourcePath, weightData, prefFit);
      row.parentSessionName = string(D.parentSessionName);
      row.parentWeightRow = D.parentWeightRow;
      row.prefFitSessionIndex = D.prefFitSessionIndex;
      row.alphaParent = D.alphaParent;
      row.gPrefParent = D.gPrefParent;
      row.betaShared = D.betaShared;
      row.lapseShared = D.lapseShared;
      row.nTrials = D.nTrialsAll;
      row.nIncTrials = numel(D.correct);
      row.nCoherenceLevels = numel(unique(D.deltaC));
      row.minCoherence = min(D.deltaC);
      row.maxCoherence = max(D.deltaC);
      row.fractionCorrect = mean(D.correct);
      row.prefNoiseMean = mean(D.nPref);
      row.prefNoiseSD = std(D.nPref);
      row.probeNoiseMean = mean(D.nProbe);
      row.probeNoiseSD = std(D.nProbe);
      row.probeNoiseMin = min(D.nProbe);
      row.probeNoiseMax = max(D.nProbe);
      row.prefCohNoisePC = D.normInfo.prefCohNoisePC;
      row.probeCohNoisePC = D.normInfo.probeCohNoisePC;
      row.nYokedProbeStreams = D.normInfo.nYokedProbeStreams;
      row.combinedProbeCohNoisePC = D.normInfo.combinedProbeCohNoisePC;
      row.kernelProbeNormFactor = D.normInfo.probeNormFactor;
      row.gainNormalizationMethod = string(D.normInfo.gainMethod);
      row.fitIncluded = row.nIncTrials >= opts.MinIncTrials && row.nCoherenceLevels >= 1;

      if row.probeDirDeg == 180 && opts.Bin180Into179 
        row.probeDirDeg = 179;
      end

      if ~row.fitIncluded
        row.message = sprintf('Too few INC trials (%d) for probe gain fit.', row.nIncTrials);
      else
        fit = fitProbeGain(D.deltaC, D.nPref, D.nProbe, D.correct, ...
          D.alphaParent, D.gPrefParent, D.betaShared, D.lapseShared, opts.GainLimit);
        row.gProbe = fit.gProbe;
        row.nllProbe = fit.nll;
        row.exitflagProbe = fit.exitflag;
        row.messageProbe = string(fit.message);

        pNoProbe = probePcorrect(D.deltaC, D.nPref, D.nProbe, ...
          D.alphaParent, D.gPrefParent, 0, D.betaShared, D.lapseShared);
        row.nllNoProbe = bernoulliNLL(pNoProbe, D.correct);
        row.deltaNLLProbe = row.nllNoProbe - row.nllProbe;
        row.deltaDevianceProbe = 2 * row.deltaNLLProbe;
        row.pApproxProbe = 1 - chi2cdfLocal(row.deltaDevianceProbe, 1);
        row.nllImprovementPerTrial = row.deltaNLLProbe / row.nIncTrials;
        row.gProbeOverGPrefParent = row.gProbe / D.gPrefParent;
        row.gProbeNearBound = abs(row.gProbe) > 0.95 * opts.GainLimit;
        row.message = "ok";

        pHat = probePcorrect(D.deltaC, D.nPref, D.nProbe, ...
          D.alphaParent, D.gPrefParent, row.gProbe, D.betaShared, D.lapseShared);
        Ttrial = table();
        Ttrial.probeSessionIndex = repmat(numel(allRows)+1, numel(D.correct), 1);
        Ttrial.probeSessionName = repmat(row.probeSessionName, numel(D.correct), 1);
        Ttrial.parentSessionName = repmat(row.parentSessionName, numel(D.correct), 1);
        Ttrial.probeDirDeg = repmat(row.probeDirDeg, numel(D.correct), 1);
        Ttrial.deltaC = D.deltaC;
        Ttrial.nPref = D.nPref;
        Ttrial.nProbe = D.nProbe;
        Ttrial.correct = logical(D.correct);
        Ttrial.pHat = pHat;
        Ttrial.deltaCEff = max(D.deltaC + D.gPrefParent .* D.nPref + row.gProbe .* D.nProbe, 0);
        allTrialTables{end+1,1} = Ttrial; %#ok<AGROW>
      end
    catch ME
      row.fitIncluded = false;
      row.message = string(ME.message);
      if opts.Verbose
        fprintf('  failed: %s\n', ME.message);
      end
    end


    allRows(end+1,1) = row; %#ok<AGROW>
  end
end

if isempty(allRows)
  error('makeIDRGainProbeOffsetFits:NoProbeFiles', ...
    'No probe-session files were found or selected.');
end

probeTable = struct2table(allRows);
probeTable.probeSessionIndex = (1:height(probeTable))';
probeTable = movevars(probeTable, 'probeSessionIndex', 'Before', 1);

if isempty(allTrialTables)
  trialTable = table();
else
  trialTable = vertcat(allTrialTables{:});
end

offsetSummary = summarizeByProbeDirection(probeTable);

population = summarizePopulation(probeTable);

metadata = struct();
metadata.version = 4;
metadata.analysisName = 'IDRGainProbeOffsetFits';
metadata.method = ['First-pass probe gain fits; beta/lapse and parent alpha/gPref fixed from ' ...
  'the stabilized preferred-noise fit; one gProbe per probe-session file. ' ...
  'The nominal preferred-direction step coherence is recovered by matching probe trialIdx ' ...
  'to the parent BetaAnalysis sessionNoise.trialIdx after filtering to INC rows. ' ...
  'probeNoiseByPatch is treated as effective yoked probe noise in percent coherence. ' ...
  'Thus gProbe/gPrefParent is directly comparable to the variance-normalized kernel scale; ' ...
  'the kernel squared-amplitude normalization factor is stored only for provenance, not applied to the predictor.'];
metadata.model = ['deltaC_eff = max(deltaC + gPref_parent*nPref + gProbe*nProbe, 0); ' ...
  'Pcorrect = Weibull(deltaC_eff; alpha_parent, betaShared, lapseShared)'];
metadata.prefFitPath = prefFitPath;
metadata.weightPath = weightPath;
metadata.options = opts;
metadata.createdBy = mfilename;
metadata.createdDate = datetime('now');

fitData = struct();
fitData.version = 4;
fitData.probeTable = probeTable;
fitData.offsetSummary = offsetSummary;
fitData.trialTable = trialTable;
fitData.population = population;
fitData.metadata = metadata;

save(matPath, 'fitData', '-v7.3');

if opts.Verbose
  fprintf('\nProbe-gain fit summary:\n');
  fprintf('  probe sessions fit:      %d/%d\n', sum(probeTable.fitIncluded), height(probeTable));
  fprintf('  median gProbe:           %.4g\n', population.medianGProbe);
  fprintf('  median gProbe/gPref:     %.4g\n', population.medianGProbeOverGPrefParent);
  fprintf('  median delta NLL:        %.4g\n', population.medianDeltaNLLProbe);
  fprintf('  delta NLL > 0:           %d/%d\n', population.nDeltaNLLPositive, population.nFitProbeSessions);
  fprintf('  near gain bound:         %d/%d\n', population.nNearGainBound, population.nFitProbeSessions);
  fprintf('Saved MAT: %s\n', matPath);
end
end

% -------------------------------------------------------------------------
function row = emptyProbeRow()
row = struct();
row.animal = "";
row.probeTag = "";
row.probeSessionIndex = NaN;
row.probeSessionName = "";
row.fileName = "";
row.filePath = "";
row.parentFileName = "";
row.parentSessionName = "";
row.probeDirDeg = NaN;
row.parentWeightRow = NaN;
row.prefFitSessionIndex = NaN;
row.alphaParent = NaN;
row.gPrefParent = NaN;
row.betaShared = NaN;
row.lapseShared = NaN;
row.nTrials = NaN;
row.nIncTrials = NaN;
row.nCoherenceLevels = NaN;
row.minCoherence = NaN;
row.maxCoherence = NaN;
row.fractionCorrect = NaN;
row.prefNoiseMean = NaN;
row.prefNoiseSD = NaN;
row.probeNoiseMean = NaN;
row.probeNoiseSD = NaN;
row.probeNoiseMin = NaN;
row.probeNoiseMax = NaN;
row.prefCohNoisePC = NaN;
row.probeCohNoisePC = NaN;
row.nYokedProbeStreams = NaN;
row.combinedProbeCohNoisePC = NaN;
row.kernelProbeNormFactor = NaN;
row.gainNormalizationMethod = "";
row.fitIncluded = false;
row.gProbe = NaN;
row.nllProbe = NaN;
row.exitflagProbe = NaN;
row.messageProbe = "";
row.nllNoProbe = NaN;
row.deltaNLLProbe = NaN;
row.deltaDevianceProbe = NaN;
row.pApproxProbe = NaN;
row.nllImprovementPerTrial = NaN;
row.gProbeOverGPrefParent = NaN;
row.gProbeNearBound = false;
row.message = "";
end

% -------------------------------------------------------------------------
function D = loadProbeGainData(sourcePath, W, prefFit)
S = load(sourcePath);
required = {'sessionHeader','sessionProbeHeader','prefNoiseByPatch', ...
  'probeNoiseByPatch','trialOutcomesAll','changeSidesAll','changeIndicesAll'};
missing = required(~isfield(S, required));
if ~isempty(missing)
  error('Probe-session file is missing required fields: %s', strjoin(missing, ', '));
end

normInfo = normalizationInfoFromHeadersLocal(S.sessionHeader, S.sessionProbeHeader);

[parentWeightRow, parentFileName] = matchParentWeightRowLocal(S.sessionHeader, S.sessionProbeHeader, W);
weights = double(W.leaveOneOutWeights(parentWeightRow, :));
stepTMS = double(W.stepTMS(:)');

[prefStep, prefStepTMS] = extractStepFramesLocal(S.prefNoiseByPatch, S.sessionHeader);
[probeStep, probeStepTMS] = extractStepFramesLocal(S.probeNoiseByPatch, S.sessionHeader);
if numel(weights) ~= size(prefStep,2) || numel(weights) ~= size(probeStep,2)
  error('Weight count (%d) does not match step-frame count (%d).', numel(weights), size(prefStep,2));
end
if max(abs(prefStepTMS - stepTMS)) > 1e-6 || max(abs(probeStepTMS - stepTMS)) > 1e-6
  error('Probe-session step-frame times do not match BetaWeights.stepTMS.');
end
if abs(sum(weights) - 1) > 1e-12
  error('Leave-one-session-out weights do not sum to one.');
end

prefByPatch = squeeze(sum(prefStep .* reshape(weights, 1, [], 1), 2));
probeByPatch = squeeze(sum(probeStep .* reshape(weights, 1, [], 1), 2));
if isvector(prefByPatch), prefByPatch = reshape(prefByPatch, 2, []); end
if isvector(probeByPatch), probeByPatch = reshape(probeByPatch, 2, []); end

changeSides = double(S.changeSidesAll(:));
changeIndex = double(S.changeIndicesAll(:));
correctAll = double(S.trialOutcomesAll(:) == 0);
changePatch = normalizePatchIndexLocal(changeSides);
nTrials = numel(correctAll);

if numel(changeIndex) ~= nTrials || numel(changePatch) ~= nTrials
  error('Trial-level vectors have inconsistent lengths.');
end

xPref = nan(nTrials,1);
xProbe = nan(nTrials,1);
for i = 1:nTrials
  xPref(i) = prefByPatch(changePatch(i), i);
  xProbe(i) = probeByPatch(changePatch(i), i);
end

% IDR convention in existing probe regression products: 1=DEC, 2=INC.
% ProbeSessions are intentionally general and may contain both DEC and INC
% validated trials.  This gain analysis is currently restricted to INC, to
% match the parent BetaAnalysis/sessionNoise product used by the preferred
% gain fit.  Therefore filter all trial-aligned probe variables before
% matching trialIdx into the parent BetaAnalysis file.
useInc = changeIndex == 2;
if ~any(useInc)
  error('Probe-session contains no INC trials.');
end

if ~isfield(S, 'trialIdx')
  error('Probe-session file is missing trialIdx.');
end
probeTrialIdxAll = numericVectorOrEmpty(S.trialIdx);
if numel(probeTrialIdxAll) ~= nTrials
  error(['Probe-session trialIdx length (%d) does not match the number of saved ' ...
    'validated probe rows (%d). Rebuild ProbeSessions with aligned trialIdx.'], ...
    numel(probeTrialIdxAll), nTrials);
end
probeTrialIdxInc = probeTrialIdxAll(useInc);
deltaCInc = extractDeltaCFromParent(sourcePath, S, parentFileName, probeTrialIdxInc);

[prefIndex, parentSessionName] = matchParentPrefFit(parentFileName, prefFit.sessionTable);

D = struct();
D.parentWeightRow = parentWeightRow;
D.parentFileName = parentFileName;
D.parentSessionName = parentSessionName;
D.prefFitSessionIndex = prefIndex;
D.alphaParent = prefFit.sessionTable.alphaSession(prefIndex);
D.gPrefParent = prefFit.sessionTable.gPrefSession(prefIndex);
D.betaShared = prefFit.params.betaShared;
D.lapseShared = prefFit.params.lapseShared;
D.nTrialsAll = nTrials;
D.deltaC = double(deltaCInc);
D.nPref = double(xPref(useInc));
% probeNoiseByPatch is already the effective yoked probe perturbation, in
% percent coherence units.  Do not apply the kernel squared-amplitude
% normalization factor to this likelihood predictor; that factor converts
% reverse-correlation kernel amplitudes, which scale with noise variance.
D.nProbe = double(xProbe(useInc));
D.normInfo = normInfo;
D.correct = double(correctAll(useInc));

finite = isfinite(D.deltaC) & D.deltaC >= 0 & isfinite(D.nPref) & isfinite(D.nProbe) & isfinite(D.correct);
D.deltaC = D.deltaC(finite);
D.nPref = D.nPref(finite);
D.nProbe = D.nProbe(finite);
D.correct = D.correct(finite);
end

% -------------------------------------------------------------------------
function deltaC = extractDeltaCFromParent(sourcePath, S, parentFileName, probeTrialIdx)
% Recover the nominal preferred-direction INC step coherence for each
% included probe trial by matching its parent-session trialIdx back to the
% parent BetaAnalysis sessionNoise.trialIdx.
%
% There is no probe-direction step or probe-direction signal coherence.  The
% Weibull is driven by the preferred-direction step, with preferred and
% probe noise added as trialwise effective-coherence perturbations.
probeTrialIdx = numericVectorOrEmpty(probeTrialIdx);
if isempty(probeTrialIdx)
  error('No INC probe trialIdx values were supplied for parent matching.');
end

parentPath = parentBetaAnalysisPath(sourcePath, S, parentFileName);
if ~isfile(parentPath)
  error('Could not find parent BetaAnalysis file: %s', parentPath);
end
P = load(parentPath, 'sessionNoise');
if ~isfield(P, 'sessionNoise')
  error('Parent BetaAnalysis file does not contain sessionNoise: %s', parentPath);
end
SN = P.sessionNoise;
if ~isfield(SN, 'signalCohPC')
  error('Parent sessionNoise is missing signalCohPC: %s', parentPath);
end
parentDeltaC = numericVectorOrEmpty(SN.signalCohPC);
if isempty(parentDeltaC)
  error('Parent sessionNoise.signalCohPC is not numeric or is empty: %s', parentPath);
end

if isfield(SN, 'trialIdx')
  parentTrialIdx = numericVectorOrEmpty(SN.trialIdx);
elseif isfield(SN, 'trialIndex')
  parentTrialIdx = numericVectorOrEmpty(SN.trialIndex);
else
  error('Parent sessionNoise is missing trialIdx; cannot map probe trials to parent signalCohPC.');
end
if numel(parentTrialIdx) ~= numel(parentDeltaC)
  error('Parent sessionNoise.trialIdx and signalCohPC have inconsistent lengths.');
end

[tf, loc] = ismember(probeTrialIdx(:), parentTrialIdx(:));
if ~all(tf)
  nMissing = sum(~tf);
  exampleMissing = probeTrialIdx(find(~tf, 1, 'first'));
  error(['Could not match %d INC probe trialIdx values to parent sessionNoise.trialIdx. ' ...
    'Example missing trialIdx: %g. This usually means the probe-session rows are not ' ...
    'aligned, or the downstream INC filter does not match the parent BetaAnalysis selection.'], ...
    nMissing, exampleMissing);
end

deltaC = double(parentDeltaC(loc));
end

% -------------------------------------------------------------------------
function parentPath = parentBetaAnalysisPath(sourcePath, S, parentFileName)
% Parent path:
%   <domain root>/Data/FullSessions/BetaAnalysis/<parentFileName>.mat
% Use domainFolder rather than parsing sourcePath so that paths containing
% repeated "Data" components do not create duplicated path fragments.

if nargin < 3 || isempty(parentFileName)
  if isfield(S, 'sessionProbeHeader') && isfield(S.sessionProbeHeader, 'parentFileName')
    parentFileName = S.sessionProbeHeader.parentFileName;
  elseif isfield(S, 'sessionHeader') && isfield(S.sessionHeader, 'fileName')
    parentFileName = S.sessionHeader.fileName;
  else
    error('Cannot determine parentFileName for probe source: %s', sourcePath);
  end
end

parentStem = char(string(parentFileName));
[~, parentStem, parentExt] = fileparts(parentStem);
if isempty(parentExt)
  parentExt = '.mat';
end
parentBase = [parentStem parentExt];

parentPath = fullfile(domainFolder(mfilename('fullpath')), ...
  'Data', 'FullSessions', 'BetaAnalysis', parentBase);
end

% -------------------------------------------------------------------------
function v = numericVectorOrEmpty(x)
while isstruct(x) && isfield(x, 'data')
  x = x.data;
end
if isnumeric(x) || islogical(x)
  v = double(x(:));
else
  v = [];
end
end

% -------------------------------------------------------------------------
function fit = fitProbeGain(deltaC, nPref, nProbe, correct, alphaParent, gPrefParent, betaShared, lapseShared, gainLimit)
theta0 = gainToThetaLocal(0, gainLimit);
objective = @(theta) probeGainNLL(theta, deltaC, nPref, nProbe, correct, ...
  alphaParent, gPrefParent, betaShared, lapseShared, gainLimit);
options = optimset('Display', 'off', 'MaxIter', 1000, 'MaxFunEvals', 3000, 'TolX', 1e-8, 'TolFun', 1e-8);
[thetaHat, nll, exitflag, output] = fminsearch(objective, theta0, options);
fit.gProbe = gainLimit * tanh(thetaHat(1));
fit.nll = nll;
fit.exitflag = exitflag;
fit.message = output.message;
end

% -------------------------------------------------------------------------
function nll = probeGainNLL(theta, deltaC, nPref, nProbe, correct, alphaParent, gPrefParent, betaShared, lapseShared, gainLimit)
gProbe = gainLimit * tanh(theta(1));
p = probePcorrect(deltaC, nPref, nProbe, alphaParent, gPrefParent, gProbe, betaShared, lapseShared);
nll = bernoulliNLL(p, correct);
end

% -------------------------------------------------------------------------
function p = probePcorrect(deltaC, nPref, nProbe, alphaParent, gPrefParent, gProbe, betaShared, lapseShared)
deltaCEff = max(deltaC + gPrefParent .* nPref + gProbe .* nProbe, 0);
p = weibullPcorrect(deltaCEff, alphaParent, betaShared, lapseShared);
end

% -------------------------------------------------------------------------
function p = weibullPcorrect(deltaC, alpha, beta, lapse)
p = 0.5 + (0.5 - lapse) .* (1 - exp(-((deltaC ./ alpha) .^ beta)));
end

% -------------------------------------------------------------------------
function nll = bernoulliNLL(p, correct)
epsP = 1e-12;
p = min(max(p, epsP), 1 - epsP);
nll = -sum(correct .* log(p) + (1 - correct) .* log(1 - p));
if ~isfinite(nll)
  nll = realmax;
end
end

% -------------------------------------------------------------------------
function theta = gainToThetaLocal(g, gainLimit)
z = max(min(g ./ gainLimit, 0.95), -0.95);
theta = atanh(z);
end

% -------------------------------------------------------------------------
function out = summarizeByProbeDirection(T)
fitMask = T.fitIncluded & isfinite(T.probeDirDeg) & isfinite(T.gProbe);
if ~any(fitMask)
  out = table();
  return;
end
probeDirs = unique(T.probeDirDeg(fitMask));
probeDirs = sort(probeDirs(:));
rows = cell(numel(probeDirs), 1);
for i = 1:numel(probeDirs)
  d = probeDirs(i);
  m = fitMask & T.probeDirDeg == d;
  g = T.gProbe(m);
  r = T.gProbeOverGPrefParent(m);
  dnll = T.deltaNLLProbe(m);
  rows{i} = table(d, sum(m), sum(T.nIncTrials(m)), mean(g), median(g), quantileLocal(g,0.25), quantileLocal(g,0.75), ...
    mean(r), median(r), median(dnll), sum(dnll > 0), sum(T.pApproxProbe(m) < 0.05), ...
    'VariableNames', {'probeDirDeg','nProbeSessions','nIncTrials','meanGProbe','medianGProbe', ...
    'q25GProbe','q75GProbe','meanGProbeOverGPref','medianGProbeOverGPref','medianDeltaNLLProbe', ...
    'nDeltaNLLPositive','nApproxP005'});
end
out = vertcat(rows{:});
end

% -------------------------------------------------------------------------
function population = summarizePopulation(T)
fitMask = T.fitIncluded & isfinite(T.gProbe);
g = T.gProbe(fitMask);
r = T.gProbeOverGPrefParent(fitMask);
dnll = T.deltaNLLProbe(fitMask);
p = T.pApproxProbe(fitMask);
population = struct();
population.nProbeSessions = height(T);
population.nFitProbeSessions = sum(fitMask);
population.nIncTrials = sum(T.nIncTrials(fitMask));
population.meanGProbe = meanFinite(g);
population.medianGProbe = medianFinite(g);
population.q25GProbe = quantileLocal(g, 0.25);
population.q75GProbe = quantileLocal(g, 0.75);
population.iqrGProbe = population.q75GProbe - population.q25GProbe;
population.meanGProbeOverGPrefParent = meanFinite(r);
population.medianGProbeOverGPrefParent = medianFinite(r);
population.medianDeltaNLLProbe = medianFinite(dnll);
population.nDeltaNLLPositive = sum(dnll > 0);
population.nApproxP005 = sum(p < 0.05);
population.nApproxP001 = sum(p < 0.001);
population.nNearGainBound = sum(T.gProbeNearBound(fitMask));
population.totalDeltaNLLProbe = sum(dnll(isfinite(dnll)));
end

% -------------------------------------------------------------------------
function m = meanFinite(x)
x = x(isfinite(x));
if isempty(x), m = NaN; else, m = mean(x); end
end

% -------------------------------------------------------------------------
function m = medianFinite(x)
x = x(isfinite(x));
if isempty(x), m = NaN; else, m = median(x); end
end

% -------------------------------------------------------------------------
function q = quantileLocal(x, p)
x = sort(x(isfinite(x)));
if isempty(x)
  q = NaN;
  return;
end
if numel(x) == 1
  q = x;
  return;
end
pos = 1 + (numel(x)-1) * p;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
  q = x(lo);
else
  q = x(lo) + (pos-lo) * (x(hi)-x(lo));
end
end

% -------------------------------------------------------------------------
function p = chi2cdfLocal(x, df)
if df ~= 1
  p = gammainc(x/2, df/2);
else
  p = erf(sqrt(max(x,0) ./ 2));
end
end

% -------------------------------------------------------------------------
function [stepNoise, stepTMS] = extractStepFramesLocal(noiseByPatch, sessionHeader)
if ndims(noiseByPatch) ~= 3 || size(noiseByPatch,1) ~= 2
  error('Noise matrix must be 2 x frames x trials.');
end
frameRateHz = headerScalarLocal(sessionHeader, 'frameRateHz');
preStepMS = headerScalarLocal(sessionHeader, 'preStepMS');
stepMS = headerScalarLocal(sessionHeader, 'stepMS');
msPerFrame = 1000 / frameRateHz;
nFrames = size(noiseByPatch, 2);
tMS = (0:nFrames-1) * msPerFrame;
mask = tMS >= preStepMS & tMS < preStepMS + stepMS;
stepNoise = double(noiseByPatch(:, mask, :));
stepTMS = tMS(mask);
end

% -------------------------------------------------------------------------
function idx = normalizePatchIndexLocal(raw)
raw = raw(:);
vals = unique(raw(isfinite(raw)));
if all(ismember(vals, [0 1]))
  idx = raw + 1;
elseif all(ismember(vals, [1 2]))
  idx = raw;
else
  error('changeSidesAll must contain 0/1 or 1/2 patch indices.');
end
idx = double(idx(:));
end

% -------------------------------------------------------------------------
function [row, parentName] = matchParentWeightRowLocal(sessionHeader, probeHeader, W)
candidates = strings(0,1);
if isfield(probeHeader, 'parentFileName')
  candidates(end+1) = string(probeHeader.parentFileName);
end
if isfield(sessionHeader, 'fileName')
  candidates(end+1) = string(sessionHeader.fileName);
end
if isfield(probeHeader, 'sessionID')
  candidates(end+1) = string(probeHeader.sessionID);
end
weightNames = string(W.sessionFileNames(:));
weightBases = strings(size(weightNames));
for i = 1:numel(weightNames)
  [~, b, e] = fileparts(weightNames(i));
  weightBases(i) = b + e;
end
matches = false(size(weightNames));
for c = 1:numel(candidates)
  [~, b, e] = fileparts(candidates(c));
  fullCandidate = b + e;
  if strlength(e) == 0
    fullCandidate = b + ".mat";
  end
  matches = matches | strcmpi(weightBases, fullCandidate);
end
row = find(matches);
if numel(row) ~= 1
  error('Could not uniquely match probe session to BetaWeights (%d matches).', numel(row));
end
parentName = char(weightNames(row));
end

% -------------------------------------------------------------------------
function [idx, parentSessionName] = matchParentPrefFit(parentFileName, sessionTable)
[~, pb, pe] = fileparts(parentFileName);
parentBaseExt = string(pb) + string(pe);
if strlength(string(pe)) == 0
  parentBaseExt = string(pb) + ".mat";
end
fitNames = string(sessionTable.fileName(:));
fitBases = strings(size(fitNames));
for i = 1:numel(fitNames)
  [~, b, e] = fileparts(fitNames(i));
  fitBases(i) = b + e;
end
matches = strcmpi(fitBases, parentBaseExt);
idx = find(matches);
if numel(idx) ~= 1
  error('Could not uniquely match parent %s to preferred-fit sessionTable (%d matches).', parentBaseExt, numel(idx));
end
parentSessionName = char(sessionTable.sessionName(idx));
end


% -------------------------------------------------------------------------
function normInfo = normalizationInfoFromHeadersLocal(sessionHeader, sessionProbeHeader)
prefCohNoisePC  = headerScalarLocal(sessionHeader, 'prefCohNoisePC');
probeCohNoisePC = headerScalarLocal(sessionProbeHeader, 'probeCohNoisePC');

if ~isfinite(prefCohNoisePC) || ~isfinite(probeCohNoisePC) || probeCohNoisePC <= 0
  error('makeIDRGainProbeOffsetFits:BadNoiseAmplitude', ...
    'Invalid pref/probe coherence noise amplitudes: pref=%g, probe=%g.', prefCohNoisePC, probeCohNoisePC);
end

nYokedProbeStreams = probeStreamCountFromSessionProbeHeaderLocal(sessionProbeHeader);
combinedProbeCohNoisePC = nYokedProbeStreams * probeCohNoisePC;

normInfo = struct();
normInfo.prefCohNoisePC = prefCohNoisePC;
normInfo.probeCohNoisePC = probeCohNoisePC;
normInfo.nYokedProbeStreams = nYokedProbeStreams;
normInfo.combinedProbeCohNoisePC = combinedProbeCohNoisePC;
normInfo.probeNormFactor = (prefCohNoisePC / combinedProbeCohNoisePC)^2;
normInfo.kernelMethod = ...
  'probe kernels multiplied by (prefCohNoisePC/(nYokedProbeStreams*probeCohNoisePC))^2 before normalized kernel scales';
normInfo.gainMethod = ...
  'gain fit uses effective yoked probe noise in percent coherence; gProbe/gPrefParent is directly comparable to normalized kernel scale';
end

% -------------------------------------------------------------------------
function n = probeStreamCountFromSessionProbeHeaderLocal(sessionProbeHeader)
if ~isfield(sessionProbeHeader, 'probeDirDeg')
  error('makeIDRGainProbeOffsetFits:MissingProbeDir', 'sessionProbeHeader.probeDirDeg is required.');
end
probeDirDeg = abs(double(sessionProbeHeader.probeDirDeg));

if probeDirDeg > 0 && probeDirDeg < 180
  n = 2;
elseif abs(probeDirDeg - 180) < 1e-9
  n = 1;
else
  error('makeIDRGainProbeOffsetFits:UnsupportedProbeDir', ...
    'Unsupported probeDirDeg for probe normalization: %g.', probeDirDeg);
end
end

% -------------------------------------------------------------------------
function value = headerScalarLocal(H, fieldName)
if ~isfield(H, fieldName)
  error('sessionHeader.%s is required.', fieldName);
end
value = H.(fieldName);
while isstruct(value) && isfield(value, 'data')
  value = value.data;
end
if ~(isnumeric(value) || islogical(value)) || ~isscalar(value)
  error('sessionHeader.%s must be a numeric scalar.', fieldName);
end
value = double(value);
end

% -------------------------------------------------------------------------
function stem = stripExtension(fileName)
[~, stem] = fileparts(char(fileName));
end
