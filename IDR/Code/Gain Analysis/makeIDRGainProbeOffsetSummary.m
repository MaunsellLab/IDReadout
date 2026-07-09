function summaryData = makeIDRGainProbeOffsetSummary(varargin)
% makeIDRGainProbeOffsetSummary
% Summarize IDR probe-noise gain fits by probe offset/direction and write
% diagnostic tables plus session-bootstrap confidence intervals for offset medians.
%
% Upstream product:
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainProbeOffsetFits_<Animal>.mat
%
% This function does not refit anything. It reads the per-probe-session
% fits, groups fitted rows by probeDirDeg, writes compact CSV summaries,
% and makes a one-page PDF diagnostic. Bootstrap CIs are computed by resampling fitted probe-session rows within each probe offset.
%
% Usage:
%   S = makeIDRGainProbeOffsetSummary();
%   S = makeIDRGainProbeOffsetSummary('Animal','All','Replace',true);
%
% Outputs:
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainProbeOffsetSummary_<Animal>.mat
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainProbeOffsetSummary_<Animal>_offsetSummary.csv
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainProbeOffsetSummary_<Animal>_population.csv
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainProbeOffsetSummary_<Animal>_failedProbeSessions.csv
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainProbeOffsetSummary_<Animal>_nearBoundProbeSessions.csv
%   Plots/AcrossProbes/GainAnalysis/IDRGainProbeOffsetSummary_<Animal>.pdf

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'Replace', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Bootstrap', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'NBootstrap', 10000, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'BootstrapAlpha', 0.05, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x < 1);
addParameter(p, 'BootstrapSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x)));
parse(p, varargin{:});
opts = p.Results;

root = domainFolder(mfilename('fullpath'));
animalTag = char(string(opts.Animal));
outStem = sprintf('IDRGainProbeOffsetSummary_%s', animalTag);
gainFolder = validFolder(fullfile(root, 'Data', 'AcrossOffsetSummaries', 'GainAnalysis'));
probeFitPath = fullfile(gainFolder, sprintf('IDRGainProbeOffsetFits_%s.mat', animalTag));
plotFolder = validFolder(fullfile(root, 'Plots', 'AcrossProbes', 'GainAnalysis'));
pdfPath = fullfile(plotFolder, [outStem '.pdf']);
matPath = fullfile(gainFolder, [outStem '.mat']);

if isfile(matPath) && ~opts.Replace
  if opts.Verbose
    fprintf('Probe offset summary already exists; loading %s\n', matPath);
  end
  X = load(matPath, 'summaryData');
  summaryData = X.summaryData;
  return;
end

X = load(probeFitPath, 'fitData');
fitData = X.fitData;
probeTable = fitData.probeTable;
fitMask = probeTable.fitIncluded & isfinite(probeTable.probeDirDeg) & isfinite(probeTable.gProbe);
if ~any(fitMask)
  error('No fitted probe sessions were available for offset summary.');
end

offsetSummary = summarizeOffsets(probeTable, fitMask);
bootstrapByOffset = bootstrapOffsets(probeTable, fitMask, opts.NBootstrap, opts.BootstrapAlpha, opts.BootstrapSeed, opts.Bootstrap);
failedProbeSessions = summarizeFailedProbeSessions(probeTable, fitMask);
nearBoundProbeSessions = summarizeNearBoundProbeSessions(probeTable, fitMask);
offsetSummary = mergeBootstrapIntoOffsetSummary(offsetSummary, bootstrapByOffset);
population = summarizePopulation(probeTable, fitMask, offsetSummary, failedProbeSessions, nearBoundProbeSessions);
populationTable = struct2table(population, 'AsArray', true);
populationTall = rowsFromScalarStruct(population);

metadata = struct();
metadata.version = 4;
metadata.analysisName = 'IDRGainProbeOffsetSummary';
metadata.method = ['Groups makeIDRGainProbeOffsetFits probe-session rows by probeDirDeg. ' ...
  'No refitting is performed. gProbeOverGPrefParent uses each probe row''s parent-session preferred gain. ' ...
  'v4 adds a session-level bootstrap over fitted probe-session gains to estimate offset-level ' ...
  'confidence intervals, SEs, and variances.'];
metadata.probeFitPath = probeFitPath;
metadata.options = opts;
metadata.createdBy = mfilename;
metadata.createdDate = datetime('now');
metadata.bootstrapStatistic = 'medianGProbeOverGPref';
metadata.bootstrapUnit = 'probeSession';
metadata.bootstrapVarianceField = 'bootstrapVarGProbeOverGPref';
metadata.bootstrapSEField = 'bootstrapSEGProbeOverGPref';
metadata.readoutFitMeasurementField = 'medianGProbeOverGPref';
metadata.readoutFitOffsetField = 'probeDirDeg';
metadata.anchorPolicy = ...
  '0 deg is an analytical normalization anchor: value=1, variance=0, excluded from fitAcrossOffsetReadout.';

summaryData = struct();
summaryData.version = 4;
summaryData.probeTable = probeTable;
summaryData.offsetSummary = offsetSummary;
summaryData.failedProbeSessions = failedProbeSessions;
summaryData.nearBoundProbeSessions = nearBoundProbeSessions;
summaryData.bootstrapByOffset = bootstrapByOffset;
summaryData.population = population;
summaryData.populationTable = populationTable;
summaryData.populationTall = populationTall;
summaryData.metadata = metadata;

save(matPath, 'summaryData', '-v7.3');
makeSummaryPDF(probeTable, fitMask, offsetSummary, population, pdfPath, animalTag);
makeIDRGainAcrossOffsetReadoutFits('Animal', opts.Animal);

if opts.Verbose
  fprintf('\nProbe gain by offset summary:\n');
  fprintf('  fitted probe sessions:       %d/%d\n', population.nFitProbeSessions, population.nProbeSessions);
  fprintf('  failed/omitted sessions:     %d\n', population.nFailedProbeSessions);
  fprintf('  offsets:                     %d\n', population.nOffsets);
  fprintf('  INC trials:                  %d\n', population.nIncTrials);
  fprintf('  median gProbe:               %.4g\n', population.medianGProbe);
  fprintf('  median gProbe/gPref:         %.4g\n', population.medianGProbeOverGPrefParent);
  fprintf('  median offset median ratio:  %.4g\n', population.medianOffsetMedianGProbeOverGPref);
  fprintf('  median Delta NLL probe:      %.4g\n', population.medianDeltaNLLProbe);
  fprintf('  total Delta NLL probe:       %.4g\n', population.totalDeltaNLLProbe);
  fprintf('  Delta NLL > 0:               %d/%d\n', population.nDeltaNLLPositive, population.nFitProbeSessions);
  fprintf('  approx p < 0.05:             %d/%d\n', population.nApproxP005, population.nFitProbeSessions);
  fprintf('  approx p < 0.001:            %d/%d\n', population.nApproxP001, population.nFitProbeSessions);
  fprintf('  near gain bound:             %d/%d\n', population.nNearGainBound, population.nFitProbeSessions);
  fprintf('Saved MAT: %s\n', matPath);
    fprintf('Saved PDF: %s\n', pdfPath);
end


end

%% -------------------------------------------------------------------------
function out = summarizeOffsets(T, fitMask)
probeDirs = unique(T.probeDirDeg(fitMask));
probeDirs = sort(probeDirs(:));
rows = cell(numel(probeDirs), 1);
for i = 1:numel(probeDirs)
  d = probeDirs(i);
  m = fitMask & T.probeDirDeg == d;
  g = T.gProbe(m);
  r = T.gProbeOverGPrefParent(m);
  dnll = T.deltaNLLProbe(m);
  p = T.pApproxProbe(m);
  npt = T.nIncTrials(m);
  psd = T.probeNoiseSD(m);
  near = T.gProbeNearBound(m);

  rows{i} = table( ...
    d, sum(m), sum(npt), medianFinite(npt), ...
    meanFinite(g), medianFinite(g), quantileLocal(g,0.25), quantileLocal(g,0.75), ...
    quantileLocal(g,0.75)-quantileLocal(g,0.25), ...
    meanFinite(r), medianFinite(r), quantileLocal(r,0.25), quantileLocal(r,0.75), ...
    quantileLocal(r,0.75)-quantileLocal(r,0.25), ...
    meanFinite(dnll), medianFinite(dnll), sum(dnll(isfinite(dnll))), ...
    sum(dnll > 0), sum(p < 0.05), sum(p < 0.001), sum(near), ...
    medianFinite(psd), ...
    'VariableNames', {'probeDirDeg','nProbeSessions','nIncTrials','medianIncTrialsPerSession', ...
    'meanGProbe','medianGProbe','q25GProbe','q75GProbe','iqrGProbe', ...
    'meanGProbeOverGPref','medianGProbeOverGPref','q25GProbeOverGPref','q75GProbeOverGPref', ...
    'iqrGProbeOverGPref','meanDeltaNLLProbe','medianDeltaNLLProbe','totalDeltaNLLProbe', ...
    'nDeltaNLLPositive','nApproxP005','nApproxP001','nNearGainBound','medianProbeNoiseSD'});
end
out = vertcat(rows{:});
out = sortrows(out, 'probeDirDeg');
end

% -------------------------------------------------------------------------
function failed = summarizeFailedProbeSessions(T, fitMask)
failedMask = ~fitMask;
if ~any(failedMask)
  failed = emptyFailedTable();
  return;
end
failed = selectExistingVariables(T(failedMask,:), { ...
  'probeSessionIndex','animal','probeTag','probeDirDeg','probeSessionName','fileName','parentFileName', ...
  'nTrials','nIncTrials','nCoherenceLevels','prefFitSessionIndex','parentWeightRow', ...
  'fitIncluded','exitflagProbe','message','messageProbe','filePath'});
failed = addFailureReason(failed);
failed = sortrowsIfPresent(failed, {'probeDirDeg','probeSessionName'});
end

% -------------------------------------------------------------------------
function near = summarizeNearBoundProbeSessions(T, fitMask)
if ~ismember('gProbeNearBound', T.Properties.VariableNames)
  near = emptyNearBoundTable();
  return;
end
nearMask = fitMask & T.gProbeNearBound;
if ~any(nearMask)
  near = emptyNearBoundTable();
  return;
end
near = selectExistingVariables(T(nearMask,:), { ...
  'probeSessionIndex','animal','probeTag','probeDirDeg','probeSessionName','fileName','parentFileName', ...
  'nIncTrials','gProbe','gPrefParent','gProbeOverGPrefParent','deltaNLLProbe','pApproxProbe', ...
  'probeNoiseSD','prefNoiseSD','nYokedProbeStreams','probeCohNoisePC','combinedProbeCohNoisePC', ...
  'message','messageProbe','filePath'});
near = sortrowsIfPresent(near, {'probeDirDeg','probeSessionName'});
end


%% -------------------------------------------------------------------------
function B = bootstrapOffsets(T, fitMask, nBoot, alpha, seed, doBootstrap)
probeDirs = unique(T.probeDirDeg(fitMask));
probeDirs = sort(probeDirs(:));
rows = cell(numel(probeDirs), 1);

if ~doBootstrap || nBoot <= 0
  for i = 1:numel(probeDirs)
    d = probeDirs(i);
    m = fitMask & T.probeDirDeg == d;
    r = T.gProbeOverGPrefParent(m);
    g = T.gProbe(m);
    rows{i} = table(d, sum(m), 0, alpha, ...
      medianFinite(r), NaN, NaN, NaN, NaN, ...
      medianFinite(g), NaN, NaN, NaN, NaN, ...
      'VariableNames', {'probeDirDeg','nProbeSessions','nBootstrap','bootstrapAlpha', ...
      'medianGProbeOverGPref','ciLowGProbeOverGPref','ciHighGProbeOverGPref', ...
      'bootstrapSEGProbeOverGPref','bootstrapVarGProbeOverGPref', ...
      'medianGProbe','ciLowGProbe','ciHighGProbe', ...
      'bootstrapSEGProbe','bootstrapVarGProbe'});
  end
  B = vertcat(rows{:});
  return;
end

oldRng = rng;
cleanupObj = onCleanup(@() rng(oldRng));
if ~isempty(seed)
  rng(seed);
end

loP = alpha/2;
hiP = 1 - alpha/2;
for i = 1:numel(probeDirs)
  d = probeDirs(i);
  m = fitMask & T.probeDirDeg == d;
  r = T.gProbeOverGPrefParent(m);
  g = T.gProbe(m);
  r = r(isfinite(r));
  g = g(isfinite(g));
  n = numel(r);
  if n == 0
    bootR = NaN(nBoot,1);
  elseif n == 1
    bootR = repmat(r, nBoot, 1);
  else
    bootR = nan(nBoot,1);
    for b = 1:nBoot
      idx = randi(n, n, 1);
      bootR(b) = median(r(idx));
    end
  end
  ng = numel(g);
  if ng == 0
    bootG = NaN(nBoot,1);
  elseif ng == 1
    bootG = repmat(g, nBoot, 1);
  else
    bootG = nan(nBoot,1);
    for b = 1:nBoot
      idx = randi(ng, ng, 1);
      bootG(b) = median(g(idx));
    end
  end
  seR = stdFinite(bootR);
  seG = stdFinite(bootG);

  rows{i} = table(d, sum(m), nBoot, alpha, ...
    medianFinite(r), quantileLocal(bootR,loP), quantileLocal(bootR,hiP), seR, seR.^2, ...
    medianFinite(g), quantileLocal(bootG,loP), quantileLocal(bootG,hiP), seG, seG.^2, ...
    'VariableNames', {'probeDirDeg','nProbeSessions','nBootstrap','bootstrapAlpha', ...
    'medianGProbeOverGPref','ciLowGProbeOverGPref','ciHighGProbeOverGPref', ...
    'bootstrapSEGProbeOverGPref','bootstrapVarGProbeOverGPref', ...
    'medianGProbe','ciLowGProbe','ciHighGProbe', ...
    'bootstrapSEGProbe','bootstrapVarGProbe'});
end
B = vertcat(rows{:});
B = sortrows(B, 'probeDirDeg');
end

% -------------------------------------------------------------------------
function O = mergeBootstrapIntoOffsetSummary(O, B)
if isempty(O) || isempty(B)
  return;
end
[tf, loc] = ismember(O.probeDirDeg, B.probeDirDeg);
ciLowR = nan(height(O),1);
ciHighR = nan(height(O),1);
seR = nan(height(O),1);
varR = nan(height(O),1);
ciLowG = nan(height(O),1);
ciHighG = nan(height(O),1);
seG = nan(height(O),1);
varG = nan(height(O),1);
nBoot = nan(height(O),1);
alpha = nan(height(O),1);
ciLowR(tf) = B.ciLowGProbeOverGPref(loc(tf));
ciHighR(tf) = B.ciHighGProbeOverGPref(loc(tf));
seR(tf) = B.bootstrapSEGProbeOverGPref(loc(tf));

ciLowG(tf) = B.ciLowGProbe(loc(tf));
ciHighG(tf) = B.ciHighGProbe(loc(tf));
seG(tf) = B.bootstrapSEGProbe(loc(tf));

if ismember('bootstrapVarGProbeOverGPref', B.Properties.VariableNames)
  varR(tf) = B.bootstrapVarGProbeOverGPref(loc(tf));
else
  varR(tf) = seR(tf).^2;
end

if ismember('bootstrapVarGProbe', B.Properties.VariableNames)
  varG(tf) = B.bootstrapVarGProbe(loc(tf));
else
  varG(tf) = seG(tf).^2;
end

nBoot(tf) = B.nBootstrap(loc(tf));
alpha(tf) = B.bootstrapAlpha(loc(tf));

O.nBootstrap = nBoot;
O.bootstrapAlpha = alpha;

O.ciLowGProbeOverGPref = ciLowR;
O.ciHighGProbeOverGPref = ciHighR;
O.bootstrapSEGProbeOverGPref = seR;
O.bootstrapVarGProbeOverGPref = varR;

O.ciLowGProbe = ciLowG;
O.ciHighGProbe = ciHighG;
O.bootstrapSEGProbe = seG;
O.bootstrapVarGProbe = varG;
end

% -------------------------------------------------------------------------
function population = summarizePopulation(T, fitMask, O, failed, near)
g = T.gProbe(fitMask);
r = T.gProbeOverGPrefParent(fitMask);
dnll = T.deltaNLLProbe(fitMask);
p = T.pApproxProbe(fitMask);
nearMask = false(height(T),1);
if ismember('gProbeNearBound', T.Properties.VariableNames)
  nearMask = fitMask & T.gProbeNearBound;
end
probeDirs = unique(T.probeDirDeg(fitMask));

% Offset-level median ratios, giving each offset equal weight.
medByOffset = nan(numel(probeDirs),1);
for i = 1:numel(probeDirs)
  m = fitMask & T.probeDirDeg == probeDirs(i);
  medByOffset(i) = medianFinite(T.gProbeOverGPrefParent(m));
end

population = struct();
population.nProbeSessions = height(T);
population.nFitProbeSessions = sum(fitMask);
population.nFailedProbeSessions = height(failed);
population.nOffsets = numel(probeDirs);
population.nIncTrials = sum(T.nIncTrials(fitMask));
population.meanGProbe = meanFinite(g);
population.medianGProbe = medianFinite(g);
population.q25GProbe = quantileLocal(g,0.25);
population.q75GProbe = quantileLocal(g,0.75);
population.iqrGProbe = population.q75GProbe - population.q25GProbe;
population.meanGProbeOverGPrefParent = meanFinite(r);
population.medianGProbeOverGPrefParent = medianFinite(r);
population.q25GProbeOverGPrefParent = quantileLocal(r,0.25);
population.q75GProbeOverGPrefParent = quantileLocal(r,0.75);
population.iqrGProbeOverGPrefParent = population.q75GProbeOverGPrefParent - population.q25GProbeOverGPrefParent;
population.medianOffsetMedianGProbeOverGPref = medianFinite(medByOffset);
population.medianDeltaNLLProbe = medianFinite(dnll);
population.q25DeltaNLLProbe = quantileLocal(dnll,0.25);
population.q75DeltaNLLProbe = quantileLocal(dnll,0.75);
population.iqrDeltaNLLProbe = population.q75DeltaNLLProbe - population.q25DeltaNLLProbe;
population.totalDeltaNLLProbe = sum(dnll(isfinite(dnll)));
population.nDeltaNLLPositive = sum(dnll > 0);
population.fracDeltaNLLPositive = safeFrac(population.nDeltaNLLPositive, population.nFitProbeSessions);
population.nApproxP005 = sum(p < 0.05);
population.nApproxP001 = sum(p < 0.001);
population.nNearGainBound = sum(nearMask);
population.nNearBoundRowsWritten = height(near);

% Additional offset-level bookkeeping useful for status reports.
if isempty(O)
  population.minOffsetMedianGProbeOverGPref = NaN;
  population.maxOffsetMedianGProbeOverGPref = NaN;
else
  population.minOffsetMedianGProbeOverGPref = min(O.medianGProbeOverGPref);
  population.maxOffsetMedianGProbeOverGPref = max(O.medianGProbeOverGPref);
end
end

% -------------------------------------------------------------------------
function makeSummaryPDF(T, fitMask, O, P, pdfPath, animalTag)
fig = figure('Color','w', 'Visible','off', 'Units','inches', 'Position',[1 1 11 8.5]);
tiledlayout(fig, 2, 2, 'Padding','compact', 'TileSpacing','compact');

x = T.probeDirDeg(fitMask);
g = T.gProbe(fitMask);
r = T.gProbeOverGPrefParent(fitMask);
dnll = T.deltaNLLProbe(fitMask);

nexttile;
hold on;
plot(x, g, 'o', 'LineStyle','none', 'MarkerSize',4);
plot(O.probeDirDeg, O.medianGProbe, 'k-', 'LineWidth',1.2);
plot(O.probeDirDeg, O.q25GProbe, 'k:', 'LineWidth',0.9);
plot(O.probeDirDeg, O.q75GProbe, 'k:', 'LineWidth',0.9);
yline(0, ':');
yline(1, '--');
xlabel('Probe direction / offset (deg)');
ylabel('gProbe');
title('Probe gain by offset');
box off;

nexttile;
hold on;
plot(x, r, 'o', 'LineStyle','none', 'MarkerSize',4);
plot(O.probeDirDeg, O.medianGProbeOverGPref, 'k-', 'LineWidth',1.2);
plot(O.probeDirDeg, O.q25GProbeOverGPref, 'k:', 'LineWidth',0.9);
plot(O.probeDirDeg, O.q75GProbeOverGPref, 'k:', 'LineWidth',0.9);
if ismember('ciLowGProbeOverGPref', O.Properties.VariableNames)
  errorbar(O.probeDirDeg, O.medianGProbeOverGPref, ...
    O.medianGProbeOverGPref - O.ciLowGProbeOverGPref, ...
    O.ciHighGProbeOverGPref - O.medianGProbeOverGPref, ...
    'LineStyle','none', 'Color','k', 'LineWidth',0.8, 'CapSize',6, 'HandleVisibility','off');
end
yline(0, ':');
yline(1, '--');
xlabel('Probe direction / offset (deg)');
ylabel('gProbe / gPrefParent');
title('Probe gain relative to parent preferred gain');
box off;

nexttile;
hold on;
plot(x, dnll, 'o', 'LineStyle','none', 'MarkerSize',4);
plot(O.probeDirDeg, O.medianDeltaNLLProbe, 'k-', 'LineWidth',1.2);
yline(0, ':');
xlabel('Probe direction / offset (deg)');
ylabel('Delta NLL probe');
title('Likelihood improvement by offset');
box off;

nexttile;
axis off;
lines = {
  sprintf('Animal: %s', animalTag)
  sprintf('Probe sessions fit: %d/%d', P.nFitProbeSessions, P.nProbeSessions)
  sprintf('Failed/omitted sessions: %d', P.nFailedProbeSessions)
  sprintf('Offsets: %d', P.nOffsets)
  sprintf('INC trials: %d', P.nIncTrials)
  sprintf('Median gProbe: %.3f  [IQR %.3f, %.3f]', P.medianGProbe, P.q25GProbe, P.q75GProbe)
  sprintf('Median gProbe/gPref: %.3f  [IQR %.3f, %.3f]', P.medianGProbeOverGPrefParent, P.q25GProbeOverGPrefParent, P.q75GProbeOverGPrefParent)
  sprintf('Median offset-median gProbe/gPref: %.3f', P.medianOffsetMedianGProbeOverGPref)
  sprintf('Median Delta NLL: %.3f', P.medianDeltaNLLProbe)
  sprintf('Total Delta NLL: %.1f', P.totalDeltaNLLProbe)
  sprintf('Delta NLL > 0: %d/%d', P.nDeltaNLLPositive, P.nFitProbeSessions)
  sprintf('Approx p < 0.05: %d/%d', P.nApproxP005, P.nFitProbeSessions)
  sprintf('Approx p < 0.001: %d/%d', P.nApproxP001, P.nFitProbeSessions)
  sprintf('Near gain bound: %d/%d', P.nNearGainBound, P.nFitProbeSessions)
  };
text(0.02, 0.98, strjoin(lines, '\n'), 'Units','normalized', ...
  'VerticalAlignment','top', 'FontName','Menlo', 'FontSize',10);
title('Summary');

% try
  exportgraphics(fig, pdfPath, 'ContentType','vector');
% catch
%   print(fig, pdfPath, '-dpdf', '-bestfit');
% end
close(fig);
end

%%----------------------------------------------------------------
function fitSummary = makeIDRGainAcrossOffsetReadoutFits(varargin)
% makeIDRGainAcrossOffsetReadoutFits
%
% Adapter from IDR Weibull probe-gain offset summary to the generic
% fitAcrossOffsetReadout.m readout-function fitter.
%
% Required table variables:
%   probeDirDeg
%   medianGProbeOverGPref
%
% Preferred uncertainty variable:
%   bootstrapVarGProbeOverGPref
%
% Also accepted, if bootstrap variance is not yet present:
%   varGProbeOverGPref
%   semGProbeOverGPref
%   ci95LowGProbeOverGPref / ci95HighGProbeOverGPref
%   q25GProbeOverGPref / q75GProbeOverGPref  [fallback only]
%
% The 0-deg normalization anchor is added if absent:
%   measurement = 1
%   variance    = 0
%   nSessions   = Inf
% The 0-deg point is retained in the output but is ignored by
% fitAcrossOffsetReadout.m during fitting.

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'Bounds', struct(), @isstruct);
addParameter(p, 'SigmaMTDeg', 37.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'PhiDeg', -180:1:179, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)));


%% REMOVE THIS KLUDGE ONCE WE ARE EXCLUDING 1° PROBES
addParameter(p, 'MinProbeSessionsForFit', 2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);

parse(p, varargin{:});
opts = p.Results;

% -------------------------------------------------------------------------
% Load summary and find offsetSummary table
% -------------------------------------------------------------------------
gainSummaryFile = fullfile(domainFolder(mfilename('fullpath')), 'Data', 'AcrossOffsetSummaries', 'GainAnalysis', ...
  sprintf('IDRGainProbeOffsetSummary_%s.mat', opts.Animal));
L = load(gainSummaryFile);
offsetSummary = L.summaryData.offsetSummary;
vars = offsetSummary.Properties.VariableNames;

% -------------------------------------------------------------------------
% Extract measurements
% -------------------------------------------------------------------------
offsetsDeg   = offsetSummary.probeDirDeg(:)';
measurements = offsetSummary.medianGProbeOverGPref(:)';
nSessions = offsetSummary.nProbeSessions(:)';
variances = offsetSummary.bootstrapVarGProbeOverGPref(:)';

% -------------------------------------------------------------------------
% Add or validate 0-deg normalization anchor
% -------------------------------------------------------------------------
anchorIdx = find(abs(offsetsDeg) < 1e-9);

if isempty(anchorIdx)
    offsetsDeg   = [0, offsetsDeg];
    measurements = [1, measurements];
    variances    = [0, variances];
    nSessions    = [Inf, nSessions];

elseif isscalar(anchorIdx)
    % Treat the anchor as analytical. The fitter excludes it, but keeping it
    % exactly normalized avoids ambiguity in downstream plots/summaries.
    measurements(anchorIdx) = 1;
    variances(anchorIdx)    = 0;

else
    error('makeIDRGainAcrossOffsetReadoutFits:DuplicateAnchor', ...
        'More than one 0-deg row found in offsetSummary.');
end

% Sort for clean summaries and plots.
[offsetsDeg, sortIdx] = sort(offsetsDeg);
measurements = measurements(sortIdx);
variances    = variances(sortIdx);
nSessions    = nSessions(sortIdx);

% Basic validation. The 0-deg variance may be zero because it is not fit.
isAnchor = abs(offsetsDeg) < 1e-9;
badNonAnchor = ~isAnchor & ...
    (~isfinite(offsetsDeg) | ~isfinite(measurements) | ...
     ~isfinite(variances) | variances <= 0 | nSessions <= 0);

assert(~any(badNonAnchor), ...
    'makeIDRGainAcrossOffsetReadoutFits:BadInputRows', ...
    'One or more non-anchor offset rows have invalid measurement, variance, or nSessions.');

% Exclude offsets with too few contributing probe sessions from the readout fit,
% while preserving their measured values in the saved summary.
fitMeasurements = measurements;
fitVariances = variances;
excludeFromFit = ~isAnchor & nSessions < opts.MinProbeSessionsForFit;
fitMeasurements(excludeFromFit) = NaN;
fitVariances(excludeFromFit) = NaN;

% Call generic readout fitter
fitSummary = fitAcrossOffsetReadout(fitMeasurements, fitVariances, offsetsDeg, ...
    'NSessions', nSessions, 'Bounds', opts.Bounds, 'SigmaMTDeg', opts.SigmaMTDeg, ...
    'PhiDeg', opts.PhiDeg, 'SourceMeasureType', 'WeibullGain_gProbeOverGPref', ...
    'SourceSideType', 'probe', 'SourceStepType', 'INC', 'SourceMode', 'sessionBootstrapMedian');

fitSummary.type = 'Beta';
fitSummary.source = struct();
fitSummary.source.gainSummaryFile = gainSummaryFile;
fitSummary.source.rawMeasurements = struct( ...               %% this one entry is a kludge to be removed
  'offsetsDeg', offsetsDeg, ...
  'values', measurements, ...
  'variances', variances, ...
  'nSessions', nSessions, ...
  'excludeFromFit', excludeFromFit, ...
  'exclusionReason', sprintf('nSessions < %g', opts.MinProbeSessionsForFit));
fitSummary.source.adapterFunction = mfilename;
fitSummary.source.offsetSummaryVariable = 'medianGProbeOverGPref';
fitSummary.source.varianceInterpretation = ...
  'Variance of bootstrap distribution for median gProbe/gPrefParent by offset.';
fitSummary.source.anchorNote = ...
  '0-deg anchor set to measurement=1 and excluded from DOG fitting by fitAcrossOffsetReadout.';

% -------------------------------------------------------------------------
% Add compatibility fields for makeAcrossOffsetPlots.
% These fields mimic the historical acrossOffsetSummary layout used by
% updateAcrossOffsetSummaries, but they are derived from the gain-summary
% offset table rather than from kernel bootstraps.
% -------------------------------------------------------------------------

raw = fitSummary.source.rawMeasurements;

plotMask = ~isAnchor & ~excludeFromFit & isfinite(raw.values) & isfinite(raw.variances) & raw.variances > 0;
plotOffsets = raw.offsetsDeg(plotMask);
plotValues  = raw.values(plotMask);
plotVars    = raw.variances(plotMask);
plotNSess   = raw.nSessions(plotMask);

% Use bootstrap percentile CI fields from offsetSummary if present.
ciLow = nan(size(plotOffsets));
ciHigh = nan(size(plotOffsets));

if all(ismember({'ciLowGProbeOverGPref','ciHighGProbeOverGPref'}, vars))
  [tfCI, locCI] = ismember(plotOffsets, offsetSummary.probeDirDeg(:)');
  ciLow(tfCI)  = offsetSummary.ciLowGProbeOverGPref(locCI(tfCI));
  ciHigh(tfCI) = offsetSummary.ciHighGProbeOverGPref(locCI(tfCI));
else
  % Fallback: approximate 95% CI from bootstrap variance.
  ciLow  = plotValues - 1.96 .* sqrt(plotVars);
  ciHigh = plotValues + 1.96 .* sqrt(plotVars);
end

% Trial counts are not used deeply by the plots except for annotation.
% Preserve total INC trials where available, placing them in column 2 because
% stepType is INC.
nTrials = nan(size(plotOffsets));
if ismember('nIncTrials', vars)
  [tfN, locN] = ismember(plotOffsets, offsetSummary.probeDirDeg(:)');
  nTrials(tfN) = offsetSummary.nIncTrials(locN(tfN));
end

empirical = repmat(struct('probeOffsetDeg', NaN, 'stepType', 2, 'pooledScale', NaN, 'meanSessionScale', NaN, ...
  'medianSessionScale', NaN, 'semSessionScale', NaN, 'boot95', [NaN NaN], 'nSessions', 0, 'nTrials', NaN, ...
  'nTrialsByStep', [NaN NaN], 'sessionWeighting', 'session_bootstrap_median_gain_ratio', ...
  'distributionSkew', NaN, 'meanScale', NaN, 'medianScale', NaN, 'semScale', NaN), numel(plotOffsets), 1);

offsetBootstrap = repmat(struct('probeOffsetDeg', NaN, 'bootScale', [], 'bootMean', NaN, ...
  'bootMedian', NaN, 'boot68', [NaN NaN], 'boot95', [NaN NaN]), numel(plotOffsets), 1);

for k = 1:numel(plotOffsets)
  empirical(k).probeOffsetDeg = plotOffsets(k);
  empirical(k).stepType = 2;  % INC
  empirical(k).pooledScale = plotValues(k);
  empirical(k).meanSessionScale = NaN;
  empirical(k).medianSessionScale = plotValues(k);
  empirical(k).semSessionScale = sqrt(plotVars(k));
  empirical(k).boot95 = [ciLow(k), ciHigh(k)];
  empirical(k).nSessions = plotNSess(k);
  empirical(k).nTrials = nTrials(k);
  empirical(k).nTrialsByStep = [NaN, nTrials(k)];

  % Backward-compatible aliases.
  empirical(k).meanScale = empirical(k).pooledScale;
  empirical(k).medianScale = empirical(k).medianSessionScale;
  empirical(k).semScale = empirical(k).semSessionScale;

  offsetBootstrap(k).probeOffsetDeg = plotOffsets(k);
  offsetBootstrap(k).bootScale = [];
  offsetBootstrap(k).bootMean = NaN;
  offsetBootstrap(k).bootMedian = plotValues(k);
  offsetBootstrap(k).boot68 = [NaN NaN];
  offsetBootstrap(k).boot95 = [ciLow(k), ciHigh(k)];
end

fitSummary.empirical = empirical;

fitSummary.bootstrap = struct();
fitSummary.bootstrap.offsetBootstrap = offsetBootstrap;
fitSummary.bootstrap.fitBootstrap = struct( ...
  'modelName', 'gain_offset_variance_only', ...
  'offsetsDeg', plotOffsets, ...
  'offsetFitVar', plotVars, ...
  'offsetFitWeights', 1 ./ plotVars, ...
  'fitWeights', 1 ./ plotVars, ...
  'note', ['Compatibility field for makeAcrossOffsetPlots. Variances are ' ...
           'bootstrap variances of offset-level median gProbe/gPrefParent.']);

fitSummary.meta = struct();
fitSummary.meta.analysisDate = datetime('now');
fitSummary.meta.scaleMetric = 'WeibullGain_gProbeOverGPref';
fitSummary.meta.ciType = 'session_bootstrap_percentile';
fitSummary.meta.nBoot = NaN;
fitSummary.meta.bootstrapType = 'probe_session_bootstrap';
fitSummary.meta.angleUnits = 'deg';
fitSummary.meta.normalization = 'gProbeOverGPrefParent_with_0deg_anchor';
fitSummary.meta.offsetKeysDeg = plotOffsets;
fitSummary.meta.notes = 'Plot-compatible summary generated by makeIDRGainAcrossOffsetReadoutFits.';

opts.PlotDir = fullfile(domainFolder(mfilename('fullpath')), 'Plots', 'AcrossProbes', 'ReadoutFits', 'Beta');
opts.NBoot = max(offsetSummary.nBootstrap, [], 'omitnan');
makeAcrossOffsetPlots(fitSummary, opts);

outputFile = fullfile(domainFolder(mfilename('fullpath')), 'Data', 'AcrossOffsetSummaries', 'GainAnalysis', ...
  sprintf('IDRGainAcrossOffsetReadoutFits_%s.mat', opts.Animal));
save(outputFile, 'fitSummary');

end

% -------------------------------------------------------------------------
% function T = csvSafeDiagnosticTable(T)
% % Make diagnostic tables safe and readable in CSV output.
% % Some upstream fields can be cell arrays or multi-element strings inside a
% % table cell (for example parentFileName), and optimizer messages can contain
% % embedded newlines.  The MAT output preserves the full table; the CSV output
% % should be compact and human-scannable.
% if isempty(T)
%   return;
% end
% vars = T.Properties.VariableNames;
% for v = 1:numel(vars)
%   name = vars{v};
%   x = T.(name);
%   if iscell(x)
%     y = strings(height(T),1);
%     for i = 1:height(T)
%       y(i) = compactCellValue(x{i});
%     end
%     T.(name) = y;
%   elseif isstring(x)
%     T.(name) = compactStringVector(x);
%   elseif ischar(x)
%     T.(name) = compactStringVector(string(x));
%   end
% end
% end
% 
% % -------------------------------------------------------------------------
% function y = compactStringVector(x)
% x = string(x);
% y = strings(size(x));
% for i = 1:numel(x)
%   y(i) = compactOneString(x(i));
% end
% end
% 
% % -------------------------------------------------------------------------
% function s = compactCellValue(c)
% if isempty(c)
%   s = "";
% elseif ischar(c) || isstring(c)
%   x = string(c);
%   x = x(:);
%   if isscalar(x)
%     s = compactOneString(x);
%   else
%     % Prefer first element for fields such as parentFileName where a scalar
%     % identity was intended but a vector was accidentally stored upstream.
%     s = compactOneString(x(1));
%   end
% elseif isnumeric(c) || islogical(c)
%   if isscalar(c)
%     s = string(c);
%   else
%     s = sprintf('[%s %s]', class(c), mat2str(size(c)));
%   end
% else
%   s = sprintf('[%s]', class(c));
% end
% end
% 
% % -------------------------------------------------------------------------
% function s = compactOneString(s)
% s = string(s);
% if ismissing(s)
%   s = "";
%   return;
% end
% s = replace(s, newline, " ");
% s = replace(s, char(13), " ");
% s = regexprep(s, '\\s+', ' ');
% s = strtrim(s);
% end

% -------------------------------------------------------------------------
function T = rowsFromScalarStruct(S)
fields = fieldnames(S);
metric = strings(numel(fields),1);
value = nan(numel(fields),1);
for i = 1:numel(fields)
  metric(i) = string(fields{i});
  v = S.(fields{i});
  if isnumeric(v) || islogical(v)
    value(i) = double(v);
  else
    value(i) = NaN;
  end
end
T = table(metric, value);
end

% -------------------------------------------------------------------------
function T = selectExistingVariables(Tin, wanted)
vars = Tin.Properties.VariableNames;
keep = wanted(ismember(wanted, vars));
if isempty(keep)
  T = table();
else
  T = Tin(:, keep);
end
end

% -------------------------------------------------------------------------
function T = addFailureReason(T)
if isempty(T)
  return;
end
n = height(T);
failureReason = strings(n,1);
for i = 1:n
  msg = "";
  if ismember('message', T.Properties.VariableNames)
    msg = string(T.message(i));
  end
  msgProbe = "";
  if ismember('messageProbe', T.Properties.VariableNames)
    msgProbe = string(T.messageProbe(i));
  end
  combined = strtrim(strjoin([msg msgProbe], " | "));
  if strlength(combined) == 0 || ismissing(combined)
    combined = "fitIncluded false, no message";
  end
  failureReason(i) = combined;
end
T.failureReason = failureReason;
T = movevars(T, 'failureReason', 'Before', 1);
end

% -------------------------------------------------------------------------
function T = sortrowsIfPresent(T, vars)
if isempty(T)
  return;
end
vars = vars(ismember(vars, T.Properties.VariableNames));
if ~isempty(vars)
  T = sortrows(T, vars);
end
end

% -------------------------------------------------------------------------
function T = emptyFailedTable()
T = table(strings(0,1), strings(0,1), nan(0,1), strings(0,1), strings(0,1), ...
  'VariableNames', {'failureReason','animal','probeDirDeg','probeSessionName','message'});
end

% -------------------------------------------------------------------------
function T = emptyNearBoundTable()
T = table(strings(0,1), nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
  'VariableNames', {'probeSessionName','probeDirDeg','gProbe','gProbeOverGPrefParent','deltaNLLProbe'});
end

% -------------------------------------------------------------------------
function f = safeFrac(n, d)
if d <= 0
  f = NaN;
else
  f = n ./ d;
end
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
function s = stdFinite(x)
x = x(isfinite(x));
if numel(x) <= 1
  s = NaN;
else
  s = std(x);
end
end

% -------------------------------------------------------------------------
function q = quantileLocal(x, p)
x = sort(x(isfinite(x)));
if isempty(x)
  q = NaN;
  return;
end
if isscalar(x)
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
