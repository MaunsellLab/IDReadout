function summaryData = makeIDRGainAcrossSessionPrefSummary(varargin)
% makeIDRGainAcrossSessionPrefSummary_v2
% Summarize and plot the stabilized IDR preferred-noise gain fit.
%
% Reads:
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainAcrossSessionPrefFit_<Animal>_v1.mat
%
% Created by:
%   makeIDRGainAcrossSessionPrefFit_v1.m
%
% The no-gain reference keeps the fitted shared beta/lapse fixed and refits
% alpha separately for each session with gPref fixed at zero. This provides a
% clean likelihood improvement for adding one session-specific gain term.
%
% Outputs:
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainAcrossSessionPrefSummary_<Animal>_v2.mat
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainAcrossSessionPrefSummary_<Animal>_v2_sessionTable.csv
%   Plots/AcrossProbes/GainAnalysis/IDRGainAcrossSessionPrefSummary_<Animal>_v2.pdf
%
% Usage:
%   S = makeIDRGainAcrossSessionPrefSummary_v2();
%   S = makeIDRGainAcrossSessionPrefSummary_v2('Animal','Meetz');
%   S = makeIDRGainAcrossSessionPrefSummary_v2('Replace',true);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'Replace', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'MakePlot', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opts = p.Results;

root = domainFolder(mfilename('fullpath'));
animalTag = char(string(opts.Animal));
gainFolder = validFolder(fullfile(root, 'Data', 'AcrossOffsetSummaries', 'GainAnalysis'));
plotFolder = validFolder(fullfile(root, 'Plots', 'AcrossProbes', 'GainAnalysis'));

fitPath = fullfile(gainFolder, sprintf('IDRGainAcrossSessionPrefFit_%s_v1.mat', animalTag));
summaryPath = fullfile(gainFolder, sprintf('IDRGainAcrossSessionPrefSummary_%s_v2.mat', animalTag));
csvPath = fullfile(gainFolder, sprintf('IDRGainAcrossSessionPrefSummary_%s_v2_sessionTable.csv', animalTag));
popCsvPath = fullfile(gainFolder, sprintf('IDRGainAcrossSessionPrefSummary_%s_v2_population.csv', animalTag));
plotPath = fullfile(plotFolder, sprintf('IDRGainAcrossSessionPrefSummary_%s_v2.pdf', animalTag));

if isfile(summaryPath) && ~opts.Replace
  if opts.Verbose
    fprintf('IDR preferred-gain summary already exists; loading %s\n', summaryPath);
  end
  S = load(summaryPath, 'summaryData');
  summaryData = S.summaryData;
  return;
end

if ~isfile(fitPath)
  error('makeIDRGainAcrossSessionPrefSummary:MissingFit', ...
    'Could not find %s. Run makeIDRGainAcrossSessionPrefFit_v1 first.', fitPath);
end
F = load(fitPath, 'fitData');
fitData = F.fitData;

requiredFitFields = {'params', 'sessionTable', 'trialTable'};
missing = requiredFitFields(~isfield(fitData, requiredFitFields));
if ~isempty(missing)
  error('makeIDRGainAcrossSessionPrefSummary:BadFitData', ...
    'fitData is missing fields: %s', strjoin(missing, ', '));
end

P = fitData.params;
Ttrial = fitData.trialTable;
Tsession = fitData.sessionTable;

betaShared = double(P.betaShared);
lapseShared = double(P.lapseShared);
nSessions = height(Tsession);

alphaNoGain = nan(nSessions, 1);
nllNoGain = nan(nSessions, 1);
exitflagNoGain = nan(nSessions, 1);
messageNoGain = strings(nSessions, 1);

for iSession = 1:nSessions
  m = Ttrial.sessionIndex == iSession;
  f = fitNoGainAlpha(Ttrial.deltaC(m), logical(Ttrial.correct(m)), ...
    betaShared, lapseShared, Tsession.alphaSession(iSession));
  alphaNoGain(iSession) = f.alpha;
  nllNoGain(iSession) = f.nll;
  exitflagNoGain(iSession) = f.exitflag;
  messageNoGain(iSession) = string(f.message);
end

sessionSummary = Tsession;
sessionSummary.alphaNoGain = alphaNoGain;
sessionSummary.nllNoGain = nllNoGain;
sessionSummary.exitflagNoGain = exitflagNoGain;
sessionSummary.messageNoGain = messageNoGain;
sessionSummary.deltaNLLPref = nllNoGain - sessionSummary.nllContribution;
sessionSummary.deltaDeviancePref = 2 * sessionSummary.deltaNLLPref;
sessionSummary.pApproxGain = erfc(sqrt(max(sessionSummary.deltaNLLPref, 0)));
sessionSummary.nllImprovementPerTrial = sessionSummary.deltaNLLPref ./ sessionSummary.nValidTrials;
sessionSummary.alphaGainOverNoGain = sessionSummary.alphaSession ./ sessionSummary.alphaNoGain;

if isfield(fitData, 'metadata') && isfield(fitData.metadata, 'options') && ...
    isfield(fitData.metadata.options, 'GainLimit')
  gainLimit = fitData.metadata.options.GainLimit;
else
  gainLimit = NaN;
end
if isfinite(gainLimit)
  sessionSummary.gPrefNearBound = abs(sessionSummary.gPrefSession) > 0.95 * gainLimit;
else
  sessionSummary.gPrefNearBound = false(nSessions, 1);
end

population = struct();
population.nSessions = nSessions;
population.nTrials = height(Ttrial);
population.betaShared = betaShared;
population.lapseShared = lapseShared;
population.totalNLLGain = fitData.nll;
population.totalNLLNoGain = sum(nllNoGain);
population.totalDeltaNLLPref = population.totalNLLNoGain - population.totalNLLGain;
population.totalDeltaDeviancePref = 2 * population.totalDeltaNLLPref;
% Approximate aggregate p value for adding nSessions gain terms. For large
% df this is only a rough diagnostic unless the model regularity assumptions
% are accepted. Avoid requiring the Statistics Toolbox.
population.meanGPref = mean(sessionSummary.gPrefSession, 'omitnan');
population.medianGPref = median(sessionSummary.gPrefSession, 'omitnan');
population.sdGPref = std(sessionSummary.gPrefSession, 'omitnan');
population.fracPositiveGPref = mean(sessionSummary.gPrefSession > 0, 'omitnan');
population.meanDeltaNLLPref = mean(sessionSummary.deltaNLLPref, 'omitnan');
population.medianDeltaNLLPref = median(sessionSummary.deltaNLLPref, 'omitnan');
population.fracPositiveDeltaNLLPref = mean(sessionSummary.deltaNLLPref > 0, 'omitnan');
population.nNearGainBound = sum(sessionSummary.gPrefNearBound);
population.nPositiveGPref = sum(sessionSummary.gPrefSession > 0, 'omitnan');
population.nNegativeGPref = sum(sessionSummary.gPrefSession < 0, 'omitnan');
population.nZeroOrNegativeGPref = sum(sessionSummary.gPrefSession <= 0, 'omitnan');
population.q25GPref = localQuantile(sessionSummary.gPrefSession, 0.25);
population.q75GPref = localQuantile(sessionSummary.gPrefSession, 0.75);
population.iqrGPref = population.q75GPref - population.q25GPref;
population.minGPref = min(sessionSummary.gPrefSession, [], 'omitnan');
population.maxGPref = max(sessionSummary.gPrefSession, [], 'omitnan');
population.nDeltaNLLPositive = sum(sessionSummary.deltaNLLPref > 0, 'omitnan');
population.nDeltaNLLNegative = sum(sessionSummary.deltaNLLPref < 0, 'omitnan');
population.nApproxP005 = sum(sessionSummary.pApproxGain < 0.05, 'omitnan');
population.nApproxP001 = sum(sessionSummary.pApproxGain < 0.01, 'omitnan');
population.q25DeltaNLLPref = localQuantile(sessionSummary.deltaNLLPref, 0.25);
population.q75DeltaNLLPref = localQuantile(sessionSummary.deltaNLLPref, 0.75);
population.iqrDeltaNLLPref = population.q75DeltaNLLPref - population.q25DeltaNLLPref;
population.minDeltaNLLPref = min(sessionSummary.deltaNLLPref, [], 'omitnan');
population.maxDeltaNLLPref = max(sessionSummary.deltaNLLPref, [], 'omitnan');
population.medianNLLImprovementPerTrial = median(sessionSummary.nllImprovementPerTrial, 'omitnan');
population.medianAlphaGainOverNoGain = median(sessionSummary.alphaGainOverNoGain, 'omitnan');
population.iqrAlphaGainOverNoGain = localQuantile(sessionSummary.alphaGainOverNoGain, 0.75) - ...
  localQuantile(sessionSummary.alphaGainOverNoGain, 0.25);

populationTable = makePopulationTable(population);

summaryData = struct();
summaryData.version = 1;
summaryData.analysisName = 'IDRGainAcrossSessionPrefSummary';
summaryData.fitPath = fitPath;
summaryData.summaryPath = summaryPath;
summaryData.csvPath = csvPath;
summaryData.populationCsvPath = popCsvPath;
summaryData.plotPath = plotPath;
summaryData.population = population;
summaryData.populationTable = populationTable;
summaryData.sessionSummary = sessionSummary;
summaryData.fitDataMetadata = fitData.metadata;
summaryData.createdBy = mfilename;
summaryData.createdDate = datetime('now');

save(summaryPath, 'summaryData', '-v7.3');
writetable(sessionSummary, csvPath);
writetable(populationTable, popCsvPath);

if opts.MakePlot
  fig = plotSummary(summaryData);
  exportgraphics(fig, plotPath, 'ContentType', 'vector');
  close(fig);
end

if opts.Verbose
  fprintf('\nPreferred-gain stabilized fit summary:\n');
  fprintf('  sessions:               %d\n', population.nSessions);
  fprintf('  trials:                 %d\n', population.nTrials);
  fprintf('  beta shared:            %.5g\n', population.betaShared);
  fprintf('  lapse shared:           %.5g\n', population.lapseShared);
  fprintf('  median gPref:           %.5g\n', population.medianGPref);
  fprintf('  fraction gPref > 0:     %.3f\n', population.fracPositiveGPref);
  fprintf('  median delta NLL:       %.5g\n', population.medianDeltaNLLPref);
  fprintf('  fraction delta NLL > 0: %.3f\n', population.fracPositiveDeltaNLLPref);
  fprintf('  total delta NLL:        %.5g\n', population.totalDeltaNLLPref);
  fprintf('Saved MAT: %s\n', summaryPath);
  fprintf('Saved CSV: %s\n', csvPath);
  if opts.MakePlot
    fprintf('Saved PDF: %s\n', plotPath);
  end
end
end

% -------------------------------------------------------------------------
function q = localQuantile(x, p)
x = sort(x(isfinite(x(:))));
if isempty(x)
  q = NaN;
  return;
end
if numel(x) == 1
  q = x;
  return;
end
pos = 1 + (numel(x) - 1) * p;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
  q = x(lo);
else
  q = x(lo) + (pos - lo) * (x(hi) - x(lo));
end
end

% -------------------------------------------------------------------------
function T = makePopulationTable(P)
metric = strings(0,1);
value = nan(0,1);
append('nSessions', P.nSessions);
append('nTrials', P.nTrials);
append('betaShared', P.betaShared);
append('lapseShared', P.lapseShared);
append('meanGPref', P.meanGPref);
append('medianGPref', P.medianGPref);
append('q25GPref', P.q25GPref);
append('q75GPref', P.q75GPref);
append('iqrGPref', P.iqrGPref);
append('minGPref', P.minGPref);
append('maxGPref', P.maxGPref);
append('nPositiveGPref', P.nPositiveGPref);
append('fracPositiveGPref', P.fracPositiveGPref);
append('medianDeltaNLLPref', P.medianDeltaNLLPref);
append('q25DeltaNLLPref', P.q25DeltaNLLPref);
append('q75DeltaNLLPref', P.q75DeltaNLLPref);
append('iqrDeltaNLLPref', P.iqrDeltaNLLPref);
append('minDeltaNLLPref', P.minDeltaNLLPref);
append('maxDeltaNLLPref', P.maxDeltaNLLPref);
append('nDeltaNLLPositive', P.nDeltaNLLPositive);
append('fracPositiveDeltaNLLPref', P.fracPositiveDeltaNLLPref);
append('nApproxP005', P.nApproxP005);
append('nApproxP001', P.nApproxP001);
append('nNearGainBound', P.nNearGainBound);
append('totalNLLGain', P.totalNLLGain);
append('totalNLLNoGain', P.totalNLLNoGain);
append('totalDeltaNLLPref', P.totalDeltaNLLPref);
append('totalDeltaDeviancePref', P.totalDeltaDeviancePref);
append('medianNLLImprovementPerTrial', P.medianNLLImprovementPerTrial);
append('medianAlphaGainOverNoGain', P.medianAlphaGainOverNoGain);
append('iqrAlphaGainOverNoGain', P.iqrAlphaGainOverNoGain);
T = table(metric, value);

  function append(name, val)
    metric(end+1,1) = string(name); %#ok<AGROW>
    value(end+1,1) = double(val); %#ok<AGROW>
  end
end

% -------------------------------------------------------------------------
function fit = fitNoGainAlpha(deltaC, correct, betaShared, lapseShared, alpha0)
theta0 = log(max(double(alpha0), eps));
objective = @(theta) noGainAlphaNLL(theta, deltaC, correct, betaShared, lapseShared);
options = optimset('Display', 'off', 'MaxIter', 1000, 'MaxFunEvals', 3000, ...
  'TolX', 1e-8, 'TolFun', 1e-8);
[thetaHat, nll, exitflag, output] = fminsearch(objective, theta0, options);
fit = struct();
fit.alpha = exp(thetaHat);
fit.nll = nll;
fit.exitflag = exitflag;
fit.message = output.message;
end

% -------------------------------------------------------------------------
function nll = noGainAlphaNLL(theta, deltaC, correct, betaShared, lapseShared)
alpha = exp(theta(1));
p = weibullPcorrect(deltaC, alpha, betaShared, lapseShared);
nll = bernoulliNLL(p, correct);
end

% -------------------------------------------------------------------------
function p = weibullPcorrect(deltaC, alpha, beta, lapse)
p = 0.5 + (0.5 - lapse) .* (1 - exp(-((deltaC ./ alpha) .^ beta)));
end

% -------------------------------------------------------------------------
function nll = bernoulliNLL(p, correct)
epsP = 1e-12;
p = min(max(p, epsP), 1 - epsP);
correct = double(correct(:));
nll = -sum(correct .* log(p(:)) + (1 - correct) .* log(1 - p(:)));
if ~isfinite(nll)
  nll = realmax;
end
end

% -------------------------------------------------------------------------
function fig = plotSummary(S)
T = S.sessionSummary;
P = S.population;
fig = figure('Color', 'w', 'WindowStyle', 'docked', 'Position', [100 100 1100 850]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
histogram(T.gPrefSession, 'BinMethod', 'fd');
xline(0, ':');
xlabel('gPrefSession');
ylabel('Sessions');
title(sprintf('Preferred gain, median %.3g', P.medianGPref));
box off;

nexttile;
histogram(T.deltaNLLPref, 'BinMethod', 'fd');
xline(0, ':');
xlabel('NLL(no gain) - NLL(gain)');
ylabel('Sessions');
title(sprintf('Likelihood improvement, median %.3g', P.medianDeltaNLLPref));
box off;

nexttile;
scatter(T.nPreferredNoiseTrials, T.gPrefSession, 20, 'filled');
yline(0, ':');
xlabel('Preferred-noise trials');
ylabel('gPrefSession');
title('Gain versus preferred-noise trial count');
box off;

nexttile;
scatter(T.prefNoiseSD, T.gPrefSession, 20, 'filled');
yline(0, ':');
xlabel('SD of weighted preferred noise (% coh)');
ylabel('gPrefSession');
title('Gain versus predictor spread');
box off;

nexttile;
scatter(T.alphaNoGain, T.alphaSession, 20, 'filled');
hold on;
lims = axis;
lo = min(lims([1 3]));
hi = max(lims([2 4]));
plot([lo hi], [lo hi], ':');
xlim([lo hi]); ylim([lo hi]);
xlabel('alpha, no-gain reference');
ylabel('alpha, gain fit');
title('Alpha stability');
box off;

nexttile;
axis off;
text(0.02, 0.95, sprintf('sessions: %d', P.nSessions), 'Units', 'normalized');
text(0.02, 0.85, sprintf('trials: %d', P.nTrials), 'Units', 'normalized');
text(0.02, 0.75, sprintf('betaShared: %.4g', P.betaShared), 'Units', 'normalized');
text(0.02, 0.65, sprintf('lapseShared: %.4g', P.lapseShared), 'Units', 'normalized');
text(0.02, 0.55, sprintf('mean gPref: %.4g', P.meanGPref), 'Units', 'normalized');
text(0.02, 0.45, sprintf('median gPref: %.4g [%.4g, %.4g]', ...
  P.medianGPref, P.q25GPref, P.q75GPref), 'Units', 'normalized');
text(0.02, 0.35, sprintf('gPref > 0: %d/%d (%.3f)', ...
  P.nPositiveGPref, P.nSessions, P.fracPositiveGPref), 'Units', 'normalized');
text(0.02, 0.25, sprintf('median delta NLL: %.4g [%.4g, %.4g]', ...
  P.medianDeltaNLLPref, P.q25DeltaNLLPref, P.q75DeltaNLLPref), 'Units', 'normalized');
text(0.02, 0.15, sprintf('p < 0.05: %d; near bound: %d', ...
  P.nApproxP005, P.nNearGainBound), 'Units', 'normalized');
text(0.02, 0.05, sprintf('total delta NLL: %.4g', P.totalDeltaNLLPref), 'Units', 'normalized');
title('Summary');
end
