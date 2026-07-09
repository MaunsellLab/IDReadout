function betaSummary = plotAcrossOffsetBetaSummary(varargin)
% plotAcrossOffsetBetaSummary  Regenerate beta across-offset plots only.
%
% Loads the saved betaSummary and writes plots without refitting or
% repeating the hierarchical bootstrap.
%
% Name-value options:
%   'SummaryFile'      saved IDR_acrossOffsetBetaSummary.mat file
%   'PlotDir'          output directory
%   'MakeRatioPlot'    regenerate BetaRatiosByOffset.pdf (default true)
%   'MakeReadoutPlot'  regenerate BetaReadoutFit.pdf (default true)
%   'MakeDiagnosticPlots' generate shared readout diagnostics (default true)
%   'Visible'          figure visibility: 'on' or 'off' (default 'on')
%
% Example:
%   plotAcrossOffsetBetaSummary();

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'MakeRatioPlot', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'MakeReadoutPlot', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'MakeDiagnosticPlots', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
parse(p, varargin{:});
opts = p.Results;
opts.Visible = char(lower(string(opts.Visible)));

baseFolder = domainFolder(mfilename('fullpath'));
summaryFile = fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries', sprintf('BetaSummary_%s.mat', opts.Animal));
load(summaryFile, 'betaSummary');

% required = {'offsetFits','measurements','readoutFitSummary','meta'};
% missing = required(~isfield(betaSummary, required));
% if ~isempty(missing)
%   error('plotAcrossOffsetBetaSummary:MissingFields', 'betaSummary is missing: %s', strjoin(missing, ', '));
% end

plotDir = validFolder(fullfile(baseFolder, 'Plots', 'AcrossProbes', 'ReadoutFits', 'Beta'));
if opts.MakeRatioPlot
  fprintf('  %s\n', fullfile(plotDir, sprintf('BetaRatiosByOffset_%s.pdf', opts.Animal)));
  plotBetaRatiosByOffset(betaSummary, fullfile(plotDir, sprintf('BetaRatiosByOffset_%s.pdf', opts.Animal)), ...
        opts.Visible);
end
if opts.MakeReadoutPlot
  filePath = fullfile(plotDir, sprintf('BetaFits_%s.pdf', opts.Animal));
  fprintf('  %s\n', filePath);
  plotBetaReadoutFit(betaSummary, filePath, opts.Visible);
end
if opts.MakeDiagnosticPlots
  R = betaSummary.readoutFitSummary.readoutModels;
  if isfield(R, 'signedDOG')
    fprintf('  %s\n', fullfile(plotDir, sprintf('BetaReadoutFunctions_Signed_%s.pdf', opts.Animal)));
    plotReadoutDiagnostics(R.signedDOG, 'NBoot', bootstrapCount(betaSummary), 'PlotDir', plotDir, ...
      'FileName', sprintf('BetaReadoutFunctions_Signed_%s.pdf', opts.Animal),  'Visible', opts.Visible, ...
      'TitlePrefix', 'Beta', 'FigureNumber', 312);
  end
  if isfield(R, 'rectifiedDOG')
    fprintf('  %s\n', fullfile(plotDir, sprintf('BetaReadoutFunctions_Rectified_%s.pdf', opts.Animal)));
      plotReadoutDiagnostics(R.rectifiedDOG, 'NBoot', bootstrapCount(betaSummary), 'PlotDir', plotDir, ...
      'FileName', sprintf('BetaReadoutFunctions_Rectified_%s.pdf', opts.Animal), 'TitlePrefix', 'Beta', ...
      'FigureNumber', 313);
  end
end
end

% -------------------------------------------------------------------------
function plotBetaRatiosByOffset(S, savePath, visible)
F = S.offsetFits;
fig = figure(310); 
clf(fig);
set(fig, 'Color','w', 'Visible', visible, 'WindowStyle', 'docked');
hold on;

hSession = gobjects(1);
% hPool = gobjects(1);
for k = 1:numel(F)
  validIdx = abs(F(k).sessionBetaRatio) < 3.0;
  betaRatios = F(k).sessionBetaRatio(validIdx);
  n = numel(betaRatios);
  jitter = zeros(n,1);
  if n > 1
    jitter = linspace(-1.5,1.5,n)';
  end
  h = plot(F(k).probeOffsetDeg+jitter, betaRatios, 'o', 'LineStyle', 'none', 'MarkerSize', 4);
  h.MarkerFaceColor = h.Color;
  h.MarkerEdgeColor = 'black';
  plot(F(k).probeOffsetDeg, mean(betaRatios), 'o', 'LineStyle', 'none', 'MarkerSize', 8, ...
          'MarkerFaceColor', 'black');
  if k == 1
    hSession = h;
  end
end

yline(0,':');
yline(1,'--');
ylim([-1, 3]);
xlabel('Probe direction offset (deg)');
ylabel('\beta_{probe}/\beta_{pref}');
title(sprintf('%s session ratios and pooled shared scales (%d bootstraps)', ...
  upper(char(S.meta.stepType)), bootstrapCount(S)));
xticks([F.probeOffsetDeg]);
legend(hSession, {'Session ratio'}, 'Location','best');
box off;

exportgraphics(fig, savePath, 'ContentType', 'vector');
if strcmpi(visible, 'off')
  close(fig);
end
end

% -------------------------------------------------------------------------
function plotBetaReadoutFit(S, savePath, visible)
M = S.measurements;
R = S.readoutFitSummary.readoutModels;
F = S.offsetFits;

fig = figure(311);
clf(fig);
set(fig, 'Color', 'w', 'Position', [100 100 900 540], 'Visible', visible, 'WindowStyle', 'docked');
hold on;

ci = nan(numel(F),2);
for k = 1:numel(F)
  [ci(k,1), ci(k,2)] = offsetCI95(S, k);
end

lo = M.pooledScale(:) - ci(:,1);
hi = ci(:,2) - M.pooledScale(:);

hObs = errorbar(M.offsetsDeg(:), M.pooledScale(:), lo, hi, 'ko', ...
  'MarkerFaceColor','k', 'LineWidth',1.2, 'CapSize',8);
hh = hObs;
labels = {'Pooled beta scale (95% CI)'};

if isfield(R, 'signedDOG') && isfield(R.signedDOG, 'fit') && ~isempty(R.signedDOG.fit) && R.signedDOG.fit.fitUsable
  h = plot(R.signedDOG.plotOffsetsDeg, R.signedDOG.plotPredictedScale, '-', 'LineWidth',1.5);
  hh(end+1) = h; 
  labels{end+1} = 'Signed DOG'; 
end

if isfield(R, 'rectifiedDOG') && isfield(R.rectifiedDOG,'fit') && ...
    ~isempty(R.rectifiedDOG.fit) && R.rectifiedDOG.fit.fitUsable
  h = plot(R.rectifiedDOG.plotOffsetsDeg, R.rectifiedDOG.plotPredictedScale, '-', 'LineWidth',1.5);
  hh(end+1) = h; 
  labels{end+1} = 'Rectified DOG'; 
  DOGFitText(0.98, 0.02, 'Rectified DOG', R.rectifiedDOG);
end

yline(0,':');
xlabel('Probe direction offset (deg)');
ylabel('Normalized beta scale');
title(sprintf('%s pooled beta scales and MT/readout fit (%d bootstraps)', ...
  upper(char(S.meta.stepType)), bootstrapCount(S)));
xlim([0 180]);
if max([lo; hi]) > 2.0
  ylim([-1, 2]);
end
legend(hh, labels, 'Location','southwest');
box off;

% Match the kernel-scale plot: pooled scale, percentile 95% CI, and trials.
textLines = cell(numel(F),1);
for k = 1:numel(F)
  textLines{k} = sprintf('%.0f°: scale %.2f, %.2f to %.2f 95%% CI (n = %d)', ...
    F(k).probeOffsetDeg, F(k).scale, ci(k,1), ci(k,2), F(k).nTrials);
end
annotation(fig, 'textbox', [0.55 0.69 0.42 0.24], 'String', textLines, 'Interpreter','none', ...
  'FitBoxToText','on', 'BackgroundColor','w', 'EdgeColor',[0.75 0.75 0.75], 'FontSize',9);
exportgraphics(fig, savePath, 'ContentType', 'vector');
end

% -------------------------------------------------------------------------
function [ciLow, ciHigh] = offsetCI95(S, offsetIndex)
F = S.offsetFits(offsetIndex);

if isfield(S, 'bootstrap') && isfield(S.bootstrap, 'bootScaleMat') && ...
    size(S.bootstrap.bootScaleMat, 2) >= offsetIndex
  x = S.bootstrap.bootScaleMat(:, offsetIndex);
  x = x(isfinite(x));
  if ~isempty(x)
    q = prctile(x, [2.5 97.5]);
    ciLow = q(1);
    ciHigh = q(2);
    return;
  end
end

if isfield(F,'boot95') && numel(F.boot95) == 2 && all(isfinite(F.boot95))
  ciLow = F.boot95(1);
  ciHigh = F.boot95(2);
elseif isfield(F,'scaleHessianSE') && isfinite(F.scaleHessianSE)
  ciLow = F.scale - 1.96*F.scaleHessianSE;
  ciHigh = F.scale + 1.96*F.scaleHessianSE;
else
  ciLow = F.scale;
  ciHigh = F.scale;
end
end

% ========================================================================
function nBoot = bootstrapCount(S)
if isfield(S,'meta') && isfield(S.meta,'nBoot') && isfinite(S.meta.nBoot)
  nBoot = S.meta.nBoot;
elseif isfield(S,'bootstrap') && isfield(S.bootstrap,'bootScaleMat')
  nBoot = size(S.bootstrap.bootScaleMat,1);
else
  nBoot = 0;
end
end

% ========================================================================
function DOGFitText(x, y, label, rm)
% One-model parameter/goodness-of-fit text block.

lines = {label};
for p = 1:min(numel(rm.params), numel(rm.paramNames))
  lines{end+1} = sprintf('  %s = %.4g', rm.paramNames{p}, rm.params(p)); %#ok<AGROW>
end
if isfield(rm, 'fit') && ~isempty(rm.fit) && isfield(rm.fit, 'goodnessOfFit') && isstruct(rm.fit.goodnessOfFit)
  g = rm.fit.goodnessOfFit;
  if isfield(g, 'weightedLoss') && isfinite(g.weightedLoss)
    lines{end+1} = sprintf('  loss = %.4g', g.weightedLoss); 
  end
  if isfield(g, 'reducedChiSq') && isfinite(g.reducedChiSq)
    lines{end+1} = sprintf('  red chi2 = %.4g', g.reducedChiSq); 
  end
  if isfield(g, 'aicc') && isfinite(g.aicc)
    lines{end+1} = sprintf('  AICc = %.4g', g.aicc); 
  end
end
txt = strjoin(lines, newline);
text(x, y, txt, 'Units', 'normalized', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8, ...
  'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 4, 'Interpreter', 'none');
end
