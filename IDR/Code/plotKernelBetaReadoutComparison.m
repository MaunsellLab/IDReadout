function fig = plotKernelBetaReadoutComparison(varargin)
% plotKernelBetaReadoutComparison  Compare kernel- and beta-derived readouts.
%
% Creates one compact 2x2 figure:
%   top row: inferred readout functions for signed and rectified templates
%   bottom row: corresponding predicted normalized scale curves

baseFolder = domainFolder(mfilename('fullpath'));
defaultKernelFile = fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries', 'IDR_KernelSummary.mat');
defaultBetaFile = fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries', 'IDR_BetaSummary.mat');
defaultPlotDir = fullfile(baseFolder, 'Plots', 'AcrossProbes', 'ReadoutFits');

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'KernelSummaryFile', defaultKernelFile, @(x) ischar(x) || isstring(x));
addParameter(p, 'BetaSummaryFile', defaultBetaFile, @(x) ischar(x) || isstring(x));
addParameter(p, 'PlotDir', defaultPlotDir, @(x) ischar(x) || isstring(x));
addParameter(p, 'FileName', 'KernelVsBetaReadoutComparison.pdf', @(x) ischar(x) || isstring(x));
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
parse(p, varargin{:});
opts = p.Results;

K = load(char(opts.KernelSummaryFile), 'acrossOffsetSummary');
B = load(char(opts.BetaSummaryFile), 'betaSummary');
if ~isfield(K, 'acrossOffsetSummary') || ~isfield(B, 'betaSummary')
  error('plotKernelBetaReadoutComparison:MissingSummary', ...
    'Expected acrossOffsetSummary and betaSummary in the supplied files.');
end

kernelModels = K.acrossOffsetSummary.readoutModels;
betaModels = B.betaSummary.readoutFitSummary.readoutModels;
modelNames = {'signedDOG', 'rectifiedDOG'};
displayNames = {'Signed template', 'Rectified template'};

fig = figure('Color', 'w', 'Position', [100 100 1000 720], ...
  'Visible', char(lower(string(opts.Visible))));
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for j = 1:2
  name = modelNames{j};
  km = kernelModels.(name);
  bm = betaModels.(name);

  nexttile(j); hold on;
  hk = plot(km.phiDeg, km.readoutPhi, '-', 'LineWidth', 2);
  hb = plot(bm.phiDeg, bm.readoutPhi, '--', 'LineWidth', 2);
  yline(0, ':');
  xlabel('\phi (deg)');
  ylabel('a(\phi)');
  title(sprintf('%s: inferred readout', displayNames{j}));
  legend([hk hb], {'Kernel scales', 'Beta scales'}, 'Location', 'best');
  box off;

  nexttile(j + 2); hold on;
  hk2 = plot(km.plotOffsetsDeg, km.plotPredictedScale, '-', 'LineWidth', 2);
  hb2 = plot(bm.plotOffsetsDeg, bm.plotPredictedScale, '--', 'LineWidth', 2);
  yline(0, ':');
  xlabel('Probe direction offset (deg)');
  ylabel('Predicted normalized scale');
  title(sprintf('%s: predicted scale', displayNames{j}));
  xlim([0 180]);
  legend([hk2 hb2], {'Kernel fit', 'Beta fit'}, 'Location', 'best');
  box off;
end

sgtitle('Kernel- and beta-derived MT readout fits');
plotDir = char(opts.PlotDir);
if ~isempty(plotDir)
  if ~isfolder(plotDir)
    mkdir(plotDir);
  end
  exportgraphics(fig, fullfile(plotDir, char(opts.FileName)), ...
    'ContentType', 'vector');
end
end
