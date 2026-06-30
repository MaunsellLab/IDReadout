function fig = plotKernelBetaReadoutComparison(varargin)
% plotKernelBetaReadoutComparison  Compare kernel- and beta-derived readouts.
%
% Creates one compact 2x2 figure:
%   top row: inferred readout functions for signed and rectified templates
%   bottom row: corresponding predicted normalized scale curves

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
parse(p, varargin{:});
opts = p.Results;

baseFolder = domainFolder(mfilename('fullpath'));
kernelFile = fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries', sprintf('KernelSummary_%s.mat', opts.Animal));
betaFile = fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries', sprintf('BetaSummary_%s.mat', opts.Animal));

K = load(char(kernelFile), 'acrossOffsetSummary');
B = load(char(betaFile), 'betaSummary');
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

plotDir = validFolder(fullfile(baseFolder, 'Plots', 'AcrossProbes', 'ReadoutFits'));
opts.FileName = sprintf('KernelVsBetaReadoutComparison_%s.pdf', opts.Animal);

exportgraphics(fig, fullfile(plotDir, char(sprintf('KernelVsBetaReadoutComparison_%s.pdf', opts.Animal))), ...
    'ContentType', 'vector');
end
