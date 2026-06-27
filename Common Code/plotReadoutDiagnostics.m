function fig = plotReadoutDiagnostics(readoutModel, varargin)
% plotReadoutDiagnostics  Plot a fitted readout, MT templates, and overlaps.
%
% fig = plotReadoutDiagnostics(readoutModel, Name, Value, ...)
%
% Required input
%   readoutModel  One fitted model summary, such as signedDOG or rectifiedDOG.
%
% Options
%   'NBoot'         bootstrap count shown in title (default NaN)
%   'PlotDir'       output folder; empty means do not save (default '')
%   'FileName'      output filename (default based on template mode)
%   'Visible'       'on' or 'off' (default 'on')
%   'FigureNumber'  optional numeric figure number (default [])
%   'TitlePrefix'   optional title prefix (default '')

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'readoutModel', @isstruct);
addParameter(p, 'NBoot', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'PlotDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'FileName', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
addParameter(p, 'FigureNumber', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'TitlePrefix', '', @(x) ischar(x) || isstring(x));
parse(p, readoutModel, varargin{:});
opts = p.Results;
opts.PlotDir = char(opts.PlotDir);
opts.FileName = char(opts.FileName);
opts.Visible = char(lower(string(opts.Visible)));
opts.TitlePrefix = char(opts.TitlePrefix);

rm = readoutModel;
if isfield(rm, 'templateMode')
  templateMode = char(rm.templateMode);
else
  templateMode = 'signed';
end
if ~isfield(rm, 'fit') || isempty(rm.fit) || ...
    ~isfield(rm.fit, 'fitSuccess') || ~rm.fit.fitSuccess
  fig = gobjects(0);
  return;
end
phiDeg = rm.phiDeg(:)';
aPhi = rm.readoutPhi(:)';
mtp = rm.mtForwardModelParams;
mtModel = makeMTReadoutForwardModel( ...
  'sigmaMTDeg', mtp.sigmaMTDeg, 'phiDeg', mtp.phiDeg);
offsetsDeg = [0, rm.fit.offsetsDeg];
nOffsets = numel(offsetsDeg);

deltaM = cell(1, nOffsets);
prodTerm = cell(1, nOffsets);
overlap = nan(1, nOffsets);
for i = 1:nOffsets
  deltaM{i} = mtReadoutTemplate(offsetsDeg(i), mtModel, ...
    'TemplateMode', templateMode);
  prodTerm{i} = aPhi .* deltaM{i};
  overlap(i) = sum(prodTerm{i});
end

if isempty(opts.FigureNumber)
  fig = figure('Color', 'w', 'Visible', opts.Visible, 'Position', [100 100 850 900], 'WindowStyle', 'docked');
else
  fig = figure(opts.FigureNumber);
  clf(fig);
  set(fig, 'Color', 'w', 'Visible', opts.Visible, 'Position', [100 100 850 900], 'WindowStyle', 'docked');
end
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
lineCol = lines(nOffsets);

% ---- Readout function ----
nexttile; hold on;
hFitReadout = plot(phiDeg, aPhi, 'k-', 'LineWidth', 2);
plot(phiDeg, zeros(size(phiDeg)), 'k:', 'LineWidth', 1);
xlabel('\phi (deg)');
ylabel('a(\phi)');
if isfinite(opts.NBoot)
  bootText = sprintf(', %d bootstraps', opts.NBoot);
else
  bootText = '';
end
if isempty(opts.TitlePrefix)
  titlePrefix = '';
else
  titlePrefix = [opts.TitlePrefix ': '];
end
title(sprintf('%sFitted DOG readout (%s template%s)', ...
  titlePrefix, templateMode, bootText), 'Interpreter', 'none');
legend(hFitReadout, {'Fitted readout a(\phi)'}, 'Location', 'northeast');

nParams = min(numel(rm.paramNames), numel(rm.params));
paramText = cell(nParams, 1);
for q = 1:nParams
  paramText{q} = sprintf('%s: %.4f', rm.paramNames{q}, rm.params(q));
end
if ~isempty(paramText)
  text(-100, 0.95, paramText, 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'Interpreter', 'none');
end
ylimits = ylim;
ylim([min(0.2, ylimits(1)), max(1.1, ylimits(2))]);
box off;

% ---- MT population templates ----
nexttile; hold on;
title(sprintf('MT population responses to probes (%s template; flat readout)', ...
  templateMode), 'Interpreter', 'none');
hTemplates = gobjects(1, nOffsets);
for i = 1:nOffsets
  hTemplates(i) = plot(phiDeg, deltaM{i}, '-', ...
    'Color', lineCol(i,:), 'LineWidth', 1.5);
end
yline(0, 'k:');
xlabel('\phi (deg)');
ylabel('\Delta m(\phi;\delta)');
legendLabels = arrayfun(@(d) sprintf('\\Delta m(\\phi;%g^\\circ)', d), ...
  offsetsDeg, 'UniformOutput', false);
legend(hTemplates, legendLabels, 'Location', 'best');
box off;

% ---- Weighted response contributions ----
nexttile; hold on;
title('Weighted population responses');
legendHandles = gobjects(1, nOffsets);
legendLabels = cell(1, nOffsets);
for i = 1:nOffsets
  legendHandles(i) = plot(phiDeg, prodTerm{i}, '-', ...
    'Color', lineCol(i,:), 'LineWidth', 2);
  legendLabels{i} = sprintf('%g^\\circ fit, <a,\\Delta m> = %.2g; S_{fit} %.2f', ...
    offsetsDeg(i), overlap(i), overlap(i) / overlap(1));
end
yline(0, 'k:');
xlabel('\phi (deg)');
ylabel('a(\phi)\Delta m(\phi;\delta)');
legend(legendHandles, legendLabels, 'Location', 'best');
box off;

if isempty(opts.FileName)
  opts.FileName = sprintf('ReadoutFunctions_%s.pdf', templateMode);
end
if ~isempty(opts.PlotDir)
  if ~isfolder(opts.PlotDir)
    mkdir(opts.PlotDir);
  end
  exportgraphics(fig, fullfile(opts.PlotDir, opts.FileName), ...
    'ContentType', 'vector');
end
end
