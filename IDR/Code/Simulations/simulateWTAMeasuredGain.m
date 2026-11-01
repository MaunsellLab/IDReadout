function results = simulateWTAMeasuredGain(varargin)
% simulateWTAMeasuredGain  Show how hard-max selection distorts kernel gain.
%
% Two RDKs (change and no-change) each contribute n independent Gaussian
% readout mechanisms. All mechanisms have the same underlying input weight,
% fixed at one. One change-side mechanism (the drift mechanism) receives a
% mean evidence advantage chosen to give a specified behavioral hit rate.
%
% The behavioral response is "change" when
%
%   max(change-side mechanisms) > max(no-change-side mechanisms).
%
% K0 is the local derivative of hit probability with respect to evidence
% added to the advantaged drift mechanism. Ktheta is the same derivative for
% one unadvantaged change-side mechanism. For weak zero-mean experimental
% noise, these derivatives are the limiting regression/kernel gains.
%
% Crucially, every underlying weight equals one. Therefore departures of
% Ktheta/K0 from one are produced entirely by WTA selection.
%
% Example:
%   R = simulateWTAMeasuredGain;
%
% Options:
%   'NValues'          mechanisms per RDK (default 2:20)
%   'HitRates'         requested hit rates (default 0.50:0.10:1.00)
%   'InternalNoiseSD'  SD of independent Gaussian internal noise (default 1)
%   'DerivativeStep'   central-difference evidence step (default 0.01 SD)
%   'NearPerfectRate'  rate used to represent an exact request of 1 (default .999)
%   'SelectedN'        n values emphasized in line plots (default [2 5 10 20])
%   'SavePDF'          save the figure (default true)
%   'PDFPath'          output PDF path (default WTAMeasuredGainSimulation.pdf)
%   'FigureVisible'    'on' or 'off' (default 'on')
%
% Exact 100%% performance has zero local gain for every mechanism, making
% Ktheta/K0 equal to 0/0. A requested hit rate of 1 is therefore evaluated
% at NearPerfectRate but displayed at 100%% on the x-axis.

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'NValues', 2:20, @(x) isnumeric(x) && isvector(x) && all(x >= 2) && all(mod(x,1)==0));
addParameter(p, 'HitRates', 0.50:0.10:1.00, @(x) isnumeric(x) && isvector(x) && all(x >= 0.5) && all(x <= 1));
addParameter(p, 'InternalNoiseSD', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'DerivativeStep', 0.01, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'NearPerfectRate', 0.999, @(x) isnumeric(x) && isscalar(x) && x > 0.99 && x < 1);
addParameter(p, 'SelectedN', [2 5 10 20], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'SavePDF', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'PDFPath', fullfile(pwd, 'WTAMeasuredGainSimulation.pdf'), @(x) ischar(x) || isstring(x));
addParameter(p, 'FigureVisible', 'on', @(x) any(strcmpi(string(x), ["on","off"])));
parse(p, varargin{:});
opts = p.Results;

nValues = unique(double(opts.NValues(:)'));
requestedHitRates = double(opts.HitRates(:)');
modelHitRates = requestedHitRates;
modelHitRates(modelHitRates == 1) = opts.NearPerfectRate;
selectedN = intersect(double(opts.SelectedN(:)'), nValues, 'stable');
sigma = double(opts.InternalNoiseSD);
h = double(opts.DerivativeStep) * sigma;

nN = numel(nValues);
nH = numel(modelHitRates);
mu0 = nan(nN, nH);
K0 = nan(nN, nH);
Ktheta = nan(nN, nH);
realizedHitRate = nan(nN, nH);

for iN = 1:nN
  n = nValues(iN);
  for iH = 1:nH
    targetHit = modelHitRates(iH);
    mu = calibrateDriftAdvantage(n, targetHit, sigma);
    mu0(iN,iH) = mu;
    realizedHitRate(iN,iH) = hitProbability(n, mu, 0, sigma);

    K0(iN,iH) = ( ...
      hitProbability(n, mu+h, 0, sigma) - ...
      hitProbability(n, mu-h, 0, sigma)) / (2*h);

    Ktheta(iN,iH) = ( ...
      hitProbability(n, mu, h, sigma) - ...
      hitProbability(n, mu, -h, sigma)) / (2*h);
  end
end

kernelRatio = Ktheta ./ K0;

% The simple occupancy approximation discussed in the accompanying analysis:
% pmax(theta)=(1-hit)/n and pmax(0)=(2*hit-1)+(1-hit)/n.
occupancyRatio = nan(nN,nH);
for iN = 1:nN
  n = nValues(iN);
  pTheta = (1-modelHitRates) ./ n;
  p0 = (2*modelHitRates-1) + pTheta;
  occupancyRatio(iN,:) = pTheta ./ p0;
end

fig = figure('Color','w', 'Visible', char(opts.FigureVisible), ...
  'Name','WTA measured-gain simulation');
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
colors = lines(max(numel(selectedN),1));
x = 100*requestedHitRates;

% A: relative measured gain.
ax1 = nexttile(tl,1);
hold(ax1,'on');
for k = 1:numel(selectedN)
  row = find(nValues == selectedN(k), 1);
  plot(ax1, x, kernelRatio(row,:), '-o', 'Color', colors(k,:), ...
    'LineWidth',1.8, 'MarkerSize',4, 'DisplayName',sprintf('n = %d',selectedN(k)));
  plot(ax1, x, occupancyRatio(row,:), ':', 'Color', colors(k,:), ...
    'LineWidth',1.1, 'HandleVisibility','off');
end
yline(ax1,1,'k--','HandleVisibility','off');
xlabel(ax1,'Hit rate (%)');
ylabel(ax1,'K_{\theta} / K_0');
title(ax1,{'Relative measured gain with true weights = 1', ...
  'solid: exact WTA; dotted: simple p_{max} approximation'});
legend(ax1,'Location','northeast');
grid(ax1,'on'); box(ax1,'off');
xlim(ax1,[50 100]);

% B: absolute local gains.
ax2 = nexttile(tl,2);
hold(ax2,'on');
for k = 1:numel(selectedN)
  row = find(nValues == selectedN(k), 1);
  plot(ax2, x, K0(row,:), '-', 'Color', colors(k,:), 'LineWidth',1.8, ...
    'DisplayName',sprintf('K_0, n=%d',selectedN(k)));
  plot(ax2, x, Ktheta(row,:), '--', 'Color', colors(k,:), 'LineWidth',1.4, ...
    'DisplayName',sprintf('K_{\theta}, n=%d',selectedN(k)));
end
xlabel(ax2,'Hit rate (%)');
ylabel(ax2,'Local probability gain per evidence unit');
title(ax2,'Absolute measured gains');
grid(ax2,'on'); box(ax2,'off');
xlim(ax2,[50 100]);
legend(ax2,'Location','eastoutside');

% C: all requested n values.
ax3 = nexttile(tl,3);
imagesc(ax3, x, nValues, kernelRatio);
set(ax3,'YDir','normal');
xlabel(ax3,'Hit rate (%)');
ylabel(ax3,'Mechanisms per RDK, n');
title(ax3,'K_{\theta} / K_0 for all n');
cb = colorbar(ax3);
cb.Label.String = 'Measured gain ratio';
clim(ax3,[0 1]);

% D: mean advantage required to create the hit rate.
ax4 = nexttile(tl,4);
hold(ax4,'on');
for k = 1:numel(selectedN)
  row = find(nValues == selectedN(k), 1);
  plot(ax4, x, mu0(row,:)./sigma, '-o', 'Color', colors(k,:), ...
    'LineWidth',1.8, 'MarkerSize',4, 'DisplayName',sprintf('n = %d',selectedN(k)));
end
xlabel(ax4,'Hit rate (%)');
ylabel(ax4,'Drift advantage / internal-noise SD');
title(ax4,'Advantage used to set performance');
grid(ax4,'on'); box(ax4,'off');
xlim(ax4,[50 100]);
legend(ax4,'Location','northwest');

title(tl, {'Hard-max selection changes measured gain despite equal underlying weights', ...
  sprintf('Two RDKs, n mechanisms each; 100%% point evaluated at %.1f%%',100*opts.NearPerfectRate)}, ...
  'FontWeight','bold');

if opts.SavePDF
  pdfPath = char(opts.PDFPath);
  pdfFolder = fileparts(pdfPath);
  if ~isempty(pdfFolder) && ~isfolder(pdfFolder)
    mkdir(pdfFolder);
  end
  exportgraphics(fig, pdfPath, 'ContentType','vector');
end

results = struct();
results.version = 1;
results.description = ['Two-pool Gaussian hard-max model with all underlying ' ...
  'weights fixed at one. K values are local derivatives of hit probability.'];
results.options = opts;
results.nValues = nValues;
results.requestedHitRates = requestedHitRates;
results.modelHitRates = modelHitRates;
results.realizedHitRates = realizedHitRate;
results.driftAdvantage = mu0;
results.K0 = K0;
results.Ktheta = Ktheta;
results.kernelRatio = kernelRatio;
results.occupancyApproximation = occupancyRatio;
results.figure = fig;
end

% -------------------------------------------------------------------------
function mu = calibrateDriftAdvantage(n, targetHit, sigma)
if abs(targetHit - 0.5) < 1e-12
  mu = 0;
  return;
end

objective = @(x) hitProbability(n, x, 0, sigma) - targetHit;
upper = 6*sigma;
while objective(upper) < 0
  upper = 2*upper;
  if upper > 100*sigma
    error('Could not bracket drift advantage for hit rate %.6f.',targetHit);
  end
end
mu = fzero(objective, [0 upper]);
end

% -------------------------------------------------------------------------
function pHit = hitProbability(n, mu0, muTheta, sigma)
% One change mechanism has mean mu0, one has mean muTheta, and the other
% n-2 change mechanisms have mean zero. All n no-change mechanisms have
% mean zero. Mechanisms are mutually independent Gaussian variables.
integrand = @(x) changeMaximumPDF(x,n,mu0,muTheta,sigma) .* ...
  normalCDF(x./sigma).^n;
pHit = integral(integrand, -Inf, Inf, 'AbsTol',1e-10, 'RelTol',1e-8);
pHit = min(max(pHit,0),1);
end

% -------------------------------------------------------------------------
function f = changeMaximumPDF(x,n,mu0,muTheta,sigma)
F0 = normalCDF(x./sigma);
Fdrift = normalCDF((x-mu0)./sigma);
Ftheta = normalCDF((x-muTheta)./sigma);
f0 = normalPDF(x./sigma)./sigma;
fdrift = normalPDF((x-mu0)./sigma)./sigma;
ftheta = normalPDF((x-muTheta)./sigma)./sigma;

basePower = F0.^(n-2);
f = fdrift .* Ftheta .* basePower + ...
    Fdrift .* ftheta .* basePower;
if n > 2
  f = f + (n-2) .* Fdrift .* Ftheta .* F0.^(n-3) .* f0;
end
end

% -------------------------------------------------------------------------
function y = normalCDF(x)
y = 0.5 .* erfc(-x ./ sqrt(2));
end

% -------------------------------------------------------------------------
function y = normalPDF(x)
y = exp(-0.5 .* x.^2) ./ sqrt(2*pi);
end
