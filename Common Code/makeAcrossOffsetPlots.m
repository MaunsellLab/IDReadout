function makeAcrossOffsetPlots(acrossOffsetSummary, opts)
% Primary plots for the DOG readout models:
%   1) observed scale values by probe offset, with signed and rectified DOG
%      predictions overlaid when fits are available
%   2) fitted readout/template diagnostics for each successful model

emp = acrossOffsetSummary.empirical;
offsets = [emp.probeOffsetDeg];
obsScale = [emp.pooledScale];
ci95 = vertcat(acrossOffsetSummary.bootstrap.offsetBootstrap.boot95);
hasAbsolute = isfield(acrossOffsetSummary, 'absolute') && ...
  isstruct(acrossOffsetSummary.absolute) && ...
  all(isfield(acrossOffsetSummary.absolute, {'offsetsDeg','values','ciLow','ciHigh','medianPreferredGain'}));

% ---- Plot 1: observed and fit scale by offset ----
fig1 = figure(300); 
clf;
set(fig1, 'Color', 'w', 'WindowStyle', 'docked');

if hasAbsolute
  set(fig1, 'Units', 'inches', 'Position', [1 1 13 6]);
  tiledlayout(fig1, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
  nexttile;
end
hold on;
hObs = errorbar([0, offsets], [1, obsScale], [0, obsScale - ci95(:,1)'], [0, ci95(:,2)' - obsScale], ...
  'ko', 'LineWidth', 1.2, 'MarkerFaceColor', 'k');
plot([0, 180], [0,0], 'k:');

% Typical MT Gaussian direction tuning, normalized to 1 at 0 deg
mtSigmaDeg = 37.5;
mtOffsetsDeg = linspace(0, 180, 361);
mtTuning = exp(-0.5 * (mtOffsetsDeg ./ mtSigmaDeg).^2);
hMT = plot(mtOffsetsDeg, mtTuning, 'k:', 'LineWidth', 1.0, 'DisplayName', 'Typical MT tuning');
legendHandles = [hObs, hMT];
legendLabels = {sprintf('Observed %s Scale (95% CI)', acrossOffsetSummary.type), 
        'Typical MT tuning (\sigma = 37.5 deg)'};
signedRM = acrossOffsetSummary.readoutModels.signedDOG;
rectRM   = acrossOffsetSummary.readoutModels.rectifiedDOG;
hasSignedFit = isfield(signedRM, 'fit') && ~isempty(signedRM.fit) && ...
  isfield(signedRM.fit, 'fitSuccess') && signedRM.fit.fitSuccess;
if hasSignedFit
  hSigned = plot(signedRM.plotOffsetsDeg, signedRM.plotPredictedScale, 'm-', 'LineWidth', 1.2);
  % DOGFitText(0.35, 0.98, 'Signed DOG', signedRM);
  legendHandles(end+1) = hSigned;
  legendLabels{end+1} = 'Fitted Signed DOG';
end
hasRectFit = isfield(rectRM, 'fit') && ~isempty(rectRM.fit) && ...
  isfield(rectRM.fit, 'fitSuccess') && rectRM.fit.fitSuccess;
if hasRectFit
  hRect = plot(rectRM.plotOffsetsDeg, rectRM.plotPredictedScale, 'b-', 'LineWidth', 1.2);
  DOGFitText(0.98, 0.02, 'Rectified DOG', rectRM);
  legendHandles(end+1) = hRect;
  legendLabels{end+1} = 'Fitted Rectified DOG';
end
legend(legendHandles, legendLabels, 'Location', 'southwest');
if hasSignedFit || hasRectFit
  title(sprintf('Normalized %s scales', acrossOffsetSummary.type));
else
  title(sprintf('Normalized %s scales (no fit)', acrossOffsetSummary.type));
end
scaleText(0.98, 0.98, offsets, obsScale, ci95, emp);
xlabel('Probe Offset (deg)');
ylabel(sprintf('Normalized %s Scale', acrossOffsetSummary.type));
xlim([0, 180]);
ylim([-0.5, 1.25]);
box off;

if hasAbsolute
  plotAbsoluteGain(acrossOffsetSummary, opts, mtOffsetsDeg, mtTuning, ...
    hasSignedFit, hasRectFit);
  sgtitle(sprintf('%s -- normalized and absolute %s gains (%d bootstraps)', ...
    opts.Animal, acrossOffsetSummary.type, opts.NBoot));
elseif hasSignedFit || hasRectFit
  title(sprintf('%s -- DOG Fits to Normalized %s Scales (%d bootstraps)', ...
    opts.Animal, acrossOffsetSummary.type, opts.NBoot));
else
  title(sprintf('%s -- Normalized %s Scales (No Fit Over %d bootstraps)', ...
    opts.Animal, acrossOffsetSummary.type, opts.NBoot));
end
saveas(fig1, fullfile(opts.PlotDir, sprintf('%sScaleFits_%s.pdf', acrossOffsetSummary.type, opts.Animal)));

% ---- Plots 2/3: fitted readout over MT preferred direction ----
if hasSignedFit
  tmpSummary = acrossOffsetSummary;
  tmpSummary.readoutModel = signedRM;
  plotReadoutDiagnostics(301, tmpSummary, opts);
end
if hasRectFit
  tmpSummary = acrossOffsetSummary;
  tmpSummary.readoutModel = rectRM;
  plotReadoutDiagnostics(302, tmpSummary, opts);
end
end

% ========================================================================
function plotAbsoluteGain(acrossOffsetSummary, ~, mtOffsetsDeg, mtTuning, hasSignedFit, hasRectFit)
% Absolute gains with normalized reference curves scaled by median gPref.

A = acrossOffsetSummary.absolute;
x = A.offsetsDeg(:)';
y = A.values(:)';
ciLow = A.ciLow(:)';
ciHigh = A.ciHigh(:)';
gPref = A.medianPreferredGain;

nexttile; hold on;
hObs = errorbar(x, y, y-ciLow, ciHigh-y, 'ko', 'LineWidth', 1.2, ...
  'MarkerFaceColor', 'k');
plot([0, 180], [0, 0], 'k:');

hMT = plot(mtOffsetsDeg, gPref .* mtTuning, 'k:', 'LineWidth', 1.0);
allY = [0, ciLow, ciHigh, gPref .* mtTuning];
legendHandles = [hObs, hMT];
legendLabels = {sprintf('Observed absolute %s gain (95%% CI)', acrossOffsetSummary.type), ...
  sprintf('MT tuning x median gPref (%.3g)', gPref)};

signedRM = acrossOffsetSummary.readoutModels.signedDOG;
rectRM = acrossOffsetSummary.readoutModels.rectifiedDOG;
if hasSignedFit
  signedAbs = gPref .* signedRM.plotPredictedScale;
  hSigned = plot(signedRM.plotOffsetsDeg, signedAbs, ...
    'm-', 'LineWidth', 1.2);
  allY = [allY, signedAbs]; %#ok<AGROW>
  legendHandles(end+1) = hSigned;
  legendLabels{end+1} = 'Normalized Signed DOG x median gPref';
end
if hasRectFit
  rectAbs = gPref .* rectRM.plotPredictedScale;
  hRect = plot(rectRM.plotOffsetsDeg, rectAbs, ...
    'b-', 'LineWidth', 1.2);
  allY = [allY, rectAbs]; %#ok<AGROW>
  legendHandles(end+1) = hRect;
  legendLabels{end+1} = 'Normalized Rectified DOG x median gPref';
end

legend(legendHandles, legendLabels, 'Location', 'southwest');
title(sprintf('Absolute %s gains', acrossOffsetSummary.type));
absoluteGainText(0.98, 0.98, x, y, ciLow, ciHigh);
xlabel('Probe Offset (deg)');
ylabel(sprintf('Absolute %s Gain', acrossOffsetSummary.type));
xlim([0, 180]);
ylim([-0.50, 1.25]);
% allY = allY(isfinite(allY));
% yl = [min(allY), max(allY)];
% if all(isfinite(yl)) && yl(2) > yl(1)
%   pad = 0.08 * (yl(2) - yl(1));
%   ylim([yl(1)-pad, yl(2)+pad]);
% end
box off;
end

% ========================================================================
function absoluteGainText(xText, yText, offsets, gains, ciLow, ciHigh)
lines = cell(1, numel(offsets));
for i = 1:numel(offsets)
  if abs(offsets(i)) < 1e-9
    label = 'preferred';
  else
    label = sprintf('%g deg', offsets(i));
  end
  lines{i} = sprintf('%s: gain %.2f, %.2f-%.2f 95%% CI', ...
    label, gains(i), ciLow(i), ciHigh(i));
end
text(xText, yText, strjoin(lines, newline), 'Units', 'normalized', ...
  'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 8, ...
  'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 4, ...
  'Interpreter', 'none');
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

% ========================================================================
function scaleText(x, y, offsets, obsScale, ci95, emp)
% One-model parameter/goodness-of-fit text block.

C = arrayfun(@(x) valueOrNaN(x.nTrialsByStep), emp, 'UniformOutput', false);
nTrialsByStep = vertcat(C{:});
stepTypes = [emp.stepType];
lines = {};
for index = 1:numel(offsets)
  if isnan(obsScale(index))
    continue;
  end
  nTrials = nTrialsByStep(index, stepTypes(index));
  lines{end+1} = sprintf('%3d°: scale %.2f, %.2f-%.2f 95%% CI (n = %5d)', ...
    offsets(index), obsScale(index), ci95(index, 1),  ci95(index, 2), nTrials); %#ok<AGROW>
end
txt = strjoin(lines, newline);
text(x, y, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 8, ...
  'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 4, 'Interpreter', 'none');
end

% ========================================================================
function y = valueOrNaN(x)
if isempty(x)
  y = [NaN NaN];
else
  y = double(x);
end
end

% ========================================================================
function plotReadoutDiagnostics(figNum, acrossOffsetSummary, opts)
% Plot fitted readout, MT templates, and their products to visualize how
% overlap determines predicted normalized scale.

  rm = acrossOffsetSummary.readoutModel;
  if isfield(rm, 'templateMode')
      templateMode = rm.templateMode;
      templateMode(1) = upper(templateMode(1));
  else
      templateMode = 'Signed';
  end
  if ~isfield(rm, 'fit') || isempty(rm.fit) || ~rm.fit.fitSuccess
      return;
  end

  phiDeg = rm.phiDeg(:)';
  aPhi   = rm.readoutPhi(:)';   % normalized display readout, a(0)=1
  mtp = rm.mtForwardModelParams;
  mtModel = makeMTReadoutForwardModel('sigmaMTDeg', mtp.sigmaMTDeg, 'phiDeg', mtp.phiDeg);
  offsetsDeg = [0, rm.fit.offsetsDeg];  
  nOffsets = numel(offsetsDeg);
  
  deltaM   = cell(1, nOffsets);
  prodTerm = cell(1, nOffsets);
  overlap  = nan(1, nOffsets);
  posPart  = nan(1, nOffsets);
  negPart  = nan(1, nOffsets);

  for i = 1:nOffsets
      deltaM{i} = mtReadoutTemplate(offsetsDeg(i), mtModel, 'TemplateMode', templateMode);
      prodTerm{i} = aPhi .* deltaM{i};
      overlap(i) = sum(prodTerm{i});
      posPart(i) = sum(max(prodTerm{i}, 0));
      negPart(i) = sum(min(prodTerm{i}, 0));
  end
  
  idx0 = find(abs(offsetsDeg) < 1e-9, 1, 'first');
  if isempty(idx0)
      warning('plotReadoutDiagnostics: OffsetsDeg does not include 0. Ratios will not be shown.');
  end
  
  % ---- Build fit-vs-flat comparison figure ----
  prodTermFit  = cell(1, nOffsets);
  prodTermFlat = cell(1, nOffsets);
  for i = 1:nOffsets
      prodTermFit{i}  = aPhi .* deltaM{i};
      prodTermFlat{i} = ones(size(aPhi)) .* deltaM{i};   % flat readout = 1
  end

  % -- set up figure to plot three panels
  fig = figure(figNum);
  clf(fig);
  tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
  lineCol = lines(nOffsets);
  
  % ---- top panel: readout function ----
  nexttile; hold on;
  hFitReadout  = plot(phiDeg, aPhi, 'k-', 'LineWidth', 2);
  plot(phiDeg, zeros(size(phiDeg)), 'k:', 'LineWidth', 1);
  xlabel('\phi (deg)');
  ylabel('a(\phi)');
  title(sprintf('%s -- %s Fitted DOG readout (%s template, %d bootstraps)', acrossOffsetSummary.type, ...
        opts.Animal, templateMode, opts.NBoot));
  legend(hFitReadout, {'Fitted readout a(\phi)'}, 'Location', 'northeast');
  paramText = cell(rm.nFreeParams, 1);
  for p = 1:rm.nFreeParams
    paramText{p} = sprintf('%s: %.4f', rm.paramNames{p}, rm.params(p));
  end
  text(-100, 0.95, paramText, 'horizontalAlignment', 'right', 'VerticalAlignment', 'top');
  ylimits = ylim();
  ylim([min(0.2, ylimits(1)), max(1.1, ylimits(2))]);
  box off;
  
  % ---- middle panel: MT populations responses ----
  nexttile; hold on;
  title(sprintf('MT Population Responses to Probes (%s template; Flat Readout)', templateMode));
  % MT templates, same colors used in the lower panel
  hTemplates = gobjects(1, nOffsets);
  for i = 1:nOffsets
      hTemplates(i) = plot(phiDeg, deltaM{i}, '-', 'Color', lineCol(i,:), 'LineWidth', 1.5);
  end
  yline(0, 'k:');
  xlabel('\phi (deg)');
  ylabel('\Delta m(\phi;\delta)');
  legend(hTemplates, [arrayfun(@(d) sprintf('\\Delta m(\\phi;%g^\\circ)', d), offsetsDeg, ...
       'UniformOutput', false)], 'Location', 'best');
  box off;
  
  % ---- Bottom panel: overlap contribution functions ----
  nexttile; hold on;
  title('Weighted Population Responses (Fit)');
  legendHandles = gobjects(0);
  legendLabels = {};
  
  for i = 1:nOffsets
      h1 = plot(phiDeg, prodTermFit{i}, '-', 'Color', lineCol(i,:), 'LineWidth', 2);
      legendHandles(end+1) = h1; %#ok<AGROW>
      legendLabels{end+1} = sprintf('%g^\\circ fit, <a,\\Deltam> = %.2g; S_{fit} %.2f', ...
        offsetsDeg(i), overlap(i), overlap(i)/overlap(1)); %#ok<AGROW>
  end
  yline(0, 'k:');
  xlabel('\phi (deg)');
  ylabel('a(\phi)\Delta m(\phi;\delta)');
  legend(legendHandles, legendLabels, 'Location', 'best');
  box off;

  saveas(fig, fullfile(opts.PlotDir, sprintf('%sReadout_%s_%s.pdf', ...
              acrossOffsetSummary.type, templateMode, opts.Animal)));
end
