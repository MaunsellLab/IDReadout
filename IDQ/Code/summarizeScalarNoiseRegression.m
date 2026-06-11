function summaryTable = summarizeScalarNoiseRegression(varargin)
% summarizeScalarNoiseRegression  Summarize scalar noise regressions.
%
%   summaryTable = summarizeScalarNoiseRegression()
%
% Searches Data/Probe*/Regression/*_scalarNoiseRegression.mat and returns a
% table with the main raw and standardized coefficients.

cleanupObj = initProjectPath(); %#ok<NASGU>
p = inputParser;
p.addParameter('RootFolder', folderPath(), @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

rootFolder = char(p.Results.RootFolder);
dataRoot = fullfile(rootFolder, 'Data');
files = dir(fullfile(dataRoot, 'Probe*', 'Regression', '*_scalarNoiseRegression.mat'));

rows = struct([]);

for iFile = 1:numel(files)
    regPath = fullfile(files(iFile).folder, files(iFile).name);
    S = load(regPath, 'reg');
    reg = S.reg;

    r = struct();
    r.regPath = string(regPath);
    r.sessionID = string(getFieldOrDefault(reg.sessionProbeHeader, 'sessionID', ''));
    r.probeTag = string(getFieldOrDefault(reg.sessionProbeHeader, 'probeTag', ''));
    r.probeDirDeg = getFieldOrDefault(reg.sessionProbeHeader, 'probeDirDeg', NaN);
    r.nTrials = reg.nTrials;
    r.meanCorrect = reg.meanCorrect;

    r.beta_change_raw = reg.effectSummary.bothPref.DchangePref.beta;
    r.se_change_raw = reg.effectSummary.bothPref.DchangePref.se;
    r.p_change_raw = reg.effectSummary.bothPref.DchangePref.p;

    r.beta_noChange_raw = reg.effectSummary.bothPref.DnoChangePref.beta;
    r.se_noChange_raw = reg.effectSummary.bothPref.DnoChangePref.se;
    r.p_noChange_raw = reg.effectSummary.bothPref.DnoChangePref.p;

    r.beta_change_z = reg.effectSummary.zBothPref.ZchangePref.beta;
    r.se_change_z = reg.effectSummary.zBothPref.ZchangePref.se;
    r.p_change_z = reg.effectSummary.zBothPref.ZchangePref.p;

    r.beta_noChange_z = reg.effectSummary.zBothPref.ZnoChangePref.beta;
    r.se_noChange_z = reg.effectSummary.zBothPref.ZnoChangePref.se;
    r.p_noChange_z = reg.effectSummary.zBothPref.ZnoChangePref.p;

    r.beta_diff_z = reg.effectSummary.zDiffPref.ZdiffPref.beta;
    r.se_diff_z = reg.effectSummary.zDiffPref.ZdiffPref.se;
    r.p_diff_z = reg.effectSummary.zDiffPref.ZdiffPref.p;

    r.beta_changeMinusNoChange_z = r.beta_change_z - r.beta_noChange_z;

    r.beta_changeProbe_raw = reg.effectSummary.changePrefProbe_DchangeProbe.beta;
    r.se_changeProbe_raw   = reg.effectSummary.changePrefProbe_DchangeProbe.se;
    r.p_changeProbe_raw    = reg.effectSummary.changePrefProbe_DchangeProbe.p;

    r.beta_changePref_raw_joint = reg.effectSummary.changePrefProbe_DchangePref.beta;
    r.se_changePref_raw_joint   = reg.effectSummary.changePrefProbe_DchangePref.se;
    r.p_changePref_raw_joint    = reg.effectSummary.changePrefProbe_DchangePref.p;

    r.beta_changeProbe_z = reg.effectSummary.zChangePrefProbe_ZchangeProbe.beta;
    r.se_changeProbe_z   = reg.effectSummary.zChangePrefProbe_ZchangeProbe.se;
    r.p_changeProbe_z    = reg.effectSummary.zChangePrefProbe_ZchangeProbe.p;

    r.beta_changePref_z_joint = reg.effectSummary.zChangePrefProbe_ZchangePref.beta;
    r.se_changePref_z_joint   = reg.effectSummary.zChangePrefProbe_ZchangePref.se;
    r.p_changePref_z_joint    = reg.effectSummary.zChangePrefProbe_ZchangePref.p;

    r.ratio_changeProbeOverPref_raw = reg.coeffRatios.changeProbeOverPref_raw;
    r.ratio_changeProbeOverPref_z   = reg.coeffRatios.changeProbeOverPref_z;

    if isempty(rows)
        rows = r;
    else
        rows(end+1) = r; %#ok<AGROW>
    end
end

minTrials = 100;
if isempty(rows)
    summaryTable = table();
    T = table();
else
    summaryTable = struct2table(rows);
    T = summaryTable(summaryTable.nTrials >= minTrials, :);
end

% Probe-level summaries.
summaryVars = {'beta_change_z', 'beta_noChange_z', 'beta_diff_z', ...
    'beta_changeMinusNoChange_z', ...
    'beta_changePref_raw_joint', 'beta_changeProbe_raw', ...
    'beta_changePref_z_joint', 'beta_changeProbe_z', ...
    'ratio_changeProbeOverPref_raw', 'ratio_changeProbeOverPref_z', ...
    'nTrials', 'meanCorrect'};
[Gidx, probeTags] = findgroups(T.probeTag);

probeSummary = table();
probeSummary.probeTag = probeTags;
probeSummary.nSessions = splitapply(@numel, T.nTrials, Gidx);
for iVar = 1:numel(summaryVars)
    v = summaryVars{iVar};
    x = T.(v);
    probeSummary.(['mean_' v]) = splitapply(@localMeanFinite, x, Gidx);
    probeSummary.(['std_'  v]) = splitapply(@localStdFinite,  x, Gidx);
    probeSummary.(['sem_'  v]) = splitapply(@localSemFinite,  x, Gidx);
end
probeSummary.ratioOfMeans_changeProbeOverPref_raw = ...
  probeSummary.mean_beta_changeProbe_raw ./ ...
  probeSummary.mean_beta_changePref_raw_joint;
probeSummary.ratioOfMeans_changeProbeOverPref_z = ...
  probeSummary.mean_beta_changeProbe_z ./ ...
  probeSummary.mean_beta_changePref_z_joint;

ratioMean = nan(height(probeSummary), 1);
ratioCI = nan(height(probeSummary), 2);
for iProbe = 1:height(probeSummary)
  idx = T.probeTag == probeTags(iProbe);
  xProbe = T.beta_changeProbe_raw(idx);
  xPref  = T.beta_changePref_raw_joint(idx);
  [ratioMean(iProbe), ratioCI(iProbe, :)] = bootstrapRatioOfMeans(xProbe, xPref, 10000);
end
probeSummary.ratioOfMeans_changeProbeOverPref_raw_boot = ratioMean;
probeSummary.ratioOfMeans_changeProbeOverPref_raw_ciLow = ratioCI(:, 1);
probeSummary.ratioOfMeans_changeProbeOverPref_raw_ciHigh = ratioCI(:, 2);

disp(probeSummary);

summary = struct();
summary.minTrials = minTrials;
summary.bySession = summaryTable;
summary.byProbe = probeSummary;

outPath = fullfile(dataRoot, 'scalarNoiseRegression_summaryTable.mat');
save(outPath, 'T', 'summary');

% Probe-level paired tests for change versus noChange coefficients.
probeTags = unique(T.probeTag);

for iProbe = 1:numel(probeTags)
  thisProbe = probeTags(iProbe);
  idx = T.probeTag == thisProbe;

  d = T.beta_changeMinusNoChange_z(idx);
  d = d(isfinite(d));

  if numel(d) >= 3
    [~, pT, ~, statsT] = ttest(d, 0);
    pSR = signrank(d, 0);
  else
    pT = NaN;
    pSR = NaN;
    statsT = struct('tstat', NaN, 'df', NaN, 'sd', NaN);
  end

  fprintf('%s: n=%d, mean diff=%.4f, SEM=%.4f, t=%.3f, p_t=%.4g, p_signrank=%.4g\n', ...
    thisProbe, numel(d), mean(d), std(d)/sqrt(numel(d)), statsT.tstat, pT, pSR);
end

fprintf('Loaded %d scalar noise regression files.\n', height(summaryTable));
fprintf('Saved summary table: %s\n', outPath);
disp(summaryTable);
plotScalarNoiseCoeffByOffset(T);
plotProbeOverPrefRatioByOffset(T);
plotProbeOverPrefRatioOfMeansByOffset(probeSummary);
end

%% plotScalarNoiseCoeffByOffset
function fig = plotScalarNoiseCoeffByOffset(summaryTable, varargin)
% plotScalarNoiseCoeffByOffset
%
% One summary plot:
%   x-axis: probe offset in degrees
%   y-axis: standardized logistic coefficient
%
% Two series:
%   beta_change_z
%   beta_noChange_z
%
% Error bars:
%   SEM by default, or approximate 95% CI.

p = inputParser;
p.addParameter('MinTrials', 100, @isscalar);
p.addParameter('ErrorType', 'sem', @(x) ismember(lower(string(x)), ["sem", "ci95"]));
p.addParameter('UseAbsOffset', false, @islogical);
p.parse(varargin{:});

minTrials = p.Results.MinTrials;
errorType = lower(string(p.Results.ErrorType));

T = summaryTable(summaryTable.nTrials >= minTrials, :);

if p.Results.UseAbsOffset
    T.plotProbeDirDeg = abs(T.probeDirDeg);
else
    T.plotProbeDirDeg = T.probeDirDeg;
end

probeDirs = unique(T.plotProbeDirDeg);
probeDirs = sort(probeDirs(:));

nProbe = numel(probeDirs);

meanChange = nan(nProbe, 1);
meanNoChange = nan(nProbe, 1);
errChange = nan(nProbe, 1);
errNoChange = nan(nProbe, 1);
nSessions = nan(nProbe, 1);

for iProbe = 1:nProbe

    thisDir = probeDirs(iProbe);
    idx = T.plotProbeDirDeg == thisDir;

    xChange = T.beta_change_z(idx);
    xNoChange = T.beta_noChange_z(idx);

    xChange = xChange(isfinite(xChange));
    xNoChange = xNoChange(isfinite(xNoChange));

    nSessions(iProbe) = numel(xChange);

    meanChange(iProbe) = mean(xChange);
    meanNoChange(iProbe) = mean(xNoChange);

    semChange = std(xChange) ./ sqrt(numel(xChange));
    semNoChange = std(xNoChange) ./ sqrt(numel(xNoChange));

    switch errorType
        case "sem"
            errChange(iProbe) = semChange;
            errNoChange(iProbe) = semNoChange;

        case "ci95"
            % Approximate t-based 95% CI half-width.
            if numel(xChange) > 1
                errChange(iProbe) = tinv(0.975, numel(xChange)-1) .* semChange;
            end

            if numel(xNoChange) > 1
                errNoChange(iProbe) = tinv(0.975, numel(xNoChange)-1) .* semNoChange;
            end
    end
end

fig = figure('Color', 'w', 'Position', [100 100 750 450]);
hold on;

errorbar(probeDirs, meanChange, errChange, '-o', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto', ...
    'DisplayName', 'Change side');

errorbar(probeDirs, meanNoChange, errNoChange, '-o', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto', ...
    'DisplayName', 'No-change side');

yline(0, ':');

xlabel('Probe offset (deg)');
ylabel('Standardized logistic coefficient');
title(sprintf('Preferred-noise coefficients by probe offset, n_{min} = %d', minTrials));

legend('Location', 'best');
box off;

% Add n labels near the bottom of the plot.
yl = ylim;
yText = yl(1) + 0.05 * range(yl);

for iProbe = 1:nProbe
    text(probeDirs(iProbe), yText, sprintf('n=%d', nSessions(iProbe)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 8);
end
end

%% =========================================================================
function fig = plotProbeOverPrefRatioByOffset(summaryTable, varargin)

p = inputParser;
p.addParameter('MinTrials', 100, @isscalar);
p.addParameter('ErrorType', 'ci95', @(x) ismember(lower(string(x)), ["sem", "ci95"]));
p.parse(varargin{:});

T = summaryTable(summaryTable.nTrials >= p.Results.MinTrials, :);

probeDirs = sort(unique(T.probeDirDeg));
nProbe = numel(probeDirs);

mu = nan(nProbe, 1);
err = nan(nProbe, 1);
n = nan(nProbe, 1);

for i = 1:nProbe
  idx = T.probeDirDeg == probeDirs(i);
  x = T.ratio_changeProbeOverPref_raw(idx);
  x = x(isfinite(x));

  n(i) = numel(x);
  mu(i) = mean(x);

  sem = std(x) ./ sqrt(numel(x));

  if lower(string(p.Results.ErrorType)) == "ci95" && numel(x) > 1
    err(i) = tinv(0.975, numel(x)-1) .* sem;
  else
    err(i) = sem;
  end
end

fig = figure('Color', 'w', 'Position', [100 100 700 450]);
hold on;

errorbar(probeDirs, mu, err, '-o', ...
  'LineWidth', 1.5, ...
  'MarkerFaceColor', 'auto');

yline(0, ':');
yline(1, ':');

xlabel('Probe offset (deg)');
ylabel('Probe / preferred coefficient ratio');
title(sprintf('Scalar-regression probe/pref ratio, n_{min} = %d', p.Results.MinTrials));
box off;

yl = ylim;
yText = yl(1) + 0.05 * range(yl);

for i = 1:nProbe
  text(probeDirs(i), yText, sprintf('n=%d', n(i)), ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 8);
end
end

%% =========================================================================
function fig = plotProbeOverPrefRatioOfMeansByOffset(probeSummary)

if ismember('probeDirDeg', probeSummary.Properties.VariableNames)
  x = probeSummary.probeDirDeg;
else
  x = nan(height(probeSummary), 1);
  for i = 1:height(probeSummary)
    x(i) = sscanf(char(probeSummary.probeTag(i)), 'Probe%d');
  end
end

y = probeSummary.ratioOfMeans_changeProbeOverPref_raw_boot;
lo = probeSummary.ratioOfMeans_changeProbeOverPref_raw_ciLow;
hi = probeSummary.ratioOfMeans_changeProbeOverPref_raw_ciHigh;

[x, order] = sort(x);
y = y(order);
lo = lo(order);
hi = hi(order);

errLow = y - lo;
errHigh = hi - y;

fig = figure('Color', 'w', 'Position', [100 100 700 450]);
hold on;

errorbar(x, y, errLow, errHigh, '-o', ...
  'LineWidth', 1.5, ...
  'MarkerFaceColor', 'auto');

yline(0, ':');
yline(1, ':');

xlabel('Probe offset (deg)');
ylabel('Probe / preferred coefficient ratio');
title('Scalar-regression ratio of mean raw coefficients');
box off;

end

%% =========================================================================
function [ratioMean, ci95] = bootstrapRatioOfMeans(xProbe, xPref, nBoot)

if nargin < 3
  nBoot = 10000;
end

valid = isfinite(xProbe) & isfinite(xPref);
xProbe = xProbe(valid);
xPref  = xPref(valid);

n = numel(xProbe);

if n == 0
  ratioMean = NaN;
  ci95 = [NaN NaN];
  return;
end

ratioMean = mean(xProbe) ./ mean(xPref);

bootRatio = nan(nBoot, 1);

for b = 1:nBoot
  idx = randi(n, n, 1);
  bootRatio(b) = mean(xProbe(idx)) ./ mean(xPref(idx));
end

ci95 = prctile(bootRatio, [2.5 97.5]);

end

%% =========================================================================
function value = getFieldOrDefault(S, fieldName, defaultValue)

if isstruct(S) && isfield(S, fieldName)
    value = S.(fieldName);
else
    value = defaultValue;
end
end

function m = localMeanFinite(x)
x = x(isfinite(x));
if isempty(x)
  m = NaN;
else
  m = mean(x);
end
end

function s = localStdFinite(x)
x = x(isfinite(x));
if numel(x) < 2
  s = NaN;
else
  s = std(x);
end
end

function se = localSemFinite(x)
x = x(isfinite(x));
if numel(x) < 2
  se = NaN;
else
  se = std(x) ./ sqrt(numel(x));
end
end