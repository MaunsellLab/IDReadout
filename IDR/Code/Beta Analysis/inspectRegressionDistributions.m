function T = inspectRegressionDistributions()
% inspectRegressionDistributions
% Collect and inspect INC betaPref, betaProbe, and betaRatio across sessions.

dataRoot = fullfile(domainFolder(mfilename('fullpath')), 'Data');
files = dir(fullfile(dataRoot, 'Probe*', 'Regression', ...
  '*_scalarNoiseRegression.mat'));

rows = repmat(struct( ...
  'fileName', '', ...
  'probeOffsetDeg', NaN, ...
  'betaPref', NaN, ...
  'betaProbe', NaN, ...
  'betaRatio', NaN, ...
  'betaPrefSE', NaN, ...
  'betaProbeSE', NaN, ...
  'betaRatioSE', NaN, ...
  'nTrials', NaN, ...
  'gradientInfNorm', NaN, ...
  'exitflag', NaN), 0, 1);

for iFile = 1:numel(files)
  filePath = fullfile(files(iFile).folder, files(iFile).name);
  S = load(filePath, 'reg');

  if ~isfield(S, 'reg') || ...
      ~isfield(S.reg, 'analysisName') || ...
      ~strcmp(S.reg.analysisName, 'kernelWeightedProbeRegression') || ...
      ~isfield(S.reg, 'fitByStep') || ...
      ~isfield(S.reg.fitByStep, 'inc')
    continue
  end

  F = S.reg.fitByStep.inc;

  row = rowsTemplate();
  row.fileName = files(iFile).name;

  if isfield(S.reg, 'sessionProbeHeader') && ...
      isfield(S.reg.sessionProbeHeader, 'probeDirDeg')
    row.probeOffsetDeg = double(S.reg.sessionProbeHeader.probeDirDeg);
  end

  row.betaPref = fieldOrNaN(F, 'betaPref');
  row.betaProbe = fieldOrNaN(F, 'betaProbe');
  row.betaRatio = fieldOrNaN(F, 'betaRatio');
  row.betaPrefSE = fieldOrNaN(F, 'betaPrefSE');
  row.betaProbeSE = fieldOrNaN(F, 'betaProbeSE');
  row.betaRatioSE = fieldOrNaN(F, 'betaRatioSE');
  row.nTrials = fieldOrNaN(F, 'nTrials');
  row.gradientInfNorm = fieldOrNaN(F, 'gradientInfNorm');
  row.exitflag = fieldOrNaN(F, 'exitflag');

  rows(end+1,1) = row; %#ok<AGROW>
end

T = struct2table(rows);

validPref = isfinite(T.betaPref);
validProbe = isfinite(T.betaProbe);
validRatio = isfinite(T.betaRatio);

fprintf('\nINC regression distributions\n');
fprintf('  files:              %d\n', height(T));
printSummary('betaPref',  T.betaPref(validPref));
printSummary('betaProbe', T.betaProbe(validProbe));
printSummary('betaRatio', T.betaRatio(validRatio));

fprintf('\nPotentially unstable ratios\n');
unstable = ...
  ~isfinite(T.betaRatio) | ...
  T.betaPref <= 0 | ...
  abs(T.betaPref) < 2 .* T.betaPrefSE | ...
  T.betaRatioSE > 1;

fprintf('  betaPref <= 0:                 %d\n', sum(T.betaPref <= 0, 'omitnan'));
fprintf('  |betaPref| < 2 SE:             %d\n', ...
  sum(abs(T.betaPref) < 2 .* T.betaPrefSE, 'omitnan'));
fprintf('  betaRatio SE > 1:              %d\n', ...
  sum(T.betaRatioSE > 1, 'omitnan'));
fprintf('  any instability criterion:     %d\n', sum(unstable));

disp(T(unstable, {'fileName','probeOffsetDeg','betaPref', ...
  'betaPrefSE','betaProbe','betaProbeSE','betaRatio','betaRatioSE'}));

figure('Color','w','Position',[100 100 1100 360]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile;
histogram(T.betaPref(validPref));
xlabel('\beta_{pref}');
ylabel('Sessions');
title('INC preferred beta');
xline(0, ':');

nexttile;
histogram(T.betaProbe(validProbe));
xlabel('\beta_{probe}');
ylabel('Sessions');
title('INC probe beta');
xline(0, ':');

nexttile;
histogram(T.betaRatio(validRatio));
xlabel('\beta_{probe}/\beta_{pref}');
ylabel('Sessions');
title('INC beta ratio');
xline(0, ':');
xline(1, '--');

end

function row = rowsTemplate()
row = struct( ...
  'fileName', '', ...
  'probeOffsetDeg', NaN, ...
  'betaPref', NaN, ...
  'betaProbe', NaN, ...
  'betaRatio', NaN, ...
  'betaPrefSE', NaN, ...
  'betaProbeSE', NaN, ...
  'betaRatioSE', NaN, ...
  'nTrials', NaN, ...
  'gradientInfNorm', NaN, ...
  'exitflag', NaN);
end

function value = fieldOrNaN(S, fieldName)
if isfield(S, fieldName)
  value = double(S.(fieldName));
else
  value = NaN;
end
end

function printSummary(name, x)
x = x(isfinite(x));

fprintf('  %-10s n=%3d  median=%8.4g  mean=%8.4g  SD=%8.4g', ...
  name, numel(x), median(x), mean(x), std(x));

if ~isempty(x)
  q = prctile(x, [2.5 25 75 97.5]);
  fprintf('  [2.5 25 75 97.5%%]=[%g %g %g %g]\n', q);
else
  fprintf('\n');
end
end