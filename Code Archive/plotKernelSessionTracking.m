function plotKernelSessionTracking(varargin)
% plotKernelSessionTracking
%
% Plot session tracking summaries for one selected sideType / stepType.
% Uses session index on the x-axis, with sparse date labels.

P = inputParser;
addParameter(P, 'summaryDir', '/Users/Shared/Data/IDReadout/Data/KernelSummaries', @ischar);
addParameter(P, 'plotCI', 68, @(x) isnumeric(x) && isscalar(x));
addParameter(P, 'savePlot', false, @islogical);
addParameter(P, 'plotFile', '', @ischar);
addParameter(P, 'tickStep', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parse(P, varargin{:});
R = P.Results;

files = dir(fullfile(R.summaryDir, '*_kernelSummary.mat'));
if isempty(files)
  warning('No summary files found.');
  return
end

S = [];
matFileDir = '/Users/Shared/Data/IDReadout/Data/Converted';
for i = 1:numel(files)
  baseName = erase(files(i).name, '_kernelSummary.mat');
  h = load(fullfile(matFileDir, [baseName, '.mat']), 'header');
  % [~, baseName, ~] = fileparts(erase(files(i).name, '_kernelSummary.mat'));
  % h.fileName = [erase(files(i).name, '_kernelSummary.mat'), '.dat'];
  if ~excludeFile(h.header)
    tmp = load(fullfile(files(i).folder, files(i).name), 'summary');
    S = [S; tmp.summary];
  end
end

dates = [S.date]';
scale = arrayfun(@(x) x.scale.estimate, S)';
prefEnergy = arrayfun(@(x) x.pref.energy, S)';
nTrials = arrayfun(@(x) x.metrics.nTrials, S)';
pCorrect = arrayfun(@(x) x.metrics.pCorrect, S)';

if R.plotCI == 68
  ciLo = arrayfun(@(x) x.scale.ci68(1), S)';
  ciHi = arrayfun(@(x) x.scale.ci68(2), S)';
  ciWidth = arrayfun(@(x) x.scale.ci68Width, S)';
else
  ciLo = arrayfun(@(x) x.scale.ci95(1), S)';
  ciHi = arrayfun(@(x) x.scale.ci95(2), S)';
  ciWidth = arrayfun(@(x) x.scale.ci95Width, S)';
end

% Keep only entries with valid dates
keep = ~isnat(dates);
dates = dates(keep);
scale = scale(keep);
prefEnergy = prefEnergy(keep);
nTrials = nTrials(keep);
pCorrect = pCorrect(keep);
ciLo = ciLo(keep);
ciHi = ciHi(keep);
ciWidth = ciWidth(keep);
keep = ~arrayfun(@(x) x.flags.excluded, S);
S = S(keep);

% Sort by date
[dates, ord] = sort(dates);
scale = scale(ord);
prefEnergy = prefEnergy(ord);
nTrials = nTrials(ord);
pCorrect = pCorrect(ord);
ciLo = ciLo(ord);
ciHi = ciHi(ord);
ciWidth = ciWidth(ord);
S = S(ord);

% Session index axis
x = 1:numel(scale);

% Sparse date ticks
tickStep = max(1, round(R.tickStep));
tickIdx = unique([1, tickStep:tickStep:numel(x), numel(x)]);
tickLabels = cell(size(tickIdx));
for i = 1:numel(tickIdx)
  tickLabels{i} = datestr(dates(tickIdx(i)), 'yyyy-mm-dd');
end

figure;

subplot(5,1,1); hold on
for i = 1:numel(scale)
  if isfinite(ciLo(i)) && isfinite(ciHi(i))
    plot([x(i) x(i)], [ciLo(i) ciHi(i)], 'k-');
  end
end
plot(x, scale, 'ko-');
ylabel('Scale');
xlim([1 max(1, numel(x))]);

subplot(5,1,2);
plot(x, prefEnergy, 'ko-');
ylabel('Pref energy');
xlim([1 max(1, numel(x))]);

subplot(5,1,3);
plot(x, nTrials, 'ko-');
ylabel('Trials');
xlim([1 max(1, numel(x))]);

subplot(5,1,4);
plot(x, pCorrect, 'ko-');
ylabel('% correct');
xlim([1 max(1, numel(x))]);

subplot(5,1,5); hold on
validCI = isfinite(ciWidth);
plot(x(validCI), ciWidth(validCI), 'ko-');
ylabel(sprintf('%d%% CI width', R.plotCI));
xlabel('Session');
xlim([1 max(1, numel(x))]);

% Apply shared x ticks to all subplots
ax = findall(gcf, 'Type', 'axes');
for i = 1:numel(ax)
  set(ax(i), 'XTick', tickIdx, 'XTickLabel', tickLabels);
end
xtickangle(45)

if ~isempty(S)
  sgtitle(sprintf('Session tracking: %s, %s', ...
    S(1).track.sideLabel, S(1).track.stepLabel));
end

if R.savePlot && ~isempty(R.plotFile)
  saveas(gcf, R.plotFile);
end

x = [];
y = [];

for i = 1:numel(S)
    if ~S(i).scale.valid
        continue
    end
    if ~isfinite(S(i).pref.energy) || ~isfinite(S(i).scale.estimate)
        continue
    end
    x(end+1) = S(i).pref.energy;
    y(end+1) = S(i).scale.estimate;
end

figure;
plot(x, y, 'ko');
xlabel('Pref energy');
ylabel('Session scale');

[rPearson, pPearson] = corr(x(:), y(:), 'type', 'Pearson');
[rSpearman, pSpearman] = corr(x(:), y(:), 'type', 'Spearman');

title(sprintf('r = %.2f (p=%.3g), rho = %.2f (p=%.3g)', ...
    rPearson, pPearson, rSpearman, pSpearman));

s = scale(:);
e = prefEnergy(:);

ok = isfinite(s) & isfinite(e) & e > 0;

simpleMean = mean(s(ok))
energyWeightedMean = sum(e(ok) .* s(ok)) / sum(e(ok))

end