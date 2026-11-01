function [summary, fig] = summarizeIDLModelComparison(source, varargin)
% summarizeIDLModelComparison  Provisional session-level IDL model comparison.
%
% [summary, fig] = summarizeIDLModelComparison(result)
% [summary, fig] = summarizeIDLModelComparison('ModelRecovery.mat')
% [summary, fig] = summarizeIDLModelComparison(..., 'PDFPath', 'summary.pdf')
%
% Reads the result saved by fitIDLModelRecoveryAcrossSessions. Within each
% session it compares predictions on the SAME divided-allocation trials.
% The main session-level quantity is 1000 * (LLmax - LLalternative) / nTrials.
% Animals are summarized separately; the combined row treats sessions as
% units, and should not be interpreted as an animal-level test.
% This function does not re-estimate kernels: run the upstream fit within
% each animal when animal-specific training is intended.
%
% The one-sided t tests and two-sided t intervals below are provisional:
% they treat session rates as independent observations. They do not include
% uncertainty from estimating the shared leave-one-session-out kernels or
% correct for inspecting the dataset repeatedly during collection.

p = inputParser;
p.addRequired('source', @(x) isstruct(x) || ischar(x) || isstring(x));
p.addParameter('PDFPath', '', @(x) ischar(x) || isstring(x));
p.addParameter('Visible', 'on', @(x) ischar(x) || isstring(x));
p.parse(source, varargin{:});

if isstruct(source)
  result = source;
else
  S = load(char(source), 'result');
  assert(isfield(S, 'result'), 'The MAT file must contain a result struct.');
  result = S.result;
end

assert(isfield(result, 'sessionTable') && isfield(result, 'heldOutLogLikelihood'), ...
  'Expected the output of fitIDLModelRecoveryAcrossSessions.');
T = result.sessionTable;
required = {'sourceSessionID', 'nHeldOutTrials'};
assert(all(ismember(required, T.Properties.VariableNames)), ...
  'The sessionTable is missing session IDs or held-out trial counts.');
names = string(result.modelNames);
indexL1 = find(names == "L1", 1);
indexL2 = find(names == "L2", 1);
indexMax = find(names == "Linf", 1);
assert(~isempty(indexL1) && ~isempty(indexL2) && ~isempty(indexMax), ...
  'Expected the L1, L2, and Linf models.');
LL = double(result.heldOutLogLikelihood);
n = double(T.nHeldOutTrials(:));
assert(size(LL,1) == height(T) && all(n > 0) && ...
  all(isfinite(LL(:))), 'Invalid session scores or held-out trial counts.');

ids = string(T.sourceSessionID(:));
animal = strings(height(T), 1);
for k = 1:height(T)
  token = regexp(char(ids(k)), '_([^_]+)_[0-9]{8}$', 'tokens', 'once');
  assert(~isempty(token), 'Cannot parse animal from session ID %s.', ids(k));
  animal(k) = string(token{1});
end

dL1 = LL(:,indexMax) - LL(:,indexL1);
dL2 = LL(:,indexMax) - LL(:,indexL2);
rL1 = 1000 * dL1 ./ n;
rL2 = 1000 * dL2 ./ n;
sessionTable = table(ids, animal, n, dL1, dL2, rL1, rL2, ...
  'VariableNames', {'sourceSessionID','animal','nHeldOutTrials', ...
  'deltaLLMaxVsL1','deltaLLMaxVsL2', ...
  'deltaLLPer1000MaxVsL1','deltaLLPer1000MaxVsL2'});

animals = unique(animal, 'stable');
group = strings(0,1);
comparison = strings(0,1);
nSessions = zeros(0,1);
nDividedTrials = zeros(0,1);
nPositive = zeros(0,1);
meanPer1000 = zeros(0,1);
sdPer1000 = zeros(0,1);
lower95 = zeros(0,1);
upper95 = zeros(0,1);
lowerOneSided95 = zeros(0,1);
oneSidedP = zeros(0,1);
totalDeltaLL = zeros(0,1);
pooledPer1000 = zeros(0,1);

groups = animals;
if numel(animals) > 1
  groups(end+1) = "All sessions";
end
for g = 1:numel(groups)
  if groups(g) == "All sessions"
    keep = true(height(T),1);
  else
    keep = animal == groups(g);
  end
  for m = 1:2
    if m == 1
      d = dL2(keep); r = rL2(keep); label = "Max - L2";
    else
      d = dL1(keep); r = rL1(keep); label = "Max - L1";
    end
    [mu, sd, lo, hi, loOneSided, pValue] = provisionalT(r);
    group(end+1,1) = groups(g);
    comparison(end+1,1) = label;
    nSessions(end+1,1) = sum(keep);
    nDividedTrials(end+1,1) = sum(n(keep));
    nPositive(end+1,1) = sum(d > 0);
    meanPer1000(end+1,1) = mu;
    sdPer1000(end+1,1) = sd;
    lower95(end+1,1) = lo;
    upper95(end+1,1) = hi;
    lowerOneSided95(end+1,1) = loOneSided;
    oneSidedP(end+1,1) = pValue;
    totalDeltaLL(end+1,1) = sum(d);
    pooledPer1000(end+1,1) = 1000 * sum(d) / sum(n(keep));
  end
end
groupTable = table(group, comparison, nSessions, nDividedTrials, ...
  nPositive, meanPer1000, sdPer1000, lower95, upper95, ...
  lowerOneSided95, oneSidedP, ...
  totalDeltaLL, pooledPer1000);

summary = struct();
summary.sourceModelResult = result.createdBy;
summary.sessionTable = sessionTable;
summary.groupTable = groupTable;
summary.inferenceNote = ['One-sided session-level t test of mean LL advantage ' ...
  'per 1000 divided trials > 0; provisional, no correction for repeated looks.'];

fig = figure('Color', 'w', 'Visible', char(p.Results.Visible));
layout = tiledlayout(fig, numel(animals), 2, ...
  'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, 'IDL held-out model comparison by session');
for g = 1:numel(animals)
  keep = animal == animals(g);
  x = 1:sum(keep);
  labels = ids(keep);
  for j = 1:numel(labels)
    date = regexp(char(labels(j)), '[0-9]{8}$', 'match', 'once');
    labels(j) = string(date(5:6)) + "/" + string(date(7:8));
  end

  ax = nexttile(layout);
  plot(ax, x, rL2(keep), 'o-', 'LineWidth', 1.5, ...
    'DisplayName', 'Max - L2');
  hold(ax, 'on');
  plot(ax, x, rL1(keep), 's-', 'LineWidth', 1.2, ...
    'DisplayName', 'Max - L1');
  yline(ax, 0, ':', 'HandleVisibility', 'off');
  set(ax, 'XTick', x, 'XTickLabel', labels);
  xlabel(ax, 'Session date');
  ylabel(ax, 'LL advantage per 1,000 divided trials');
  title(ax, sprintf('%s: individual sessions', char(animals(g))));
  grid(ax, 'on'); box(ax, 'off'); legend(ax, 'Location', 'best');

  ax = nexttile(layout);
  plot(ax, x, cumsum(dL2(keep)), 'o-', 'LineWidth', 1.5, ...
    'DisplayName', 'Max - L2');
  hold(ax, 'on');
  plot(ax, x, cumsum(dL1(keep)), 's-', 'LineWidth', 1.2, ...
    'DisplayName', 'Max - L1');
  yline(ax, 0, ':', 'HandleVisibility', 'off');
  set(ax, 'XTick', x, 'XTickLabel', labels);
  xlabel(ax, 'Session date');
  ylabel(ax, 'Cumulative held-out LL advantage');
  title(ax, sprintf('%s: cumulative scores (descriptive)', char(animals(g))));
  grid(ax, 'on'); box(ax, 'off'); legend(ax, 'Location', 'best');
end

if strlength(string(p.Results.PDFPath)) > 0
  exportgraphics(fig, char(p.Results.PDFPath), 'ContentType', 'vector');
end
end

function [mu, sd, lo, hi, loOneSided, pValue] = provisionalT(x)
% Student-t calculations using base MATLAB's betainc and fzero.
x = double(x(:));
n = numel(x);
mu = mean(x);
if n < 2
  sd = NaN; lo = NaN; hi = NaN;
  loOneSided = NaN; pValue = NaN;
  return;
end
sd = std(x, 0);
se = sd / sqrt(n);
if se == 0
  lo = mu; hi = mu; loOneSided = mu;
  pValue = double(mu < 0) + 0.5 * double(mu == 0);
  return;
end
df = n - 1;
t = mu / se;
pValue = tUpperTail(t, df);
tCritical = fzero(@(z) tUpperTail(z, df) - 0.025, [0 1000]);
lo = mu - tCritical * se;
hi = mu + tCritical * se;
tOneSided = fzero(@(z) tUpperTail(z, df) - 0.05, [0 1000]);
loOneSided = mu - tOneSided * se;
end

function p = tUpperTail(t, df)
halfTail = 0.5 * betainc(df / (df + t.^2), df/2, 0.5);
if t >= 0
  p = halfTail;
else
  p = 1 - halfTail;
end
end
