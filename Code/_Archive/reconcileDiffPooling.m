function out = reconcileDiffPooling(summary, stepType)
% reconcileDiffPooling
%
% Reconcile pooled diff kernels under:
%   1) common session weights
%   2) kernelAverage-style separate weights by sideType x streamType
%
% Inputs:
%   summary   : table from collectKernelSessionSummary()
%   stepType  : 1=inc, 2=dec   (default = 1)
%
% Output:
%   out : struct with pooled vectors, scales, and reconstruction errors
%
% Purpose:
%   Test whether:
%       pooled(diff) == pooled(change) - pooled(noChange)
%   under different pooling rules.
%
% Key point:
%   With a single common session weight, the identity should hold.
%   With kernelAverage-style separate weights, it generally need not hold.

if nargin < 2 || isempty(stepType)
  stepType = 1;
end

% --- pull rows ---
rowsDiff   = summary(summary.sideType == 1 & summary.stepType == stepType, :);
rowsChange = summary(summary.sideType == 2 & summary.stepType == stepType, :);
rowsNoCh   = summary(summary.sideType == 3 & summary.stepType == stepType, :);

if isempty(rowsDiff) || isempty(rowsChange) || isempty(rowsNoCh)
  error('reconcileDiffPooling:MissingRows', ...
    'Missing diff/change/noChange rows for stepType=%d.', stepType);
end

% Sort to align sessions
rowsDiff   = sortrows(rowsDiff,   'sessionIndex');
rowsChange = sortrows(rowsChange, 'sessionIndex');
rowsNoCh   = sortrows(rowsNoCh,   'sessionIndex');

if any(rowsDiff.sessionIndex ~= rowsChange.sessionIndex) || ...
   any(rowsDiff.sessionIndex ~= rowsNoCh.sessionIndex)
  error('reconcileDiffPooling:SessionMismatch', ...
    'Session indices do not align across diff/change/noChange.');
end

nSessions = height(rowsDiff);

out = struct;
out.stepType = stepType;
out.nSessions = nSessions;

% ============================================================
% A. COMMON-WEIGHT pooling
% ============================================================
commonModes = {'equal', 'prefenergy', 'ntrials'};

for iMode = 1:numel(commonModes)
  mode = commonModes{iMode};

  w = getCommonWeights(rowsChange, mode);
  valid = isfinite(w) & w > 0;

  rd = rowsDiff(valid,:);
  rc = rowsChange(valid,:);
  rn = rowsNoCh(valid,:);
  w  = w(valid);

  w = w / sum(w);

  pD_direct = weightedAverageVecCell(rd.prefVec,   w);
  qD_direct = weightedAverageVecCell(rd.probeVec,  w);

  pC = weightedAverageVecCell(rc.prefVec, w);
  qC = weightedAverageVecCell(rc.probeVec, w);

  pN = weightedAverageVecCell(rn.prefVec, w);
  qN = weightedAverageVecCell(rn.probeVec, w);

  pD_fromParts = pC - pN;
  qD_fromParts = qC - qN;

  S = struct;
  S.weights = w;
  S.nValid = numel(w);

  S.vectors.diff_direct.pref   = pD_direct;
  S.vectors.diff_direct.probe  = qD_direct;
  S.vectors.change.pref        = pC;
  S.vectors.change.probe       = qC;
  S.vectors.noChange.pref      = pN;
  S.vectors.noChange.probe     = qN;
  S.vectors.diff_fromParts.pref  = pD_fromParts;
  S.vectors.diff_fromParts.probe = qD_fromParts;

  S.errors.prefMaxAbs  = max(abs(pD_direct - pD_fromParts));
  S.errors.probeMaxAbs = max(abs(qD_direct - qD_fromParts));

  S.scales.diff_direct    = safeScale(pD_direct, qD_direct);
  S.scales.change         = safeScale(pC, qC);
  S.scales.noChange       = safeScale(pN, qN);
  S.scales.diff_fromParts = safeScale(pD_fromParts, qD_fromParts);

  S.cosines.diff_direct      = safeCosine(pD_direct, qD_direct);
  S.cosines.diff_fromParts   = safeCosine(pD_fromParts, qD_fromParts);
  S.cosines.change           = safeCosine(pC, qC);
  S.cosines.noChange         = safeCosine(pN, qN);

  out.common.(mode) = S;
end

% ============================================================
% B. kernelAverage-style separate weights
%    weights depend on sideType and streamType
% ============================================================
% This mirrors the fact that poolSessionKernels() weights each
% sideType x stepType x streamType separately. :contentReference[oaicite:1]{index=1}

% Direct pooled diff
wd_pref = safeNormalize(1 ./ rowsDiff.prefVar);
wd_probe = safeNormalize(1 ./ rowsDiff.probeVar);

pD_direct = weightedAverageVecCell(rowsDiff.prefVec,  wd_pref);
qD_direct = weightedAverageVecCell(rowsDiff.probeVec, wd_probe);

% Direct pooled change
wc_pref = safeNormalize(1 ./ rowsChange.prefVar);
wc_probe = safeNormalize(1 ./ rowsChange.probeVar);

pC = weightedAverageVecCell(rowsChange.prefVec, wc_pref);
qC = weightedAverageVecCell(rowsChange.probeVec, wc_probe);

% Direct pooled noChange
wn_pref = safeNormalize(1 ./ rowsNoCh.prefVar);
wn_probe = safeNormalize(1 ./ rowsNoCh.probeVar);

pN = weightedAverageVecCell(rowsNoCh.prefVec, wn_pref);
qN = weightedAverageVecCell(rowsNoCh.probeVec, wn_probe);

pD_fromParts = pC - pN;
qD_fromParts = qC - qN;

K = struct;
K.weights.diff.pref     = wd_pref;
K.weights.diff.probe    = wd_probe;
K.weights.change.pref   = wc_pref;
K.weights.change.probe  = wc_probe;
K.weights.noChange.pref = wn_pref;
K.weights.noChange.probe= wn_probe;

K.vectors.diff_direct.pref   = pD_direct;
K.vectors.diff_direct.probe  = qD_direct;
K.vectors.change.pref        = pC;
K.vectors.change.probe       = qC;
K.vectors.noChange.pref      = pN;
K.vectors.noChange.probe     = qN;
K.vectors.diff_fromParts.pref  = pD_fromParts;
K.vectors.diff_fromParts.probe = qD_fromParts;

K.errors.prefMaxAbs  = max(abs(pD_direct - pD_fromParts));
K.errors.probeMaxAbs = max(abs(qD_direct - qD_fromParts));

K.scales.diff_direct    = safeScale(pD_direct, qD_direct);
K.scales.change         = safeScale(pC, qC);
K.scales.noChange       = safeScale(pN, qN);
K.scales.diff_fromParts = safeScale(pD_fromParts, qD_fromParts);

K.cosines.diff_direct      = safeCosine(pD_direct, qD_direct);
K.cosines.diff_fromParts   = safeCosine(pD_fromParts, qD_fromParts);
K.cosines.change           = safeCosine(pC, qC);
K.cosines.noChange         = safeCosine(pN, qN);

out.kernelAverageStyle = K;

% ============================================================
% Print summary
% ============================================================
stepLabel = ternary(stepType == 1, 'inc', 'dec');

fprintf('\n=== reconcileDiffPooling: %s ===\n', stepLabel);
fprintf('nSessions = %d\n', nSessions);

fprintf('\n-- Common weights --\n');
for iMode = 1:numel(commonModes)
  mode = commonModes{iMode};
  S = out.common.(mode);

  fprintf('%-10s : diff_direct=% .4f  diff_fromParts=% .4f  prefErr=%.3g  probeErr=%.3g\n', ...
    mode, S.scales.diff_direct, S.scales.diff_fromParts, ...
    S.errors.prefMaxAbs, S.errors.probeMaxAbs);
end

fprintf('\n-- kernelAverage-style separate weights --\n');
fprintf('diff_direct    = % .4f\n', K.scales.diff_direct);
fprintf('change         = % .4f\n', K.scales.change);
fprintf('noChange       = % .4f\n', K.scales.noChange);
fprintf('diff_fromParts = % .4f\n', K.scales.diff_fromParts);
fprintf('pref max err   = %.3g\n', K.errors.prefMaxAbs);
fprintf('probe max err  = %.3g\n', K.errors.probeMaxAbs);

fprintf('\nCosines (kernelAverage-style):\n');
fprintf('  diff_direct    = % .4f\n', K.cosines.diff_direct);
fprintf('  diff_fromParts = % .4f\n', K.cosines.diff_fromParts);
fprintf('  change         = % .4f\n', K.cosines.change);
fprintf('  noChange       = % .4f\n', K.cosines.noChange);
fprintf('=========================================\n\n');

end


function w = getCommonWeights(rows, mode)
switch lower(mode)
  case 'equal'
    w = ones(height(rows),1);

  case 'prefenergy'
    w = rows.prefEnergy;

  case 'ntrials'
    w = rows.nTrialsStep;

  otherwise
    error('getCommonWeights:BadMode', 'Unknown mode: %s', mode);
end
end


function vAvg = weightedAverageVecCell(vecCell, w)
n = numel(vecCell);
m = numel(vecCell{1});
mat = nan(m, n);

for i = 1:n
  v = vecCell{i};
  v = v(:);
  if numel(v) ~= m
    error('weightedAverageVecCell:LengthMismatch', ...
      'Vector lengths do not match.');
  end
  mat(:,i) = v;
end

vAvg = mat * w(:);
end


function w = safeNormalize(x)
x = x(:);
valid = isfinite(x) & x > 0;
w = zeros(size(x));
if any(valid)
  w(valid) = x(valid) / sum(x(valid));
else
  w(:) = nan;
end
end


function s = safeScale(p, q)
p = p(:);
q = q(:);
den = sum(p.^2);
if den > 0
  s = sum(p .* q) / den;
else
  s = nan;
end
end


function c = safeCosine(x, y)
x = x(:);
y = y(:);
nx = sqrt(sum(x.^2));
ny = sqrt(sum(y.^2));
if nx > 0 && ny > 0
  c = sum(x .* y) / (nx * ny);
else
  c = nan;
end
end


function out = ternary(cond, a, b)
if cond
  out = a;
else
  out = b;
end
end