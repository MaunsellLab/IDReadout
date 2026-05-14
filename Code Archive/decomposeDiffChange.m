function out = decomposeDiffChange(summary, stepType, poolMode)
% decomposeDiffChange
%
% Decompose pooled diff kernel into change/noChange cross-terms.
%
% Inputs:
%   summary   : table from collectKernelSessionSummary()
%   stepType  : 1=inc, 2=dec
%   poolMode  : 'equal', 'prefenergy', 'invvar', or 'ntrials'
%
% Output:
%   out       : struct containing pooled vectors, decomposition terms,
%               and resulting scales
%
% Notes:
%   sideType: 1=diff, 2=change, 3=noChange, 4=RF, 5=Opp
%   streamType in summary vectors:
%       prefVec  = stream 1
%       probeVec = stream 2

if nargin < 2 || isempty(stepType)
  stepType = 1;
end
if nargin < 3 || isempty(poolMode)
  poolMode = 'equal';
end

% --- pull rows for this stepType ---
idxChange = summary.sideType == 2 & summary.stepType == stepType;
idxNoCh   = summary.sideType == 3 & summary.stepType == stepType;
idxDiff   = summary.sideType == 1 & summary.stepType == stepType;

changeRows = summary(idxChange, :);
noChRows   = summary(idxNoCh, :);
diffRows   = summary(idxDiff, :); %#ok<NASGU>

if height(changeRows) ~= height(noChRows)
  error('decomposeDiffChange:Mismatch', ...
    'change and noChange row counts do not match.');
end

nSessions = height(changeRows);
if nSessions == 0
  error('decomposeDiffChange:NoRows', ...
    'No rows found for stepType=%d.', stepType);
end

% --- check session ordering matches ---
if any(changeRows.sessionIndex ~= noChRows.sessionIndex)
  error('decomposeDiffChange:SessionOrder', ...
    'change and noChange rows are not aligned by sessionIndex.');
end

% --- get weights ---
w = nan(nSessions,1);
for i = 1:nSessions
  switch lower(poolMode)
    case 'equal'
      w(i) = 1;

    case 'prefenergy'
      w(i) = changeRows.prefEnergy(i);

    case 'invvar'
      if isfinite(changeRows.poolWeightPref(i)) && changeRows.poolWeightPref(i) > 0
        w(i) = changeRows.poolWeightPref(i);
      else
        w(i) = nan;
      end

    case 'ntrials'
      w(i) = changeRows.nTrialsStep(i);

    otherwise
      error('decomposeDiffChange:BadPoolMode', ...
        'Unknown poolMode: %s', poolMode);
  end
end

valid = isfinite(w) & w > 0;
changeRows = changeRows(valid,:);
noChRows   = noChRows(valid,:);
w          = w(valid);

if isempty(w)
  error('decomposeDiffChange:NoValidWeights', ...
    'No valid sessions after applying poolMode=%s.', poolMode);
end

w = w / sum(w);

% --- pool kernels ---
pC = weightedAverageVecCell(changeRows.prefVec,  w);
qC = weightedAverageVecCell(changeRows.probeVec, w);

pN = weightedAverageVecCell(noChRows.prefVec,    w);
qN = weightedAverageVecCell(noChRows.probeVec,   w);

pD_fromParts = pC - pN;
qD_fromParts = qC - qN;

% Initialize direct check (always exists)
out.directCheck.prefMaxAbsErr  = nan;
out.directCheck.probeMaxAbsErr = nan;

% Optional consistency check against directly pooled diff rows
if any(summary.sideType == 1 & summary.stepType == stepType)
  diffRows = summary(summary.sideType == 1 & summary.stepType == stepType, :);
  diffRows = diffRows(valid,:);
  pD_direct = weightedAverageVecCell(diffRows.prefVec,  w);
  qD_direct = weightedAverageVecCell(diffRows.probeVec, w);

  out.directCheck.prefMaxAbsErr  = max(abs(pD_direct - pD_fromParts));
  out.directCheck.probeMaxAbsErr = max(abs(qD_direct - qD_fromParts));
else
  out.directCheck.prefMaxAbsErr  = nan;
  out.directCheck.probeMaxAbsErr = nan;
end

% --- numerator decomposition ---
term_cc = sum(pC .* qC);
term_cn = sum(pC .* qN);
term_nc = sum(pN .* qC);
term_nn = sum(pN .* qN);

num_change = term_cc;
num_diff   = term_cc - term_cn - term_nc + term_nn;

% --- denominator decomposition ---
den_cc = sum(pC .* pC);
den_cn = sum(pC .* pN);
den_nn = sum(pN .* pN);

den_change = den_cc;
den_diff   = den_cc - 2*den_cn + den_nn;

% --- scales ---
scale_change = num_change / den_change;
scale_noCh   = sum(pN .* qN) / sum(pN .* pN);
scale_diff   = num_diff / den_diff;

% --- useful normalized summaries ---
% --- useful normalized summaries ---
out = struct;
out.stepType = stepType;
out.poolMode = poolMode;
out.nSessions = numel(w);

% initialize direct-check fields so printing is always safe
out.directCheck.prefMaxAbsErr  = nan;
out.directCheck.probeMaxAbsErr = nan;
out.pooled.change.pref = pC;
out.pooled.change.probe = qC;
out.pooled.noChange.pref = pN;
out.pooled.noChange.probe = qN;
out.pooled.diff.pref = pD_fromParts;
out.pooled.diff.probe = qD_fromParts;

out.numTerms.cc = term_cc;
out.numTerms.cn = term_cn;
out.numTerms.nc = term_nc;
out.numTerms.nn = term_nn;
out.numTerms.change = num_change;
out.numTerms.diff   = num_diff;

out.denTerms.cc = den_cc;
out.denTerms.cn = den_cn;
out.denTerms.nn = den_nn;
out.denTerms.change = den_change;
out.denTerms.diff   = den_diff;

out.scales.change   = scale_change;
out.scales.noChange = scale_noCh;
out.scales.diff     = scale_diff;

out.relativeNumerator.cn_over_cc = term_cn / term_cc;
out.relativeNumerator.nc_over_cc = term_nc / term_cc;
out.relativeNumerator.nn_over_cc = term_nn / term_cc;

out.relativeDenominator.cn_over_cc = den_cn / den_cc;
out.relativeDenominator.nn_over_cc = den_nn / den_cc;

out.cosines.change_pref_vs_noCh_pref = ...
  safeCosine(pC, pN);
out.cosines.change_probe_vs_noCh_probe = ...
  safeCosine(qC, qN);
out.cosines.change_pref_vs_change_probe = ...
  safeCosine(pC, qC);
out.cosines.noCh_pref_vs_noCh_probe = ...
  safeCosine(pN, qN);
out.cosines.diff_pref_vs_diff_probe = ...
  safeCosine(pD_fromParts, qD_fromParts);

% --- compact printout ---
stepLabel = ternary(stepType == 1, 'inc', 'dec');

fprintf('\n=== decomposeDiffChange: %s, poolMode=%s ===\n', stepLabel, poolMode);
fprintf('nSessions = %d\n', out.nSessions);

fprintf('\nScales:\n');
fprintf('  change   = % .4f\n', out.scales.change);
fprintf('  noChange = % .4f\n', out.scales.noChange);
fprintf('  diff     = % .4f\n', out.scales.diff);

fprintf('\nNumerator terms:\n');
fprintf('  cc = <pC,qC>  = % .6f\n', term_cc);
fprintf('  cn = <pC,qN>  = % .6f\n', term_cn);
fprintf('  nc = <pN,qC>  = % .6f\n', term_nc);
fprintf('  nn = <pN,qN>  = % .6f\n', term_nn);
fprintf('  change num    = % .6f\n', num_change);
fprintf('  diff num      = % .6f\n', num_diff);

fprintf('\nDenominator terms:\n');
fprintf('  cc = <pC,pC>  = % .6f\n', den_cc);
fprintf('  cn = <pC,pN>  = % .6f\n', den_cn);
fprintf('  nn = <pN,pN>  = % .6f\n', den_nn);
fprintf('  change den    = % .6f\n', den_change);
fprintf('  diff den      = % .6f\n', den_diff);

fprintf('\nCosines:\n');
fprintf('  cos(pC,pN)    = % .4f\n', out.cosines.change_pref_vs_noCh_pref);
fprintf('  cos(qC,qN)    = % .4f\n', out.cosines.change_probe_vs_noCh_probe);
fprintf('  cos(pC,qC)    = % .4f\n', out.cosines.change_pref_vs_change_probe);
fprintf('  cos(pN,qN)    = % .4f\n', out.cosines.noCh_pref_vs_noCh_probe);
fprintf('  cos(pD,qD)    = % .4f\n', out.cosines.diff_pref_vs_diff_probe);

fprintf('\nDirect diff reconstruction check:\n');
fprintf('  max abs err pref  = %.3g\n', out.directCheck.prefMaxAbsErr);
fprintf('  max abs err probe = %.3g\n', out.directCheck.probeMaxAbsErr);
fprintf('=============================================\n\n');

end


function vAvg = weightedAverageVecCell(vecCell, w)
% vecCell: Nx1 cell array, each cell contains column vector of same length
n = numel(vecCell);
firstVec = vecCell{1};
m = numel(firstVec);

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