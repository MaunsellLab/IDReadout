function P = poolStrategyForensics(tblAll, varargin)
% poolStrategyForensics  Fit pooled behavioral strategy models across sessions.
%
% Usage:
%   P = poolStrategyForensics(tblAll)
%
% Input:
%   tblAll   concatenated per-trial table across sessions
%
% Required variables in tblAll for current implementation:
%   isCorrect
%   sumDM
%   longestPosRun
%   longestNegRun
%   nSwitch
%   sessionID     (recommended; created by collectStrategyTables)
%
% Optional name/value:
%   'groupVar'          session grouping variable name (default 'sessionID')
%   'fitSessionFixed'   include session fixed effects (default true)
%   'verbose'           print summary to command window (default true)
%
% Output struct P:
%   P.tblAll
%   P.models
%   P.meta
%   P.summary
%   P.compare
%
% Notes:
%   - Uses binomial GLMs via fitglm.
%   - This first implementation fits:
%         p0: isCorrect ~ 1 [+ sessionID]
%         p1: isCorrect ~ sumDM [+ sessionID]
%         p2: isCorrect ~ sumDM + longestPosRun + longestNegRun + nSwitch [+ sessionID]
%   - Deviance explained is computed relative to p0.
%
% Future:
%   - add p3 and p4
%   - optionally add Jeffreys/Firth-style penalization if you want pooled code
%     to match session code exactly

% ----------------------------
% parse inputs
% ----------------------------
p = inputParser;
p.addRequired('tblAll', @istable);
p.addParameter('groupVar', 'sessionID', @(x) ischar(x) || isstring(x));
p.addParameter('fitSessionFixed', true, @(x) islogical(x) && isscalar(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(tblAll, varargin{:});

groupVar = char(p.Results.groupVar);
fitSessionFixed = p.Results.fitSessionFixed;
verbose = p.Results.verbose;

% ----------------------------
% validate variables
% ----------------------------
requiredVars = {'isCorrect', 'sumDM', 'longestPosRun', 'longestNegRun', 'nSwitch'};
if fitSessionFixed
    requiredVars{end+1} = groupVar;
end

missingVars = requiredVars(~ismember(requiredVars, tblAll.Properties.VariableNames));
if ~isempty(missingVars)
    error('poolStrategyForensics:MissingVars', ...
        'tblAll is missing required variables: %s', strjoin(missingVars, ', '));
end

T = tblAll;

% normalize types
T.isCorrect = localBinary01(T.isCorrect);

if fitSessionFixed
    T.(groupVar) = categorical(T.(groupVar));
end

% remove rows with missing values in variables used here
varsUsed = {'isCorrect', 'sumDM', 'longestPosRun', 'longestNegRun', 'nSwitch'};
if fitSessionFixed
    varsUsed{end+1} = groupVar;
end

keep = true(height(T),1);
for i = 1:numel(varsUsed)
    v = T.(varsUsed{i});
    if iscategorical(v)
        keep = keep & ~ismissing(v);
    else
        keep = keep & ~isnan(v);
    end
end

nDroppedMissing = sum(~keep);
T = T(keep,:);

if isempty(T)
    error('poolStrategyForensics:NoRowsLeft', ...
        'No valid rows remain after excluding missing values.');
end

% ----------------------------
% formulas
% ----------------------------
if fitSessionFixed
    f0 = sprintf('isCorrect ~ 1 + %s', groupVar);
    f1 = sprintf('isCorrect ~ sumDM + %s', groupVar);
    f2 = sprintf('isCorrect ~ sumDM + longestPosRun + longestNegRun + nSwitch + %s', groupVar);
else
    f0 = 'isCorrect ~ 1';
    f1 = 'isCorrect ~ sumDM';
    f2 = 'isCorrect ~ sumDM + longestPosRun + longestNegRun + nSwitch';
end

% ----------------------------
% fit models
% ----------------------------
models = struct();
models.p0 = fitglm(T, f0, 'Distribution', 'binomial');
models.p1 = fitglm(T, f1, 'Distribution', 'binomial');
models.p2 = fitglm(T, f2, 'Distribution', 'binomial');

% ----------------------------
% descriptive summaries
% ----------------------------
summary = struct();

summary.nTrials = height(T);

if fitSessionFixed
    summary.nSessions = numel(categories(T.(groupVar)));
else
    summary.nSessions = NaN;
end

summary.pCorrect = mean(T.isCorrect);
summary.nCorrect = sum(T.isCorrect == 1);
summary.nError = sum(T.isCorrect == 0);
summary.nDroppedMissing = nDroppedMissing;

% correct-error differences for included predictors
predNames = {'sumDM', 'longestPosRun', 'longestNegRun', 'nSwitch'};
for i = 1:numel(predNames)
    vn = predNames{i};
    x = T.(vn);
    summary.([vn '_meanCorrect']) = mean(x(T.isCorrect == 1), 'omitnan');
    summary.([vn '_meanError'])   = mean(x(T.isCorrect == 0), 'omitnan');
    summary.([vn '_diffCorrectError']) = ...
        summary.([vn '_meanCorrect']) - summary.([vn '_meanError']);
end

% deviance and deviance explained
summary.deviance_p0 = models.p0.Deviance;
summary.deviance_p1 = models.p1.Deviance;
summary.deviance_p2 = models.p2.Deviance;

summary.devExplained_p1 = localDevianceExplained(models.p0, models.p1);
summary.devExplained_p2 = localDevianceExplained(models.p0, models.p2);

% coefficients of interest
summary.beta_sumDM_p1 = localCoeff(models.p1, 'sumDM');
summary.beta_sumDM_p2 = localCoeff(models.p2, 'sumDM');
summary.beta_longestPosRun_p2 = localCoeff(models.p2, 'longestPosRun');
summary.beta_longestNegRun_p2 = localCoeff(models.p2, 'longestNegRun');
summary.beta_nSwitch_p2 = localCoeff(models.p2, 'nSwitch');

% ----------------------------
% model comparisons
% ----------------------------
compare = struct();
compare.p1_vs_p0 = localCompareNested(models.p0, models.p1);
compare.p2_vs_p1 = localCompareNested(models.p1, models.p2);
compare.p2_vs_p0 = localCompareNested(models.p0, models.p2);

% ----------------------------
% meta
% ----------------------------
meta = struct();
meta.groupVar = groupVar;
meta.fitSessionFixed = fitSessionFixed;
meta.formula_p0 = f0;
meta.formula_p1 = f1;
meta.formula_p2 = f2;
meta.predictors = predNames;

% carry through a few convenient pieces if present
carryFields = {'postFrames', 'msPerVFrame', 'proj'};
for i = 1:numel(carryFields)
    fn = carryFields{i};
    if ismember(fn, T.Properties.VariableNames)
        meta.(fn) = T.(fn);
    elseif ismember(fn, tblAll.Properties.VariableNames)
        meta.(fn) = tblAll.(fn);
    end
end

% ----------------------------
% output struct
% ----------------------------
P = struct();
P.tblAll = T;
P.models = models;
P.meta = meta;
P.summary = summary;
P.compare = compare;

% ----------------------------
% print
% ----------------------------
if verbose
    fprintf('\nPooled strategy analysis\n');
    fprintf('  Trials:    %d\n', summary.nTrials);
    if fitSessionFixed
        fprintf('  Sessions:  %d\n', summary.nSessions);
    end
    fprintf('  Correct:   %.2f%%\n', 100 * summary.pCorrect);
    if summary.nDroppedMissing > 0
        fprintf('  Dropped missing rows: %d\n', summary.nDroppedMissing);
    end

    fprintf('\nDescriptive correct-error differences:\n');
    fprintf('  sumDM           %+0.6g\n', summary.sumDM_diffCorrectError);
    fprintf('  longestPosRun   %+0.6g\n', summary.longestPosRun_diffCorrectError);
    fprintf('  longestNegRun   %+0.6g\n', summary.longestNegRun_diffCorrectError);
    fprintf('  nSwitch         %+0.6g\n', summary.nSwitch_diffCorrectError);

    fprintf('\nModel deviances:\n');
    fprintf('  p0   %0.6f\n', summary.deviance_p0);
    fprintf('  p1   %0.6f   (dev explained = %.4f)\n', ...
        summary.deviance_p1, summary.devExplained_p1);
    fprintf('  p2   %0.6f   (dev explained = %.4f)\n', ...
        summary.deviance_p2, summary.devExplained_p2);

    fprintf('\nKey coefficients:\n');
    fprintf('  p1: sumDM            %+0.6g\n', summary.beta_sumDM_p1);
    fprintf('  p2: sumDM            %+0.6g\n', summary.beta_sumDM_p2);
    fprintf('  p2: longestPosRun    %+0.6g\n', summary.beta_longestPosRun_p2);
    fprintf('  p2: longestNegRun    %+0.6g\n', summary.beta_longestNegRun_p2);
    fprintf('  p2: nSwitch          %+0.6g\n', summary.beta_nSwitch_p2);

    fprintf('\nNested comparisons:\n');
    fprintf('  p1 vs p0:  dDev = %0.6f, dDF = %d, p = %.6g\n', ...
        compare.p1_vs_p0.deltaDeviance, compare.p1_vs_p0.deltaDF, compare.p1_vs_p0.pValue);
    fprintf('  p2 vs p1:  dDev = %0.6f, dDF = %d, p = %.6g\n', ...
        compare.p2_vs_p1.deltaDeviance, compare.p2_vs_p1.deltaDF, compare.p2_vs_p1.pValue);
    fprintf('  p2 vs p0:  dDev = %0.6f, dDF = %d, p = %.6g\n', ...
        compare.p2_vs_p0.deltaDeviance, compare.p2_vs_p0.deltaDF, compare.p2_vs_p0.pValue);
end

end


% ========================================================================
function y = localBinary01(y)
% convert logical/numeric/categorical binary response to 0/1 double column

if islogical(y)
    y = double(y);
elseif isnumeric(y)
    y = double(y);
elseif iscategorical(y) || isstring(y) || ischar(y)
    y = double(categorical(y));
    y = y - min(y);
else
    error('poolStrategyForensics:BadResponseType', ...
        'Unsupported type for isCorrect.');
end

y = y(:);

u = unique(y(~isnan(y)));
if ~all(ismember(u, [0 1]))
    error('poolStrategyForensics:BadResponseValues', ...
        'isCorrect must contain only binary values 0/1 (or logical equivalent).');
end
end


% ========================================================================
function devExp = localDevianceExplained(m0, m1)
% fraction of null deviance reduced by more complex model
dev0 = m0.Deviance;
dev1 = m1.Deviance;

if dev0 <= 0
    devExp = NaN;
else
    devExp = (dev0 - dev1) / dev0;
end
end


% ========================================================================
function c = localCoeff(mdl, termName)
% return coefficient estimate for termName, NaN if absent
c = NaN;
coefNames = mdl.CoefficientNames;
idx = strcmp(coefNames, termName);
if any(idx)
    c = mdl.Coefficients.Estimate(find(idx,1,'first'));
end
end


% ========================================================================
function C = localCompareNested(mSmall, mBig)
% Likelihood-ratio / deviance comparison for nested GLMs

C = struct();

C.devianceSmall = mSmall.Deviance;
C.devianceBig = mBig.Deviance;
C.dfSmall = mSmall.DFE;
C.dfBig = mBig.DFE;

C.deltaDeviance = C.devianceSmall - C.devianceBig;
C.deltaDF = C.dfSmall - C.dfBig;

if C.deltaDF > 0 && C.deltaDeviance >= 0
    C.pValue = 1 - chi2cdf(C.deltaDeviance, C.deltaDF);
else
    C.pValue = NaN;
end
end