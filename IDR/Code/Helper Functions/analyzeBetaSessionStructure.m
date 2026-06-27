function structureAnalysis = analyzeBetaSessionStructure()
% analyzeBetaSessionStructure
% Test whether session-to-session beta heterogeneity is associated with
% identifiable session characteristics.
%
% Predictors examined:
%   - recording date
%   - signal coherence used on preferred-noise trials
%   - preferred-noise amplitude
%   - fraction correct
%   - trial count
%   - error count
%   - number of probe directions in the parent session
%
% Uses both raw session beta and random-effects shrunken beta.
%
% Reads:
%   Data/FullSessions/BetaAnalysis/AcrossSessions/BetaSessionFitAnalysis.mat
%   Data/FullSessions/BetaAnalysis/SessionData/*.mat
%
% Saves:
%   Data/FullSessions/BetaAnalysis/AcrossSessions/BetaSessionStructureAnalysis.mat

% cleanupObj = initProjectPath(); %#ok<NASGU>

baseFolder = domainFolder(mfilename('fullpath'));
sessionFolder = fullfile(baseFolder, 'Data', 'FullSessions', 'BetaAnalysis');
acrossFolder = validFolder(fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries'));

inputPath = fullfile(acrossFolder, 'BetaSessionFitAnalysis.mat');
S = load(inputPath, 'sessionBetaAnalysis');

A = S.sessionBetaAnalysis;
T = A.sessionTable;

n = height(T);

dateNum = nan(n,1);
dateText = strings(n,1);
signalCohPC = nan(n,1);
prefCohNoisePC = nan(n,1);
nProbeDirections = nan(n,1);

for i = 1:n
  filePath = fullfile(sessionFolder, T.fileName{i});
  D = load(filePath, 'sessionNoise');
  N = D.sessionNoise;
  H = N.sessionHeader;

  % Date
  if isfield(H, 'date')
    d = localValue(H.date);
    try
      if isdatetime(d)
        dt = d;
      elseif isnumeric(d)
        dt = datetime(d, 'ConvertFrom', 'datenum');
      else
        dt = datetime(string(d));
      end
      dateNum(i) = datenum(dt);
      dateText(i) = string(dt, 'yyyy-MM-dd');
    catch
      dateNum(i) = NaN;
      dateText(i) = "";
    end
  end

  % Preferred-noise trials should use one signal coherence within session.
  use = logical(N.hasPreferredNoise(:));
  sig = unique(double(N.signalCohPC(use)));
  if numel(sig) == 1
    signalCohPC(i) = sig;
  end

  if isfield(H, 'prefCohNoisePC')
    prefCohNoisePC(i) = double(localValue(H.prefCohNoisePC));
  end

  if isfield(H, 'nProbeDirections')
    nProbeDirections(i) = double(localValue(H.nProbeDirections));
  end
end

T.dateNum = dateNum;
T.dateText = dateText;
T.signalCohPC = signalCohPC;
T.prefCohNoisePC = prefCohNoisePC;
T.nProbeDirections = nProbeDirections;

predictorNames = { ...
  'dateNum', ...
  'signalCohPC', ...
  'prefCohNoisePC', ...
  'fractionCorrect', ...
  'nTrials', ...
  'nError', ...
  'nProbeDirections'};

results = table('Size',[numel(predictorNames) 9], ...
  'VariableTypes', {'string','double','double','double','double', ...
                    'double','double','double','double'}, ...
  'VariableNames', {'predictor','rhoRaw','pRaw','rhoShrunken', ...
                    'pShrunken','slopeRaw','slopeRawSE', ...
                    'slopeShrunken','slopeShrunkenSE'});

for i = 1:numel(predictorNames)
  name = predictorNames{i};
  x = T.(name);

  [rhoRaw,pRaw] = corr(x,T.beta,'Type','Spearman','Rows','complete');
  [rhoShr,pShr] = corr(x,T.shrunkenBeta,'Type','Spearman','Rows','complete');

  [bRaw,seRaw] = simpleSlope(x,T.beta);
  [bShr,seShr] = simpleSlope(x,T.shrunkenBeta);

  results.predictor(i) = string(name);
  results.rhoRaw(i) = rhoRaw;
  results.pRaw(i) = pRaw;
  results.rhoShrunken(i) = rhoShr;
  results.pShrunken(i) = pShr;
  results.slopeRaw(i) = bRaw;
  results.slopeRawSE(i) = seRaw;
  results.slopeShrunken(i) = bShr;
  results.slopeShrunkenSE(i) = seShr;
end

% Multiple regression on standardized continuous predictors.
Xvars = {'dateNum','signalCohPC','prefCohNoisePC', ...
         'fractionCorrect','nTrials','nProbeDirections'};

X = [];
validNames = strings(0,1);

for i = 1:numel(Xvars)
  x = T.(Xvars{i});
  if numel(unique(x(isfinite(x)))) > 1
    X = [X, zscoreMissing(x)]; %#ok<AGROW>
    validNames(end+1,1) = string(Xvars{i}); %#ok<AGROW>
  end
end

use = all(isfinite(X),2) & isfinite(T.beta);
mdlRaw = fitlm(X(use,:), T.beta(use), ...
  'VarNames', [cellstr(validNames); {'beta'}]);

useShr = all(isfinite(X),2) & isfinite(T.shrunkenBeta);
mdlShrunken = fitlm(X(useShr,:), T.shrunkenBeta(useShr), ...
  'VarNames', [cellstr(validNames); {'shrunkenBeta'}]);

structureAnalysis = struct();
structureAnalysis.version = 1;
structureAnalysis.sessionTable = T;
structureAnalysis.univariateResults = results;
structureAnalysis.multipleRegressionRaw = mdlRaw;
structureAnalysis.multipleRegressionShrunken = mdlShrunken;
structureAnalysis.predictorNames = validNames;
structureAnalysis.createdBy = mfilename;
structureAnalysis.createdDate = datetime('now');

outputPath = fullfile(acrossFolder, 'BetaSessionStructureAnalysis.mat');
save(outputPath, 'structureAnalysis', '-v7.3');

disp(results);
fprintf('\nMultiple regression on raw beta:\n');
disp(mdlRaw);
fprintf('\nMultiple regression on shrunken beta:\n');
disp(mdlShrunken);

plotStructure(T);
end

% -------------------------------------------------------------------------
function value = localValue(x)

value = x;
while isstruct(value) && isfield(value,'data')
  value = value.data;
end

if isnumeric(value) && ~isscalar(value)
  value = value(1);
end
end

% -------------------------------------------------------------------------
function z = zscoreMissing(x)

x = double(x(:));
z = nan(size(x));
use = isfinite(x);

mu = mean(x(use));
sd = std(x(use));

if sd > 0
  z(use) = (x(use)-mu)/sd;
end
end

% -------------------------------------------------------------------------
function [slope,se] = simpleSlope(x,y)

use = isfinite(x) & isfinite(y);
x = x(use);
y = y(use);

if numel(x) < 3 || numel(unique(x)) < 2
  slope = NaN;
  se = NaN;
  return;
end

X = [ones(numel(x),1), x];
b = X\y;
resid = y-X*b;
sigma2 = sum(resid.^2)/(numel(y)-2);
covB = sigma2*inv(X'*X);

slope = b(2);
se = sqrt(covB(2,2));
end

% -------------------------------------------------------------------------
function plotStructure(T)

figure;
scatter(T.dateNum,T.beta,'filled');
hold on
yline(T.referenceBeta(1),'r--');
datetick('x','yyyy-mm');
xlabel('Recording date');
ylabel('Session beta');
title('Session beta over time');
box off

figure;
scatter(T.signalCohPC,T.beta,'filled');
hold on
yline(T.referenceBeta(1),'r--');
xlabel('Signal coherence on noise trials (%)');
ylabel('Session beta');
title('Session beta versus signal coherence');
box off

figure;
scatter(T.prefCohNoisePC,T.beta,'filled');
hold on
yline(T.referenceBeta(1),'r--');
xlabel('Preferred-noise amplitude (%)');
ylabel('Session beta');
title('Session beta versus noise amplitude');
box off

figure;
scatter(T.nProbeDirections,T.beta,'filled');
hold on
yline(T.referenceBeta(1),'r--');
xlabel('Number of probe directions');
ylabel('Session beta');
title('Session beta versus probe context');
box off
end
