function [M,fig] = runIDSGaussianProbeGainModel(varargin)
% runIDSGaussianProbeGainModel  Predict IDS probe gain from Gaussian competition.
%
% This is a restricted model with one drift mechanism and one effective
% probe competitor.  The competitor may be centered on either probe or
% on the opponent (180-degree) axis, whichever has the larger mean
% response to the two probe pedestals.  It does not use a full direction
% bank and does not fit prestep noise.
%
% The 120-degree IDS gain fit supplies three empirical anchors:
%   direct probe term  = gNP
%   positive win term  = gCP-gNP
%   preferred response = gCD
%
% MT direction tuning determines direct preferred/opponent terms.  The two
% probe pedestals sum through the same Gaussian tuning as their directions
% converge.  Equal-variance Gaussian competition converts the resulting
% drift/probe mean separation into a probe-win probability.  The unknown
% Gaussian SD is fixed by the assumed per-pedestal win probability at 120
% degrees, rather than by fitting probe gains at other directions.
%
% Important scope: with PreferredTuningPower=1, the near-drift component is
% the literal MT-width prediction.  Values greater than one provide a
% transparent way to test later whether max selection sharpens that side of
% the curve.  The initial 180-degree prediction does not depend materially
% on this parameter.
%
% Default use:
%   runIDSGaussianProbeGainModel
%   runIDSGaussianProbeGainModel('IncludeOppositeAxisProbeMechanism',false)
% compares the earlier geometry in which the competitor is always centered
% on one of the probes.  Neither version fits gains at the other offsets.
%
% Existing Data/Probe*/AcrossSessionSummaries/IDS_Probe*_... MAT files are
% loaded automatically to display every analyzed IDS probe offset.  The
% specialized IDS-only 120-degree analysis and explicit ObservedProbeGains
% rows take precedence over older across-offset summaries.
%
% Supply updated threshold estimates as [angle thresholdPC] rows:
%   runIDSGaussianProbeGainModel('AdditionalThresholdAnchors', ...
%     [45 16.2; 120 12.3; 170 11.2])

p = inputParser;
p.addParameter('SummaryFile','',@(x) ischar(x) || isstring(x));
p.addParameter('AnglesDeg',(0:180)',@(x) isnumeric(x) && isvector(x));
p.addParameter('MTSigmaDeg',37.5,@(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('PreferredTuningPower',1,@(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('IncludeOppositeAxisProbeMechanism',true, ...
  @(x) islogical(x) && isscalar(x));
p.addParameter('PedestalWinRate120',0.05,@(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1/6);
p.addParameter('PedestalWinRateRange',[0.04 0.06], ...
  @(x) isnumeric(x) && numel(x) == 2 && all(x > 0) && all(x < 1/6));
p.addParameter('AdditionalThresholdAnchors',[170 11.2], ...
  @(x) isnumeric(x) && (isempty(x) || size(x,2) == 2));
p.addParameter('ObservedProbeGains',[170 0.51 0.49 0.60], ...
  @(x) isnumeric(x) && (isempty(x) || size(x,2) == 4));
p.addParameter('ExportPDF',true,@(x) islogical(x) && isscalar(x));
p.addParameter('SaveResults',true,@(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;
codeVersion = '2026-09-18-opposite-axis-v4';

root = domainFolder(mfilename('fullpath'));
if strlength(string(opts.SummaryFile)) == 0
  summaryFile = fullfile(root,'Data','Probe120','AcrossSessionSummaries', ...
    'IDS120CN_IDS_All_AcrossSessionSummary.mat');
else
  summaryFile = char(opts.SummaryFile);
end
if ~isfile(summaryFile)
  error('runIDSGaussianProbeGainModel:MissingSummary', ...
    'IDS 120-degree summary not found: %s',summaryFile);
end
Q = load(summaryFile,'acrossSummary');
if ~isfield(Q,'acrossSummary')
  error('runIDSGaussianProbeGainModel:BadSummary', ...
    '%s does not contain acrossSummary.',summaryFile);
end
S = Q.acrossSummary;
validateSummary(S,summaryFile);

fit = S.noiseGain.sideRole;
iCD = predictorIndex(fit,'C Drift');
iCP = predictorIndex(fit,'C Probe');
iNP = predictorIndex(fit,'N Probe');
gCD = fit.gain(iCD);
gCP = fit.gain(iCP);
gNP = fit.gain(iNP);
direct120 = gNP;
win120 = gCP-gNP;
if win120 <= 0
  error('runIDSGaussianProbeGainModel:NonpositiveWinAnchor', ...
    'The IDS C Probe - N Probe anchor must be positive; found %.4g.',win120);
end

angles = double(opts.AnglesDeg(:));
if any(~isfinite(angles)) || any(angles < 0 | angles > 180)
  error('runIDSGaussianProbeGainModel:BadAngles', ...
    'AnglesDeg must contain finite values between 0 and 180 degrees.');
end
sigmaMT = double(opts.MTSigmaDeg);
pedestalPC = median(S.trialTable.cohNoisePC,'omitnan');
hitRate = mean(S.trialTable.correct,'omitnan');
threshold120PC = physicalThresholdPC(S);
offsetSummaries = loadIDSOffsetSummaries(root);
thresholdAnchors = combineThresholdAnchors(offsetSummaries.thresholdRows, ...
  threshold120PC,opts.AdditionalThresholdAnchors);
thresholdPC = thresholdAtAngles(angles,thresholdAnchors);

modelInputs = struct('gCD',gCD,'gCP',gCP,'gNP',gNP, ...
  'direct120',direct120,'win120',win120,'sigmaMT',sigmaMT, ...
  'preferredTuningPower',opts.PreferredTuningPower, ...
  'includeOppositeAxisProbeMechanism',opts.IncludeOppositeAxisProbeMechanism, ...
  'pedestalPC',pedestalPC,'threshold120PC',threshold120PC, ...
  'thresholdAnchors',thresholdAnchors,'hitRate',hitRate);

central = computeCurve(angles,thresholdPC,modelInputs,opts.PedestalWinRate120);
singleProbeInputs = modelInputs;
singleProbeInputs.includeOppositeAxisProbeMechanism = false;
singleProbeCurve = computeCurve(angles,thresholdPC,singleProbeInputs, ...
  opts.PedestalWinRate120);
rateRange = sort(double(opts.PedestalWinRateRange(:))');
lowRateCurve = computeCurve(angles,thresholdPC,modelInputs,rateRange(1));
highRateCurve = computeCurve(angles,thresholdPC,modelInputs,rateRange(2));
rateEnvelope = [min(lowRateCurve.totalGain,highRateCurve.totalGain), ...
  max(lowRateCurve.totalGain,highRateCurve.totalGain)];

[bootstrapCI,nBootstrap] = bootstrapModelCI(angles,thresholdPC,modelInputs, ...
  fit,iCD,iCP,iNP,opts.PedestalWinRate120);

observed = observedGainTable(fit,iCP,offsetSummaries,opts.ObservedProbeGains);
comparison = compareObserved(observed,modelInputs,opts.PedestalWinRate120);
downstream = downstreamNoiseDiagnostic(opts.PedestalWinRate120,hitRate);

fprintf('%s [%s]\n',mfilename,codeVersion);
fprintf('Observed IDS offsets loaded for plotting: %s deg\n', ...
  strjoin(compose('%g',observed.angleDeg'),', '));
fprintf('Observed gain minus model prediction at each offset:\n');
disp(comparison(:,{'angleDeg','observedGain','singleProbePredictedGain','predictedGain', ...
  'residual','pointPredictionOutsideObservedCI','source'}));

M = struct();
M.createdAt = datetime('now');
M.createdBy = mfilename;
M.codeVersion = codeVersion;
M.summaryFile = summaryFile;
M.inputs = modelInputs;
M.pedestalWinRate120 = opts.PedestalWinRate120;
M.pedestalWinRateRange = rateRange;
M.mechanismNoiseSDPC = central.mechanismNoiseSDPC;
M.downstreamNoise = downstream;
M.bootstrapN = nBootstrap;
M.observedProbeGains = observed;
M.observedComparison = comparison;
M.offsetSummaries = offsetSummaries;
M.prediction = table(angles,thresholdPC,central.directGain,central.winProbability, ...
  central.winGain,central.totalGain,singleProbeCurve.totalGain, ...
  central.oppositeAxisSelected,bootstrapCI(:,1),bootstrapCI(:,2), ...
  rateEnvelope(:,1),rateEnvelope(:,2), ...
  'VariableNames',{'angleDeg','thresholdPC','directGain','probeWinProbability', ...
  'winGain','totalGain','singleProbeTotalGain','oppositeAxisSelected', ...
  'bootstrapCI95Lower','bootstrapCI95Upper', ...
  'winRateRangeLower','winRateRangeUpper'});

fig = plotModel(M,central,lowRateCurve,highRateCurve);
if opts.SaveResults
  dataFolder = validFolder(fullfile(root,'Data','ProbeGainModel'));
  modelSummary = M; %#ok<NASGU>
  save(fullfile(dataFolder,'IDSGaussianProbeGainModel.mat'),'modelSummary','-v7.3');
  writetable(M.prediction,fullfile(dataFolder,'IDSGaussianProbeGainModel.csv'));
  writetable(M.observedComparison, ...
    fullfile(dataFolder,'IDSGaussianProbeGainModel_ObservedComparison.csv'));
end
if opts.ExportPDF
  plotFolder = validFolder(fullfile(root,'Plots','ProbeGainModel'));
  drawnow;
  exportgraphics(fig,fullfile(plotFolder,'IDSGaussianProbeGainModel.pdf'), ...
    'ContentType','vector','Units','inches','Width',8.5,'Height',5.5, ...
    'Padding',0.1,'PreserveAspectRatio','on');
end
end

function validateSummary(S,path)
required = {'noiseGain','trialTable','psychFit'};
missing = required(~isfield(S,required));
if ~isempty(missing)
  error('runIDSGaussianProbeGainModel:IncompleteSummary', ...
    '%s lacks required field(s): %s.',path,strjoin(missing,', '));
end
if ~isfield(S.noiseGain,'sideRole')
  error('runIDSGaussianProbeGainModel:MissingGainFit', ...
    '%s lacks noiseGain.sideRole.',path);
end
end

function index = predictorIndex(fit,name)
index = find(strcmpi(string(fit.predictorNames),name));
if numel(index) ~= 1
  error('runIDSGaussianProbeGainModel:MissingPredictor', ...
    'Could not identify predictor "%s".',name);
end
end

function thresholdPC = physicalThresholdPC(S)
T = S.trialTable;
scale = median(T.stepCoh./T.alignedCoh,'omitnan');
thresholdPC = S.psychFit.alignedWeibull.threshold*scale;
if ~isfinite(thresholdPC) || thresholdPC <= 0
  error('runIDSGaussianProbeGainModel:BadThreshold', ...
    'Could not recover the physical 120-degree threshold.');
end
end

function A = loadIDSOffsetSummaries(root)
A = struct('gainRows',zeros(0,4),'thresholdRows',zeros(0,2), ...
  'nSessions',zeros(0,1),'files',strings(0,1));
files = dir(fullfile(root,'Data','Probe*','AcrossSessionSummaries', ...
  'IDS_Probe*_AcrossSessionSummary.mat'));
for k = 1:numel(files)
  path = fullfile(files(k).folder,files(k).name);
  Q = load(path,'acrossSummary');
  if ~isfield(Q,'acrossSummary')
    warning('runIDSGaussianProbeGainModel:BadOffsetSummary', ...
      'Skipping %s because it lacks acrossSummary.',path);
    continue
  end
  X = Q.acrossSummary;
  if isfield(X,'probeOffsetDeg') && isfinite(X.probeOffsetDeg)
    offset = double(X.probeOffsetDeg);
  else
    token = regexp(files(k).name,'IDS_Probe(\d+)_','tokens','once');
    if isempty(token)
      warning('runIDSGaussianProbeGainModel:UnknownOffset', ...
        'Skipping %s because its probe offset cannot be identified.',path);
      continue
    end
    offset = str2double(token{1});
  end
  if isfield(X,'noisePredictorType') && ~strcmpi(X.noisePredictorType,'looStep')
    warning('runIDSGaussianProbeGainModel:PredictorMismatch', ...
      'Skipping %s because its predictor is %s rather than looStep.', ...
      path,X.noisePredictorType);
    continue
  end
  if ~isfield(X,'noiseGain') || ~isfield(X.noiseGain,'driftProbe') || ...
      ~isfield(X,'trialTable') || ~isfield(X,'psychFit')
    warning('runIDSGaussianProbeGainModel:IncompleteOffsetSummary', ...
      'Skipping incomplete offset summary %s.',path);
    continue
  end
  F = X.noiseGain.driftProbe;
  iProbe = predictorIndex(F,'Probe');
  gain = F.gain(iProbe);
  if ~isfinite(gain)
    warning('runIDSGaussianProbeGainModel:MissingOffsetGain', ...
      'Skipping %s because its probe gain is unavailable.',path);
    continue
  end
  ci = [nan nan];
  if isfield(F,'CI95') && size(F.CI95,1) >= iProbe && ...
      all(isfinite(F.CI95(iProbe,:)))
    ci = F.CI95(iProbe,:);
  else
    warning('runIDSGaussianProbeGainModel:MissingOffsetCI', ...
      'Displaying the gain from %s without a confidence interval.',path);
  end
  A.gainRows(end+1,:) = [offset gain ci]; %#ok<AGROW>
  try
    threshold = physicalThresholdPC(X);
    A.thresholdRows(end+1,:) = [offset threshold]; %#ok<AGROW>
  catch ME
    warning('runIDSGaussianProbeGainModel:MissingOffsetThreshold', ...
      'Could not use the threshold from %s: %s',path,ME.message);
  end
  if isfield(X,'nSessions')
    A.nSessions(end+1,1) = X.nSessions; %#ok<AGROW>
  else
    A.nSessions(end+1,1) = nan; %#ok<AGROW>
  end
  A.files(end+1,1) = string(path); %#ok<AGROW>
end
if ~isempty(A.gainRows)
  [~,order] = sort(A.gainRows(:,1));
  A.gainRows = A.gainRows(order,:);
  A.nSessions = A.nSessions(order);
  A.files = A.files(order);
  if ~isempty(A.thresholdRows)
    A.thresholdRows = sortrows(A.thresholdRows,1);
  end
else
  warning('runIDSGaussianProbeGainModel:NoOffsetSummaries', ...
    ['No across-offset IDS summaries were found. Run plotIDSSummary ' ...
    'before the model if additional observed offsets should be displayed.']);
end
end

function anchors = combineThresholdAnchors(autoAnchors,threshold120PC,additional)
% Later rows override earlier ones: specialized 120, then explicit inputs.
anchors = [double(autoAnchors); 120 threshold120PC; double(additional)];
anchors = anchors(all(isfinite(anchors),2) & anchors(:,1) >= 0 & ...
  anchors(:,1) <= 180 & anchors(:,2) > 0,:);
if isempty(anchors)
  error('runIDSGaussianProbeGainModel:NoThresholdAnchors', ...
    'At least the 120-degree threshold must be available.');
end
[~,last] = unique(anchors(:,1),'last');
anchors = sortrows(anchors(last,:),1);
end

function threshold = thresholdAtAngles(angles,anchors)
if size(anchors,1) == 1
  threshold = repmat(anchors(1,2),size(angles));
  return
end
threshold = interp1(anchors(:,1),anchors(:,2),angles,'linear','extrap');
threshold(angles < anchors(1,1)) = anchors(1,2);
threshold(angles > anchors(end,1)) = anchors(end,2);
end

function C = computeCurve(angles,thresholdPC,I,pedestalWinRate120)
preferred = gaussianTuning(angles,I.sigmaMT);
opponent = gaussianTuning(180-angles,I.sigmaMT);
preferred120 = gaussianTuning(120,I.sigmaMT);
opponent120 = gaussianTuning(60,I.sigmaMT);
preferredShape = preferred.^I.preferredTuningPower;
preferredShape120 = preferred120.^I.preferredTuningPower;

% Preferred amplitude is fixed by gCD.  The opponent amplitude is the
% value required to reproduce the empirical N-probe gain at 120 degrees.
opponentAmplitude = (I.gCD*preferredShape120-I.direct120)/opponent120;
directGain = I.gCD*preferredShape-opponentAmplitude*opponent;

separation = min(2*angles,360-2*angles);
singleProbeMean = I.pedestalPC*(1+gaussianTuning(separation,I.sigmaMT));
% Each of the two probes lies 180-angle degrees from the opponent axis.
% That axis becomes the better probe competitor as the probes converge
% toward 180 degrees.  Select by mean; this approximation does not model
% trialwise switching between the two competing probe locations.
oppositeAxisProbeMean = 2*I.pedestalPC*gaussianTuning(180-angles,I.sigmaMT);
oppositeAxisSelected = I.includeOppositeAxisProbeMechanism & ...
  oppositeAxisProbeMean > singleProbeMean;
probeMean = singleProbeMean;
probeMean(oppositeAxisSelected) = oppositeAxisProbeMean(oppositeAxisSelected);
driftMean = I.pedestalPC+thresholdPC;
probeMean120 = I.pedestalPC*(1+gaussianTuning(120,I.sigmaMT));
driftMean120 = I.pedestalPC+I.threshold120PC;
z = abs(-sqrt(2)*erfcinv(2*pedestalWinRate120));
mechanismNoiseSDPC = (driftMean120-probeMean120)/(sqrt(2)*z);
if ~isfinite(mechanismNoiseSDPC) || mechanismNoiseSDPC <= 0
  error('runIDSGaussianProbeGainModel:BadMechanismNoise', ...
    'The 120-degree means and win rate imply a nonpositive Gaussian SD.');
end
rawWinProbability = 0.5*erfc((driftMean-probeMean)/(2*mechanismNoiseSDPC));

% A distinct probe mechanism disappears as the probe converges on drift.
% The gate is normalized so the empirical 120-degree win gain is exact.
separateMechanismGate = 1-preferred;
gate120 = 1-preferred120;
winProbability = rawWinProbability.*(separateMechanismGate/gate120);
winGain = I.win120*(winProbability/pedestalWinRate120);
totalGain = directGain+winGain;

C = struct('preferredTuning',preferred,'opponentTuning',opponent, ...
  'opponentAmplitude',opponentAmplitude,'probeSeparationDeg',separation, ...
  'probeMeanPC',probeMean,'singleProbeMeanPC',singleProbeMean, ...
  'oppositeAxisProbeMeanPC',oppositeAxisProbeMean, ...
  'oppositeAxisSelected',oppositeAxisSelected,'driftMeanPC',driftMean, ...
  'mechanismNoiseSDPC',mechanismNoiseSDPC, ...
  'rawWinProbability',rawWinProbability,'winProbability',winProbability, ...
  'directGain',directGain, ...
  'winGain',winGain,'totalGain',totalGain);
end

function y = gaussianTuning(delta,sigma)
y = exp(-0.5*(double(delta)/sigma).^2);
end

function [ci,nBootstrap] = bootstrapModelCI(angles,thresholdPC,I,fit,iCD,iCP,iNP,winRate)
ci = nan(numel(angles),2);
nBootstrap = 0;
if ~isfield(fit,'bootstrapGain') || isempty(fit.bootstrapGain)
  return
end
B = fit.bootstrapGain;
valid = all(isfinite(B(:,[iCD iCP iNP])),2) & B(:,iCP) > B(:,iNP);
B = B(valid,:);
nBootstrap = size(B,1);
if nBootstrap == 0
  return
end
curves = nan(numel(angles),nBootstrap);
for b = 1:nBootstrap
  J = I;
  J.gCD = B(b,iCD);
  J.gCP = B(b,iCP);
  J.gNP = B(b,iNP);
  J.direct120 = J.gNP;
  J.win120 = J.gCP-J.gNP;
  C = computeCurve(angles,thresholdPC,J,winRate);
  curves(:,b) = C.totalGain;
end
for k = 1:numel(angles)
  ci(k,:) = empiricalCI(curves(k,:));
end
end

function ci = empiricalCI(values)
values = sort(values(isfinite(values)));
if isempty(values)
  ci = [nan nan];
  return
end
n = numel(values);
indices = max(1,min(n,round([0.025 0.975]*(n-1)+1)));
ci = values(indices);
end

function T = observedGainTable(fit,iCP,offsetSummaries,additional)
if isfield(fit,'bootstrapCI95') && all(isfinite(fit.bootstrapCI95(iCP,:)))
  ci120 = fit.bootstrapCI95(iCP,:);
else
  ci120 = fit.CI95(iCP,:);
end
autoRows = offsetSummaries.gainRows;
values = [autoRows; 120 fit.gain(iCP) ci120; double(additional)];
source = [repmat("IDS offset summary",size(autoRows,1),1); ...
  "IDS 120 C/N summary"; repmat("Explicit override",size(additional,1),1)];
nSessions = [offsetSummaries.nSessions; fit.nSessions; ...
  nan(size(additional,1),1)];
% Later rows override earlier rows at the same angle.
[~,last] = unique(values(:,1),'last');
values = values(last,:);
source = source(last);
nSessions = nSessions(last);
[~,order] = sort(values(:,1));
values = values(order,:);
source = source(order);
nSessions = nSessions(order);
T = table(values(:,1),values(:,2),values(:,3),values(:,4),source,nSessions, ...
  'VariableNames',{'angleDeg','gain','CI95Lower','CI95Upper','source','nSessions'});
end

function T = compareObserved(observed,I,winRate)
% Evaluate at the measured offsets, even when the displayed curve uses a
% coarser AnglesDeg grid.  The 120-degree point is a calibration anchor.
angles = observed.angleDeg;
threshold = thresholdAtAngles(angles,I.thresholdAnchors);
C = computeCurve(angles,threshold,I,winRate);
singleProbeInputs = I;
singleProbeInputs.includeOppositeAxisProbeMechanism = false;
singleProbeCurve = computeCurve(angles,threshold,singleProbeInputs,winRate);
hasCI = isfinite(observed.CI95Lower) & isfinite(observed.CI95Upper);
outside = hasCI & (C.totalGain < observed.CI95Lower | ...
  C.totalGain > observed.CI95Upper);
T = table(angles,observed.gain,singleProbeCurve.totalGain, ...
  C.totalGain,observed.gain-C.totalGain, ...
  C.directGain,C.winGain,C.winProbability, ...
  C.singleProbeMeanPC,C.oppositeAxisProbeMeanPC,C.oppositeAxisSelected, ...
  observed.CI95Lower,observed.CI95Upper,hasCI,outside, ...
  observed.source,observed.nSessions,angles == 120, ...
  'VariableNames',{'angleDeg','observedGain','singleProbePredictedGain', ...
  'predictedGain','residual', ...
  'predictedDirectGain','predictedWinGain','predictedWinProbability', ...
  'singleProbeMeanPC','oppositeAxisProbeMeanPC','oppositeAxisSelected', ...
  'observedCI95Lower','observedCI95Upper','hasCI','pointPredictionOutsideObservedCI', ...
  'source','nSessions','isCalibrationAnchor'});
end

function D = downstreamNoiseDiagnostic(perPedestalWinRate,hitRate)
noChangeWinRate = 3*perPedestalWinRate;
preDecisionCorrect = 1-noChangeWinRate;
denominator = 2*preDecisionCorrect-1;
flipRate = (preDecisionCorrect-hitRate)/denominator;
valid = isfinite(flipRate) && flipRate >= 0 && flipRate <= 0.5;
D = struct('perPedestalWinRate',perPedestalWinRate, ...
  'noChangeWinRate',noChangeWinRate, ...
  'preDecisionCorrect',preDecisionCorrect, ...
  'observedHitRate',hitRate,'symmetricFlipRate',flipRate,'valid',valid);
end

function fig = plotModel(M,C,lowRateCurve,highRateCurve)
fig = figure('Color','w','Units','inches','Position',[1 1 8.5 5.5]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
title(tl,'IDS Gaussian Probe-Gain Model','FontWeight','bold');
I = M.inputs;
subtitle(tl,sprintf('120° anchors: drift %.3f, N probe %.3f, probe win %.3f; MT SD %.1f°, preferred power %.1f', ...
  I.gCD,I.gNP,I.win120,I.sigmaMT,I.preferredTuningPower));

ax = nexttile(tl,1); hold(ax,'on');
x = M.prediction.angleDeg;
ci = [M.prediction.bootstrapCI95Lower M.prediction.bootstrapCI95Upper];
if all(isfinite(ci),'all')
  fill(ax,[x; flipud(x)],[ci(:,1); flipud(ci(:,2))],[.75 .82 .92], ...
    'EdgeColor','none','FaceAlpha',.55,'DisplayName','120° gain bootstrap 95% CI');
end
plot(ax,x,C.directGain,'-','Color',[.85 .35 .15],'LineWidth',1.5, ...
  'DisplayName','Direct MT term');
plot(ax,x,C.winGain,'-','Color',[.15 .45 .8],'LineWidth',1.5, ...
  'DisplayName','Probe-win term');
plot(ax,x,C.totalGain,'k-','LineWidth',2,'DisplayName','Total prediction');
if I.includeOppositeAxisProbeMechanism
  plot(ax,x,M.prediction.singleProbeTotalGain,':','Color',[.35 .35 .35], ...
    'LineWidth',1.5,'DisplayName','Probe-centered competitor only');
end
plot(ax,x,lowRateCurve.totalGain,'k--','LineWidth',.8,'HandleVisibility','off');
plot(ax,x,highRateCurve.totalGain,'k--','LineWidth',.8, ...
  'DisplayName','Pedestal-win-rate range');
O = M.observedProbeGains;
for k = 1:height(O)
  plot(ax,[O.angleDeg(k) O.angleDeg(k)],[O.CI95Lower(k) O.CI95Upper(k)], ...
    'ko-','LineWidth',1.2,'HandleVisibility','off');
end
plot(ax,O.angleDeg,O.gain,'ko','MarkerFaceColor','k','MarkerSize',5, ...
  'DisplayName','Observed');
yline(ax,0,'k:','HandleVisibility','off'); xline(ax,120,':','HandleVisibility','off');
xline(ax,170,':','HandleVisibility','off');
xlabel(ax,'Probe offset (deg)'); ylabel(ax,'Probe gain');
xlim(ax,[0 180]); xticks(ax,0:30:180); grid(ax,'on'); box(ax,'off');
legend(ax,'Location','best','FontSize',7);

ax = nexttile(tl,2); hold(ax,'on');
yyaxis(ax,'left');
hDrift = plot(ax,x,C.driftMeanPC,'k-','LineWidth',1.5,'DisplayName','Drift mean');
hProbe = plot(ax,x,C.probeMeanPC,'-','Color',[.2 .6 .75],'LineWidth',1.5, ...
  'DisplayName','Effective probe mean');
ylabel(ax,'Mechanism mean (% coherence equivalent)');
yyaxis(ax,'right');
hWin = plot(ax,x,C.winProbability,'-','Color',[.55 .25 .65],'LineWidth',1.5, ...
  'DisplayName','Probe-win probability');
ylabel(ax,'Probe-win probability'); ylim(ax,[0 max(.35,1.1*max(C.winProbability))]);
xlabel(ax,'Probe offset (deg)'); xlim(ax,[0 180]); xticks(ax,0:30:180);
grid(ax,'on'); box(ax,'off');
legend(ax,[hDrift hProbe hWin],'Location','southoutside','FontSize',7);

end
