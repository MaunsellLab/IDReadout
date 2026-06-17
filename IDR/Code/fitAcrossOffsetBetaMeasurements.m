function betaSummary = fitAcrossOffsetBetaMeasurements(varargin)
% fitAcrossOffsetBetaMeasurements  Complete across-offset beta readout fit.
%
% Loads current per-probe regression products, estimates one shared
% probe/preferred scale at each offset while allowing session-specific
% intercepts and preferred sensitivities, bootstraps those scale estimates,
% and calls fitAcrossOffsetReadout.
%
% Name-value options:
%   'StepType'          'inc' (default), 'dec', or 'combined'
%   'NBoot'             hierarchical bootstrap count (default 1000)
%   'RandomSeed'        bootstrap seed (default 1)
%   'Bin179With180'     combine 179 and 180 deg (default false)
%   'MakePlots'         default true
%   'SaveFile'          default Data/AcrossOffsetSummaries/IDR_acrossOffsetBetaSummary.mat
%   'PlotDir'           default Plots/AcrossProbes/ReadoutFits/Beta
%   'Bounds'            forwarded to fitAcrossOffsetReadout
%   'Verbose'           default true

baseFolder = domainFolder(mfilename('fullpath'));
defaultSave = fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries', 'IDR_acrossOffsetBetaSummary.mat');
defaultPlotDir = fullfile(baseFolder, 'Plots', 'AcrossProbes', ...
  'ReadoutFits', 'Beta');

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'StepType', 'inc', @(x) any(strcmpi(string(x),["inc","dec","combined"])));
addParameter(p, 'NBoot', 10, @(x) isnumeric(x) && isscalar(x) && x >= 0 && mod(x,1)==0);
addParameter(p, 'RandomSeed', 1, @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'Bin179With180', true, @(x) islogical(x) && isscalar(x));
% addParameter(p, 'MakePlots', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'SaveFile', defaultSave, @(x) ischar(x) || isstring(x));
addParameter(p, 'PlotDir', defaultPlotDir, @(x) ischar(x) || isstring(x));
addParameter(p, 'PlotOnly', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Bounds', struct(), @isstruct);
addParameter(p, 'Verbose', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opts = p.Results;

if opts.PlotOnly
  S = load(opts.SaveFile, 'betaSummary');
  betaSummary = S.betaSummary;
  if ~isfolder(opts.PlotDir)
    mkdir(opts.PlotDir);
  end
  plotBetaRatiosByOffset(betaSummary, fullfile(opts.PlotDir, 'BetaRatiosByOffset.pdf'));
  plotBetaReadoutFit(betaSummary,  fullfile(opts.PlotDir, 'BetaReadoutFit.pdf'));
  return;
end

opts.StepType = lower(char(string(opts.StepType)));
opts.SaveFile = char(opts.SaveFile);
opts.PlotDir = char(opts.PlotDir);
if ~isempty(opts.RandomSeed), rng(opts.RandomSeed); end

files = dir(fullfile(baseFolder, 'Data', 'Probe*', 'Regression', '*_scalarNoiseRegression.mat'));
if isempty(files)
  error('fitAcrossOffsetBetaMeasurements:NoRegressionFiles', ...
    'No regression files were found.');
end

sessionRecords = repmat(emptySessionRecord(),0,1);
for i = 1:numel(files)
  filePath = fullfile(files(i).folder, files(i).name);
  S = load(filePath, 'reg');
  if ~isfield(S,'reg') || ~isfield(S.reg,'analysisName') || ...
      ~strcmp(S.reg.analysisName,'kernelWeightedProbeRegression') || ...
      ~isfield(S.reg,'fitByStep') || ~isfield(S.reg.fitByStep,opts.StepType)
    continue
  end
  F = S.reg.fitByStep.(opts.StepType);
  if ~isfield(F,'fitUsable') || ~F.fitUsable || ...
      ~isfield(S.reg,'trialTable') || isempty(S.reg.trialTable)
    continue
  end
  T = S.reg.trialTable;
  switch opts.StepType
    case 'inc', use = T.changeIndex == 2;
    case 'dec', use = T.changeIndex == 1;
    otherwise, use = true(height(T),1);
  end
  if sum(use) < 3 || numel(unique(T.correct(use))) < 2
    continue
  end

  r = emptySessionRecord();
  r.fileName = files(i).name;
  r.filePath = filePath;
  r.probeOffsetDeg = headerScalar(S.reg.sessionProbeHeader,'probeDirDeg');
  r.xPref = double(T.effectivePrefNoisePC(use));
  r.xProbe = double(T.effectiveProbeNoisePC(use));
  r.correct = double(T.correct(use));
  r.betaPref = fieldOrNaN(F,'betaPref');
  r.betaProbe = fieldOrNaN(F,'betaProbe');
  r.betaRatio = fieldOrNaN(F,'betaRatio');
  r.betaPrefSE = fieldOrNaN(F,'betaPrefSE');
  r.betaProbeSE = fieldOrNaN(F,'betaProbeSE');
  r.betaRatioSE = fieldOrNaN(F,'betaRatioSE');
  sessionRecords(end+1,1) = r; %#ok<AGROW>
end

if isempty(sessionRecords)
  error('fitAcrossOffsetBetaMeasurements:NoUsableFiles', ...
    'No current usable regression products were found.');
end

% Mirror the current kernel-readout eligibility: paired probes >1 and <=179.
offsets = [sessionRecords.probeOffsetDeg];
eligible = isfinite(offsets) & offsets > 1 & offsets <= 179;
if opts.Bin179With180
  eligible = isfinite(offsets) & offsets > 1 & offsets <= 180;
  for i = 1:numel(sessionRecords)
    if sessionRecords(i).probeOffsetDeg == 180
      sessionRecords(i).probeOffsetDeg = 179;
    end
  end
end
sessionRecords = sessionRecords(eligible);

offsetKeys = unique([sessionRecords.probeOffsetDeg]);
offsetFits = repmat(emptyOffsetFit(),numel(offsetKeys),1);
bootScaleMat = nan(opts.NBoot,numel(offsetKeys));

for k = 1:numel(offsetKeys)
  thisOffset = offsetKeys(k);
  R = sessionRecords([sessionRecords.probeOffsetDeg] == thisOffset);
  D = recordsToSessionData(R);
  fit = fitSharedProbeScale(D);

  offsetFits(k).probeOffsetDeg = thisOffset;
  offsetFits(k).nSessions = numel(R);
  offsetFits(k).nTrials = sum(cellfun(@(x) numel(x.correct),D));
  offsetFits(k).scale = fit.scale;
  offsetFits(k).scaleHessianSE = fit.scaleSE;
  offsetFits(k).scaleHessianCI95 = fit.scaleCI95;
  offsetFits(k).fit = fit;
  offsetFits(k).sessionFileNames = {R.fileName};
  offsetFits(k).sessionBetaPref = [R.betaPref];
  offsetFits(k).sessionBetaProbe = [R.betaProbe];
  offsetFits(k).sessionBetaRatio = [R.betaRatio];
  offsetFits(k).sessionBetaPrefSE = [R.betaPrefSE];
  offsetFits(k).sessionBetaProbeSE = [R.betaProbeSE];
  offsetFits(k).sessionBetaRatioSE = [R.betaRatioSE];

  if opts.Verbose
    fprintf('Offset %g: %d sessions, %d trials, pooled scale %.6g\n', ...
      thisOffset, offsetFits(k).nSessions, offsetFits(k).nTrials, fit.scale);
  end

  for b = 1:opts.NBoot
    nS = numel(R);
    sampled = randi(nS,nS,1);
    Db = cell(nS,1);
    for j = 1:nS
      src = R(sampled(j));
      nT = numel(src.correct);
      idx = randi(nT,nT,1);
      Db{j} = struct('xPref',src.xPref(idx), ...
        'xProbe',src.xProbe(idx),'correct',src.correct(idx));
    end
    try
      Fb = fitSharedProbeScale(Db);
      if Fb.fitUsable
        bootScaleMat(b,k) = Fb.scale;
      end
    catch
      bootScaleMat(b,k) = NaN;
    end
  end
end

measurements = [offsetFits.scale];
variances = nan(size(measurements));
for k = 1:numel(offsetFits)
  x = bootScaleMat(:,k);
  x = x(isfinite(x));
  if numel(x) >= 2
    variances(k) = var(x,0);
    offsetFits(k).bootMean = mean(x);
    offsetFits(k).bootMedian = median(x);
    offsetFits(k).boot68 = prctile(x,[16 84]);
    offsetFits(k).boot95 = prctile(x,[2.5 97.5]);
    offsetFits(k).nValidBoot = numel(x);
  elseif opts.NBoot == 0
    variances(k) = offsetFits(k).scaleHessianSE^2;
    offsetFits(k).bootMean = NaN;
    offsetFits(k).bootMedian = NaN;
    offsetFits(k).boot68 = [NaN NaN];
    offsetFits(k).boot95 = [NaN NaN];
    offsetFits(k).nValidBoot = 0;
  end
end

% Protect the shared readout fit against exactly zero numerical variances.
finiteVar = variances(isfinite(variances) & variances > 0);
if isempty(finiteVar)
  error('fitAcrossOffsetBetaMeasurements:NoUsableVariance', ...
    'No usable offset variance estimates were obtained.');
end
varFloor = max(1e-8,0.01*median(finiteVar));
variances(isfinite(variances) & variances <= 0) = varFloor;

readoutFitSummary = fitAcrossOffsetReadout( ...
  measurements, variances, offsetKeys, ...
  'NSessions', [offsetFits.nSessions], ...
  'Bounds', opts.Bounds, ...
  'SourceMeasureType', 'pooledBetaScale', ...
  'SourceSideType', 'change', ...
  'SourceStepType', opts.StepType, ...
  'SourceMode', 'sharedScale_sessionSpecificPreferredBeta');

betaSummary = struct();
betaSummary.version = 1;
betaSummary.analysisName = 'acrossOffsetBetaMeasurements';
betaSummary.meta = struct( ...
  'createdDate', datetime('now'), 'stepType', opts.StepType, 'nBoot', opts.NBoot,  'randomSeed', opts.RandomSeed, ...
  'bin179With180', opts.Bin179With180, 'model', ['session-specific intercept and preferred beta; ' ...
    'shared probe/preferred scale within offset']);
betaSummary.sessionRecords = sessionRecords;
betaSummary.offsetFits = offsetFits;
betaSummary.measurements = struct( ...
  'offsetsDeg', offsetKeys, ...
  'pooledScale', measurements, ...
  'bootstrapVar', variances, ...
  'nSessions', [offsetFits.nSessions], ...
  'nTrials', [offsetFits.nTrials]);
betaSummary.bootstrap = struct( ...
  'bootScaleMat', bootScaleMat, ...
  'varFloor', varFloor);
betaSummary.readoutFitSummary = readoutFitSummary;
betaSummary.readoutModels = readoutFitSummary.readoutModels;
betaSummary.readoutModel = readoutFitSummary.readoutModel;
betaSummary.readoutModelComparison = readoutFitSummary.readoutModelComparison;

saveDir = fileparts(opts.SaveFile);
if ~isfolder(saveDir), mkdir(saveDir); end
save(opts.SaveFile,'betaSummary','-v7.3');
if opts.Verbose, fprintf('Saved %s\n',opts.SaveFile); end

% if opts.MakePlots
%   if ~isfolder(opts.PlotDir), mkdir(opts.PlotDir); end
%   plotBetaRatiosByOffset(betaSummary, fullfile(opts.PlotDir,'BetaRatiosByOffset.pdf'));
%   plotBetaReadoutFit(betaSummary, fullfile(opts.PlotDir,'BetaReadoutFit.pdf'));
% end
end

function D = recordsToSessionData(R)
D = cell(numel(R),1);
for i = 1:numel(R)
  D{i} = struct('xPref',R(i).xPref,'xProbe',R(i).xProbe,'correct',R(i).correct);
end
end

% function plotBetaRatiosByOffset(S, savePath)
% F = S.offsetFits;
% fig = figure('Color','w','Position',[100 100 1050 500]); hold on;
% hSession = gobjects(1); hPool = gobjects(1);
% for k = 1:numel(F)
%   n = numel(F(k).sessionBetaRatio);
%   jitter = zeros(n,1);
%   if n > 1, jitter = linspace(-1.5,1.5,n)'; end
%   h = errorbar(F(k).probeOffsetDeg+jitter, F(k).sessionBetaRatio(:), ...
%     F(k).sessionBetaRatioSE(:), 'o','LineStyle','none','MarkerSize',4,'CapSize',0);
%   if k==1, hSession = h; end
%   if all(isfinite(F(k).boot95))
%     lo = F(k).scale-F(k).boot95(1); hi = F(k).boot95(2)-F(k).scale;
%   else
%     lo = 1.96*F(k).scaleHessianSE; hi = lo;
%   end
%   hp = errorbar(F(k).probeOffsetDeg,F(k).scale,lo,hi,'ks', ...
%     'MarkerFaceColor','k','MarkerSize',8,'LineWidth',1.4,'CapSize',8);
%   if k==1, hPool = hp; end
% end
% yline(0,':'); yline(1,'--');
% xlabel('Probe direction offset (deg)');
% ylabel('\beta_{probe}/\beta_{pref}');
% title(sprintf('%s session ratios and pooled shared scales',upper(S.meta.stepType)));
% xticks([F.probeOffsetDeg]);
% legend([hSession hPool],{'Session ratio \pm SE','Pooled shared scale (95% CI)'},'Location','best');
% box off;
% exportgraphics(fig, savePath,'ContentType','vector');
% end
% 
% function plotBetaReadoutFit(S, savePath)
% M = S.measurements;
% R = S.readoutFitSummary.readoutModels;
% fig = figure('Color','w','Position',[100 100 850 520]); hold on;
% ci = nan(numel(S.offsetFits),2);
% for k=1:numel(S.offsetFits), ci(k,:)=S.offsetFits(k).boot95; end
% lo = M.pooledScale-ci(:,1)'; hi = ci(:,2)'-M.pooledScale;
% missing = ~isfinite(lo) | ~isfinite(hi);
% lo(missing)=1.96*sqrt(M.bootstrapVar(missing));
% hi(missing)=lo(missing);
% hObs = errorbar(M.offsetsDeg, M.pooledScale,lo,hi,'ko', 'MarkerFaceColor','k','LineWidth',1.2,'CapSize',8);
% hh = hObs; labels = {'Pooled beta scale (95% CI)'};
% if isfield(R.signedDOG,'fit') && ~isempty(R.signedDOG.fit) && R.signedDOG.fit.fitUsable
%   h=plot(R.signedDOG.plotOffsetsDeg,R.signedDOG.plotPredictedScale,'-','LineWidth',1.5);
%   hh(end+1)=h; labels{end+1}='Signed DOG';
% end
% if isfield(R.rectifiedDOG,'fit') && ~isempty(R.rectifiedDOG.fit) && R.rectifiedDOG.fit.fitUsable
%   h=plot(R.rectifiedDOG.plotOffsetsDeg,R.rectifiedDOG.plotPredictedScale,'-','LineWidth',1.5);
%   hh(end+1)=h; labels{end+1}='Rectified DOG';
% end
% yline(0,':');
% xlabel('Probe direction offset (deg)'); ylabel('Normalized beta scale');
% title(sprintf('%s pooled beta scales and MT/readout fit',upper(S.meta.stepType)));
% xlim([0 180]); legend(hh,labels,'Location','best'); box off;
% exportgraphics(fig,savePath,'ContentType','vector');
% end

function r = emptySessionRecord()
r = struct('fileName','','filePath','','probeOffsetDeg',NaN, ...
  'xPref',[],'xProbe',[],'correct',[], ...
  'betaPref',NaN,'betaProbe',NaN,'betaRatio',NaN, ...
  'betaPrefSE',NaN,'betaProbeSE',NaN,'betaRatioSE',NaN);
end

function f = emptyOffsetFit()
f = struct('probeOffsetDeg',NaN,'nSessions',0,'nTrials',0, ...
  'scale',NaN,'scaleHessianSE',NaN,'scaleHessianCI95',[NaN NaN], ...
  'fit',struct(),'sessionFileNames',{{}}, ...
  'sessionBetaPref',[],'sessionBetaProbe',[],'sessionBetaRatio',[], ...
  'sessionBetaPrefSE',[],'sessionBetaProbeSE',[],'sessionBetaRatioSE',[], ...
  'bootMean',NaN,'bootMedian',NaN,'boot68',[NaN NaN], ...
  'boot95',[NaN NaN],'nValidBoot',0);
end

function v = fieldOrNaN(S,name)
if isfield(S,name), v=double(S.(name)); else, v=NaN; end
end

function v = headerScalar(H,name)
if ~isfield(H,name), error('Missing header field %s.',name); end
v=H.(name); if isstruct(v)&&isfield(v,'data'),v=v.data;end
v=double(v(1));
end
