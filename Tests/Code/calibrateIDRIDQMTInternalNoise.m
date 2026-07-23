function calibration=calibrateIDRIDQMTInternalNoise(context,varargin)
% calibrateIDRIDQMTInternalNoise  Match session 75% no-noise thresholds.
%
% A fixed dense bank is used for every dataset and session. Independent
% white Gaussian MT-population noise is drawn for the two patches and then
% projected through the readout bank, inducing correlated candidate noise.

context=validateAnalysisContext(context);
if context.Mode~="synthetic"
  error('IDReadout:SyntheticContextRequired', ...
    'Internal-noise calibration requires a synthetic context.');
end
validateSyntheticManifest(context);

p=inputParser;
addParameter(p,'IDQSummaryPath',defaultSummaryPath('IDQ','IDQ_AcrossSideSummary.mat'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'IDRSummaryPath',defaultSummaryPath('IDR','IDR_SideGainAnalysis.mat'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'SigmaMTDeg',37.5,@isPositiveScalar);
addParameter(p,'SigmaReadoutDeg',20,@isPositiveScalar);
addParameter(p,'CandidateSpacingDeg',2,@isPositiveScalar);
addParameter(p,'TargetPerformance',0.75, ...
  @(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>.5&&x<1);
addParameter(p,'NMonteCarlo',2000, ...
  @(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=200&&fix(x)==x);
addParameter(p,'RectifyCandidates',true,@(x) islogical(x)&&isscalar(x));
addParameter(p,'BisectionIterations',18, ...
  @(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=8&&fix(x)==x);
addParameter(p,'Visible','off',@(x) any(strcmpi(string(x),["on","off"])));
addParameter(p,'SaveOutputs',true,@(x) islogical(x)&&isscalar(x));
parse(p,varargin{:});opts=p.Results;

sourcePaths=string({opts.IDQSummaryPath;opts.IDRSummaryPath});
sourceRecords=validateSimulationSourceFiles(context,sourcePaths);
sessions=makeIDRIDQMTCalibrationSessions(sourcePaths(1),sourcePaths(2));
mtModel=makeMTReadoutForwardModel('sigmaMTDeg',opts.SigmaMTDeg);
directions=(0:opts.CandidateSpacingDeg:(360-opts.CandidateSpacingDeg))';
if abs(numel(directions)*opts.CandidateSpacingDeg-360)>1e-9
  error('calibrateIDRIDQMTInternalNoise:SpacingMustTileCircle', ...
    'CandidateSpacingDeg must divide 360 exactly.');
end
bank=makeGaussianReadoutBank(directions,opts.SigmaReadoutDeg,mtModel);
noiseModel=factorGaussianReadoutNoise(bank);

stream=RandStream('mt19937ar','Seed',context.Seed);
zChange=randn(stream,opts.NMonteCarlo,noiseModel.rank);
zNoChange=randn(stream,opts.NMonteCarlo,noiseModel.rank);
unitChange=zChange*noiseModel.factor';
unitNoChange=zNoChange*noiseModel.factor';
preferredKernel=mtPopulationTemplate(0,mtModel)*double(bank.weightsPhi)';

n=height(sessions);sigmaMTNoise=nan(n,1);achievedPerformance=nan(n,1);
lowerBound=nan(n,1);upperBound=nan(n,1);
for k=1:n
  base=sessions.baseCohPC(k);threshold=sessions.thresholdPC(k);
  changeSignal=(base+threshold).*preferredKernel;
  noChangeSignal=base.*preferredKernel;
  objective=@(sigma) simulatedPerformance(changeSignal,noChangeSignal, ...
    unitChange,unitNoChange,sigma,opts.RectifyCandidates);
  [sigmaMTNoise(k),achievedPerformance(k),lowerBound(k),upperBound(k)]= ...
    solveNoiseScale(objective,threshold,preferredKernel, ...
      opts.TargetPerformance,opts.BisectionIterations);
end
sessions.sigmaMTNoise=sigmaMTNoise;
sessions.sigmaMTNoisePerThreshold=sigmaMTNoise./sessions.thresholdPC;
sessions.achievedPerformance=achievedPerformance;
sessions.finalLowerBound=lowerBound;
sessions.finalUpperBound=upperBound;

calibration=struct();
calibration.createdAt=datetime('now');calibration.createdBy=mfilename;
calibration.dataOrigin="synthetic_calibration_from_manifested_thresholds";
calibration.context=context;calibration.sourceRecords=sourceRecords;
calibration.sigmaMTDeg=opts.SigmaMTDeg;
calibration.sigmaReadoutDeg=opts.SigmaReadoutDeg;
calibration.candidateSpacingDeg=opts.CandidateSpacingDeg;
calibration.candidateDirectionsDeg=bank.candidateDirectionsDeg;
calibration.nCandidates=numel(bank.candidateDirectionsDeg);
calibration.noiseRank=noiseModel.rank;
calibration.targetPerformance=opts.TargetPerformance;
calibration.nMonteCarlo=opts.NMonteCarlo;
calibration.rectifyCandidates=opts.RectifyCandidates;
calibration.sessionTable=sessions;
calibration.noiseDefinition=noiseModel.definition;
calibration.calibrationDefinition=[ ...
  "Per-session white-MT noise SD chosen by common-random-number bisection " + ...
  "so dense-bank hard-max patch comparison is 75% correct at the measured " + ...
  "session no-noise threshold and median prestep coherence"];

fig=makePlot(calibration,opts.Visible);
if opts.SaveOutputs
  dataFolder=analysisPath(context,'Common Code','Data','Simulation');
  plotFolder=analysisPath(context,'Common Code','Plots','Simulation');
  if ~isfolder(dataFolder),mkdir(dataFolder);end
  if ~isfolder(plotFolder),mkdir(plotFolder);end
  matPath=fullfile(dataFolder,'IDRIDQ_MTInternalNoiseCalibration.mat');
  csvPath=fullfile(dataFolder,'IDRIDQ_MTInternalNoiseCalibration.csv');
  pdfPath=fullfile(plotFolder,'IDRIDQ_MTInternalNoiseCalibration.pdf');
  calibration.outputPaths=string({matPath;csvPath;pdfPath});
  info=whos('calibration');assertSimulationRunCapacity(context,info.bytes);
  save(matPath,'calibration','-v7.3');writetable(sessions,csvPath);
  exportgraphics(fig,pdfPath,'ContentType','vector');
  assertSimulationRunCapacity(context,0);
end
if strcmpi(opts.Visible,'off'),close(fig);end
end

function p=simulatedPerformance(changeSignal,noChangeSignal, ...
  unitChange,unitNoChange,sigma,rectify)
A=changeSignal+sigma.*unitChange;
B=noChangeSignal+sigma.*unitNoChange;
if rectify,A=max(A,0);B=max(B,0);end
p=mean(max(A,[],2)>max(B,[],2));
end

function [sigma,p,lo,hi]=solveNoiseScale(objective,threshold,kernel,target,nIter)
preferred=max(kernel);
z=sqrt(2)*erfinv(2*target-1);
analytic=max(preferred*threshold/(sqrt(2)*z),eps);
lo=0;hi=analytic;
while objective(hi)>target && hi<analytic*2^20,hi=hi*2;end
if objective(hi)>target
  error('calibrateIDRIDQMTInternalNoise:CouldNotBracket', ...
    'Could not bracket the target performance.');
end
for j=1:nIter
  mid=(lo+hi)/2;
  if objective(mid)>target,lo=mid;else,hi=mid;end
end
sigma=(lo+hi)/2;p=objective(sigma);
end

function fig=makePlot(C,visible)
T=C.sessionTable;
fig=figure('Color','w','Visible',visible,'Units','inches','Position',[1 1 10 4.5]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('Dense-bank MT internal-noise calibration (%d candidates)',C.nCandidates));
ax=nexttile(tl);hold(ax,'on');
for dataset=["IDQ","IDR"]
  use=T.dataset==dataset;
  scatter(ax,T.thresholdPC(use),T.sigmaMTNoise(use),20,'filled', ...
    'DisplayName',dataset);
end
xlabel(ax,'Session 75% threshold (% coherence)');ylabel(ax,'Calibrated white-MT noise SD');
legend(ax,'Location','best');grid(ax,'on');box(ax,'off');
ax=nexttile(tl);hold(ax,'on');
for dataset=["IDQ","IDR"]
  use=T.dataset==dataset;
  scatter(ax,T.baseToThresholdRatio(use),T.achievedPerformance(use),20,'filled', ...
    'DisplayName',dataset);
end
yline(ax,C.targetPerformance,'k--','Target','HandleVisibility','off');
xlabel(ax,'Prestep coherence / threshold');ylabel(ax,'Monte Carlo performance');
ylim(ax,[.70 .80]);legend(ax,'Location','best');grid(ax,'on');box(ax,'off');
end

function path=defaultSummaryPath(domain,file)
path=fullfile(domainFolder(mfilename('fullpath'),domain), ...
  'Data','AcrossSessionSummaries',file);
end
function tf=isPositiveScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0;
end
