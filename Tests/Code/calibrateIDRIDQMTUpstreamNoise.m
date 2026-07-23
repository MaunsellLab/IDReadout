function calibration=calibrateIDRIDQMTUpstreamNoise(context,varargin)
% calibrateIDRIDQMTUpstreamNoise  Calibrate pre-MT directional variability.

context=validateAnalysisContext(context);
if context.Mode~="synthetic"
  error('IDReadout:SyntheticContextRequired', ...
    'Upstream-noise calibration requires a synthetic context.');
end
validateSyntheticManifest(context);

p=inputParser;
addParameter(p,'IDQSummaryPath',defaultSummaryPath('IDQ','IDQ_AcrossSideSummary.mat'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'IDRSummaryPath',defaultSummaryPath('IDR','IDR_SideGainAnalysis.mat'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'SigmaMTDeg',37.5,@isPositiveScalar);
addParameter(p,'SigmaReadoutDeg',20,@isPositiveScalar);
addParameter(p,'NegativeAsymptoteMagnitude',0,@isNonnegativeScalar);
addParameter(p,'SurroundMagnitude',0,@isNonnegativeScalar);
addParameter(p,'SigmaSurroundDeg',nan,@isOptionalPositiveScalar);
addParameter(p,'CandidateSpacingDeg',2,@isPositiveScalar);
addParameter(p,'InputSpacingDeg',1,@isPositiveScalar);
addParameter(p,'IDQPreStepPedestalPC',7, ...
  @(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
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
requireCircleTiling(opts.CandidateSpacingDeg,'CandidateSpacingDeg');
requireCircleTiling(opts.InputSpacingDeg,'InputSpacingDeg');

sourcePaths=string({opts.IDQSummaryPath;opts.IDRSummaryPath});
sourceRecords=validateSimulationSourceFiles(context,sourcePaths);
sessions=makeIDRIDQMTUpstreamCalibrationSessions( ...
  sourcePaths(1),sourcePaths(2),opts.IDQPreStepPedestalPC);
mtModel=makeMTReadoutForwardModel('sigmaMTDeg',opts.SigmaMTDeg);
candidateDirections=(0:opts.CandidateSpacingDeg: ...
  (360-opts.CandidateSpacingDeg))';
inputDirections=(0:opts.InputSpacingDeg:(360-opts.InputSpacingDeg))';
bank=makeGaussianReadoutBank(candidateDirections,opts.SigmaReadoutDeg,mtModel, ...
  'NegativeAsymptoteMagnitude',opts.NegativeAsymptoteMagnitude, ...
  'SurroundMagnitude',opts.SurroundMagnitude, ...
  'SigmaSurroundDeg',opts.SigmaSurroundDeg);
noiseModel=factorUpstreamDirectionalNoise(inputDirections,bank,mtModel);

stream=RandStream('mt19937ar','Seed',context.Seed);
zChange=randn(stream,opts.NMonteCarlo,noiseModel.rank);
zNoChange=randn(stream,opts.NMonteCarlo,noiseModel.rank);
unitChange=zChange*noiseModel.factor';
unitNoChange=zNoChange*noiseModel.factor';
preferredKernel=mtPopulationTemplate(0,mtModel)*double(bank.weightsPhi)';

n=height(sessions);upstreamDirectionNoiseSD=nan(n,1);
achievedPerformance=nan(n,1);lowerBound=nan(n,1);upperBound=nan(n,1);
for k=1:n
  base=sessions.decisionBaseCohPC(k);threshold=sessions.thresholdPC(k);
  changeSignal=(base+threshold).*preferredKernel;
  noChangeSignal=base.*preferredKernel;
  objective=@(sigma) simulatedPerformance(changeSignal,noChangeSignal, ...
    unitChange,unitNoChange,sigma,opts.RectifyCandidates);
  [upstreamDirectionNoiseSD(k),achievedPerformance(k), ...
    lowerBound(k),upperBound(k)]=solveNoiseScale( ...
      objective,threshold,preferredKernel,opts.TargetPerformance, ...
      opts.BisectionIterations);
end
sessions.upstreamDirectionNoiseSD=upstreamDirectionNoiseSD;
sessions.upstreamNoisePerThreshold= ...
  upstreamDirectionNoiseSD./sessions.thresholdPC;
sessions.achievedPerformance=achievedPerformance;
sessions.finalLowerBound=lowerBound;sessions.finalUpperBound=upperBound;

correlationTable=makeMTCorrelationTable(noiseModel,mtModel);
calibration=struct();calibration.createdAt=datetime('now');
calibration.createdBy=mfilename;
calibration.dataOrigin="synthetic_calibration_from_manifested_thresholds";
calibration.context=context;calibration.sourceRecords=sourceRecords;
calibration.sigmaMTDeg=opts.SigmaMTDeg;
calibration.sigmaReadoutDeg=opts.SigmaReadoutDeg;
calibration.negativeAsymptoteMagnitude=opts.NegativeAsymptoteMagnitude;
calibration.surroundMagnitude=opts.SurroundMagnitude;
calibration.sigmaSurroundDeg=opts.SigmaSurroundDeg;
calibration.candidateSpacingDeg=opts.CandidateSpacingDeg;
calibration.inputSpacingDeg=opts.InputSpacingDeg;
calibration.idqPreStepPedestalPC=opts.IDQPreStepPedestalPC;
calibration.nCandidates=numel(candidateDirections);
calibration.nInputDirections=numel(inputDirections);
calibration.noiseRank=noiseModel.rank;
calibration.targetPerformance=opts.TargetPerformance;
calibration.nMonteCarlo=opts.NMonteCarlo;
calibration.rectifyCandidates=opts.RectifyCandidates;
calibration.sessionTable=sessions;
calibration.mtCorrelationTable=correlationTable;
calibration.noiseDefinition=noiseModel.definition;
calibration.temporalDefinition=[ ...
  "IDQ has a three-direction 7% pedestal throughout prestep on every " + ...
  "trial. It continues during hasStepNoise trials and terminates at step " + ...
  "on no-noise psychometric trials. Calibration uses the decision epoch; " + ...
  "prestep values are retained for later temporal models."];

fig=makePlot(calibration,opts.Visible);
if opts.SaveOutputs
  dataFolder=analysisPath(context,'Common Code','Data','Simulation');
  plotFolder=analysisPath(context,'Common Code','Plots','Simulation');
  if ~isfolder(dataFolder),mkdir(dataFolder);end
  if ~isfolder(plotFolder),mkdir(plotFolder);end
  matPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseCalibration.mat');
  csvPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseCalibration.csv');
  corrPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseCorrelation.csv');
  pdfPath=fullfile(plotFolder,'IDRIDQ_MTUpstreamNoiseCalibration.pdf');
  calibration.outputPaths=string({matPath;csvPath;corrPath;pdfPath});
  info=whos('calibration');assertSimulationRunCapacity(context,info.bytes);
  save(matPath,'calibration','-v7.3');writetable(sessions,csvPath);
  writetable(correlationTable,corrPath);
  exportgraphics(fig,pdfPath,'ContentType','vector');
  assertSimulationRunCapacity(context,0);
end
if strcmpi(opts.Visible,'off'),close(fig);end
end

function p=simulatedPerformance(changeSignal,noChangeSignal, ...
  unitChange,unitNoChange,sigma,rectify)
A=changeSignal+sigma.*unitChange;B=noChangeSignal+sigma.*unitNoChange;
if rectify,A=max(A,0);B=max(B,0);end
p=mean(max(A,[],2)>max(B,[],2));
end

function [sigma,p,lo,hi]=solveNoiseScale(objective,threshold,kernel,target,nIter)
preferred=max(kernel);z=sqrt(2)*erfinv(2*target-1);
analytic=max(preferred*threshold/(sqrt(2)*z),eps);
lo=0;hi=analytic;
while objective(hi)>target&&hi<analytic*2^20,hi=hi*2;end
if objective(hi)>target
  error('calibrateIDRIDQMTUpstreamNoise:CouldNotBracket', ...
    'Could not bracket target performance.');
end
for j=1:nIter
  mid=(lo+hi)/2;
  if objective(mid)>target,lo=mid;else,hi=mid;end
end
sigma=(lo+hi)/2;p=objective(sigma);
end

function T=makeMTCorrelationTable(N,mtModel)
phi=double(mtModel.phiDeg(:));i0=find(abs(phi)<1e-9,1,'first');
separationDeg=(0:180)';correlation=nan(size(separationDeg));
for k=1:numel(separationDeg)
  d=mod(phi-separationDeg(k)+180,360)-180;
  j=find(abs(d)<1e-9,1,'first');
  correlation(k)=N.mtCorrelation(i0,j);
end
T=table(separationDeg,correlation);
end

function fig=makePlot(C,visible)
T=C.sessionTable;fig=figure('Color','w','Visible',visible, ...
  'Units','inches','Position',[1 1 12 4.3]);
tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('Upstream directional-noise calibration (%d MT inputs, %d readouts)', ...
  C.nInputDirections,C.nCandidates));
ax=nexttile(tl);hold(ax,'on');plotDatasets(ax,T,'thresholdPC','upstreamDirectionNoiseSD');
xlabel(ax,'Session 75% threshold (% coherence)');ylabel(ax,'Upstream directional-noise SD');
grid(ax,'on');box(ax,'off');legend(ax,'Location','best');
ax=nexttile(tl);hold(ax,'on');plotDatasets(ax,T,'baseToThresholdRatio','upstreamNoisePerThreshold');
xlabel(ax,'Decision base coherence / threshold');ylabel(ax,'Upstream noise SD / threshold');
grid(ax,'on');box(ax,'off');legend(ax,'Location','best');
ax=nexttile(tl);plot(ax,C.mtCorrelationTable.separationDeg, ...
  C.mtCorrelationTable.correlation,'k-','LineWidth',1.5);yline(ax,0,'k:');
xlabel(ax,'MT preferred-direction separation (deg)');ylabel(ax,'MT response-noise correlation');
grid(ax,'on');box(ax,'off');
end

function plotDatasets(ax,T,xName,yName)
x=T.(xName);y=T.(yName);
for dataset=["IDQ","IDR"]
  use=T.dataset==dataset;
  scatter(ax,x(use),y(use),20,'filled','DisplayName',dataset);
end
end
function requireCircleTiling(x,name)
if abs(round(360/x)*x-360)>1e-9
  error('calibrateIDRIDQMTUpstreamNoise:SpacingMustTileCircle', ...
    '%s must divide 360 exactly.',name);
end
end
function path=defaultSummaryPath(domain,file)
path=fullfile(domainFolder(mfilename('fullpath'),domain), ...
  'Data','AcrossSessionSummaries',file);
end
function tf=isPositiveScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0;
end
function tf=isNonnegativeScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0;
end
function tf=isOptionalPositiveScalar(x)
tf=isnumeric(x)&&isscalar(x)&&(isnan(x)||(isfinite(x)&&x>0));
end
