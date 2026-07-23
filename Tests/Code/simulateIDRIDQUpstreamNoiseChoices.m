function simulation=simulateIDRIDQUpstreamNoiseChoices( ...
  context,calibration,varargin)
% simulateIDRIDQUpstreamNoiseChoices  Replicated choices on actual trials.

context=validateAnalysisContext(context);
if context.Mode~="synthetic"
  error('IDReadout:SyntheticContextRequired', ...
    'Choice simulation requires a synthetic context.');
end
validateSyntheticManifest(context);
if ~isstruct(calibration)||~isfield(calibration,'context')|| ...
    string(calibration.context.RunID)~=context.RunID
  error('simulateIDRIDQUpstreamNoiseChoices:CalibrationContextMismatch', ...
    'Calibration must belong to the supplied synthetic run.');
end
required={'sourceRecords','sessionTable','sigmaMTDeg','sigmaReadoutDeg', ...
  'candidateSpacingDeg','inputSpacingDeg','idqPreStepPedestalPC', ...
  'rectifyCandidates'};
if ~all(isfield(calibration,required))
  error('simulateIDRIDQUpstreamNoiseChoices:InvalidCalibration', ...
    'Calibration lacks required upstream-noise fields.');
end

p=inputParser;
addParameter(p,'NReplicates',20, ...
  @(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&fix(x)==x);
addParameter(p,'ChoiceSeedOffset',104729, ...
  @(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&fix(x)==x);
addParameter(p,'RectifyCandidates',calibration.rectifyCandidates, ...
  @(x) islogical(x)&&isscalar(x));
addParameter(p,'ChunkSize',3000, ...
  @(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&fix(x)==x);
addParameter(p,'Visible','off',@(x) any(strcmpi(string(x),["on","off"])));
addParameter(p,'SaveOutputs',true,@(x) islogical(x)&&isscalar(x));
parse(p,varargin{:});opts=p.Results;

sourcePaths=string({calibration.sourceRecords.Path});sourcePaths=sourcePaths(:);
validateSimulationSourceFiles(context,sourcePaths);
if numel(sourcePaths)~=2
  error('simulateIDRIDQUpstreamNoiseChoices:ExpectedTwoSources', ...
    'Expected the IDQ and IDR summary sources in the calibration.');
end
conditions=makeIDRIDQMTTrialInputs( ...
  sourcePaths(1),sourcePaths(2),calibration.idqPreStepPedestalPC);
mtModel=makeMTReadoutForwardModel('sigmaMTDeg',calibration.sigmaMTDeg);
candidateDirections=(0:calibration.candidateSpacingDeg: ...
  (360-calibration.candidateSpacingDeg))';
inputDirections=(0:calibration.inputSpacingDeg: ...
  (360-calibration.inputSpacingDeg))';
if isfield(calibration,'negativeAsymptoteMagnitude')
  negativeAsymptoteMagnitude=calibration.negativeAsymptoteMagnitude;
else
  negativeAsymptoteMagnitude=0;
end
if isfield(calibration,'surroundMagnitude')
  surroundMagnitude=calibration.surroundMagnitude;
else
  surroundMagnitude=0;
end
if isfield(calibration,'sigmaSurroundDeg')
  sigmaSurroundDeg=calibration.sigmaSurroundDeg;
else
  sigmaSurroundDeg=nan;
end
bank=makeGaussianReadoutBank(candidateDirections, ...
  calibration.sigmaReadoutDeg,mtModel, ...
  'NegativeAsymptoteMagnitude',negativeAsymptoteMagnitude, ...
  'SurroundMagnitude',surroundMagnitude, ...
  'SigmaSurroundDeg',sigmaSurroundDeg);
noiseModel=factorUpstreamDirectionalNoise(inputDirections,bank,mtModel);

choiceSeed=mod(double(context.Seed)+double(opts.ChoiceSeedOffset),2^32);
stream=RandStream('mt19937ar','Seed',choiceSeed);
choiceResults=cell(numel(conditions),1);
for k=1:numel(conditions)
  sigma=mapSessionNoise(conditions(k),calibration.sessionTable);
  choiceResults{k}=simulateUpstreamChoiceCondition(conditions(k),bank, ...
    noiseModel,mtModel,sigma,opts.NReplicates,stream, ...
    'RectifyCandidates',opts.RectifyCandidates,'ChunkSize',opts.ChunkSize);
end

[conditionSummary,replicateSummary,sessionSummary]= ...
  summarizeChoices(choiceResults,opts.NReplicates);
noiseMetricTable=makeNoiseMetricTable(calibration,noiseModel,bank,mtModel);

simulation=struct();simulation.createdAt=datetime('now');
simulation.createdBy=mfilename;
simulation.dataOrigin="synthetic_choices_from_actual_noise_trials";
simulation.context=context;simulation.calibrationCreatedAt=calibration.createdAt;
simulation.choiceSeed=choiceSeed;simulation.choiceSeedOffset=opts.ChoiceSeedOffset;
simulation.nReplicates=opts.NReplicates;
simulation.rectifyCandidates=opts.RectifyCandidates;
simulation.poolingMode="hardMax";
simulation.sigmaMTDeg=calibration.sigmaMTDeg;
simulation.sigmaReadoutDeg=calibration.sigmaReadoutDeg;
simulation.candidateSpacingDeg=calibration.candidateSpacingDeg;
simulation.inputSpacingDeg=calibration.inputSpacingDeg;
simulation.negativeAsymptoteMagnitude=negativeAsymptoteMagnitude;
simulation.surroundMagnitude=surroundMagnitude;
simulation.sigmaSurroundDeg=sigmaSurroundDeg;
simulation.choiceResults=choiceResults;
simulation.conditionSummary=conditionSummary;
simulation.replicateSummary=replicateSummary;
simulation.sessionSummary=sessionSummary;
simulation.noiseMetricTable=noiseMetricTable;
simulation.noiseSamplingDefinition=[ ...
  "Candidate perturbations are sampled from the covariance obtained by " + ...
  "passing independent upstream directional noise through MT tuning and " + ...
  "the specified Gaussian-asymptote readouts. This is an exact " + ...
  "covariance-space shortcut, not a " + ...
  "separate readout-noise source."];

fig=makePlot(simulation,opts.Visible);
if opts.SaveOutputs
  dataFolder=analysisPath(context,'Common Code','Data','Simulation');
  plotFolder=analysisPath(context,'Common Code','Plots','Simulation');
  if ~isfolder(dataFolder),mkdir(dataFolder);end
  if ~isfolder(plotFolder),mkdir(plotFolder);end
  matPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseChoices.mat');
  conditionPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseChoiceConditions.csv');
  replicatePath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseChoiceReplicates.csv');
  sessionPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseChoiceSessions.csv');
  metricPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseMetrics.csv');
  pdfPath=fullfile(plotFolder,'IDRIDQ_MTUpstreamNoiseChoices.pdf');
  paths=string({matPath;conditionPath;replicatePath;sessionPath;metricPath;pdfPath});
  if any(isfile(paths))
    error('simulateIDRIDQUpstreamNoiseChoices:OutputExists', ...
      'Choice output already exists in this run; refusing to overwrite it.');
  end
  simulation.outputPaths=paths;
  info=whos('simulation');assertSimulationRunCapacity(context,info.bytes);
  save(matPath,'simulation','-v7.3');
  writetable(conditionSummary,conditionPath);writetable(replicateSummary,replicatePath);
  writetable(sessionSummary,sessionPath);writetable(noiseMetricTable,metricPath);
  exportgraphics(fig,pdfPath,'ContentType','vector');
  assertSimulationRunCapacity(context,0);
end
if strcmpi(opts.Visible,'off'),close(fig);end
end

function sigma=mapSessionNoise(C,T)
keys=string(T.dataset)+":"+string(T.sessionIndex);
trialKeys=string(C.dataset)+":"+string(C.sessionIndex);
[tf,row]=ismember(trialKeys,keys);
if ~all(tf)
  error('simulateIDRIDQUpstreamNoiseChoices:UncalibratedSession', ...
    '%d trials lack calibrated session noise.',sum(~tf));
end
sigma=double(T.upstreamDirectionNoiseSD(row));
end

function [Csummary,Rsummary,Ssummary]=summarizeChoices(results,nRep)
conditionRows=cell(numel(results),1);repTables=cell(numel(results),1);
sessionTables=cell(numel(results),1);
for k=1:numel(results)
  R=results{k};rep=double(R.replicateCorrectFraction);
  q=localQuantile(rep,[.05 .5 .95]);
  conditionRows{k}=struct('dataset',string(R.dataset), ...
    'probeOffsetDeg',R.offsetDeg,'nTrials',numel(R.trialIndex), ...
    'nSessions',numel(unique(R.sessionIndex)), ...
    'nReplicates',nRep,'experimentalCorrectFraction', ...
    mean(double(R.experimentalCorrect)),'syntheticMeanCorrectFraction',mean(rep), ...
    'syntheticReplicateSD',std(rep),'syntheticQ05',q(1), ...
    'syntheticMedian',q(2),'syntheticQ95',q(3));
  repTables{k}=table(repmat(string(R.dataset),nRep,1), ...
    repmat(R.offsetDeg,nRep,1),(1:nRep)',rep, ...
    'VariableNames',{'dataset','probeOffsetDeg','replicate','correctFraction'});
  sessions=unique(R.sessionIndex);rows=cell(numel(sessions),1);
  for j=1:numel(sessions)
    use=R.sessionIndex==sessions(j);
    rows{j}=struct('dataset',string(R.dataset),'probeOffsetDeg',R.offsetDeg, ...
      'sessionIndex',sessions(j),'nTrials',sum(use), ...
      'experimentalCorrectFraction',mean(double(R.experimentalCorrect(use))), ...
      'syntheticMeanCorrectFraction',mean(R.trialProbabilityCorrect(use)));
  end
  sessionTables{k}=struct2table(vertcat(rows{:}));
end
Csummary=struct2table(vertcat(conditionRows{:}));
Rsummary=vertcat(repTables{:});Ssummary=vertcat(sessionTables{:});
end

function T=makeNoiseMetricTable(C,N,bank,mtModel)
S=C.sessionTable;
iMT=find(abs(double(mtModel.phiDeg))<1e-9,1,'first');
iReadout=find(abs(double(bank.candidateDirectionsDeg))<1e-9,1,'first');
unitMTSD=sqrt(N.mtCovariance(iMT,iMT));
unitReadoutSD=sqrt(N.candidateCovariance(iReadout,iReadout));
template=mtPopulationTemplate(0,mtModel);
mtSensitivity=template(iMT);
readoutSensitivity=template*double(bank.weightsPhi(iReadout,:))';
T=S(:,{'dataset','sessionIndex','thresholdPC','upstreamDirectionNoiseSD'});
T.inducedPreferredMTNoiseSD=S.upstreamDirectionNoiseSD.*unitMTSD;
T.inducedPreferredReadoutNoiseSD=S.upstreamDirectionNoiseSD.*unitReadoutSD;
T.preferredMTResponsePerPC=repmat(mtSensitivity,height(S),1);
T.preferredReadoutResponsePerPC=repmat(readoutSensitivity,height(S),1);
T.mtNoiseInOnePCUnits=T.inducedPreferredMTNoiseSD./mtSensitivity;
T.readoutNoiseInOnePCUnits=T.inducedPreferredReadoutNoiseSD./readoutSensitivity;
end

function fig=makePlot(S,visible)
C=S.conditionSummary;R=S.replicateSummary;M=S.noiseMetricTable;
fig=figure('Color','w','Visible',visible,'Units','inches','Position',[1 1 12 4.3]);
tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('Upstream-noise hard-max choices (%d replicates)',S.nReplicates));
ax=nexttile(tl);hold(ax,'on');
datasets=unique(C.dataset,'stable');colors=lines(numel(datasets));
for k=1:numel(datasets)
  use=C.dataset==datasets(k);
  scatter(ax,C.experimentalCorrectFraction(use), ...
    C.syntheticMeanCorrectFraction(use),50,colors(k,:),'filled', ...
    'DisplayName',datasets(k));
end
plot(ax,[0 1],[0 1],'k:','HandleVisibility','off');
xlabel(ax,'Experimental correct fraction');ylabel(ax,'Synthetic correct fraction');
xlim(ax,[.4 1]);ylim(ax,[.4 1]);grid(ax,'on');box(ax,'off');
legend(ax,'Location','best');
ax=nexttile(tl);hold(ax,'on');
datasets=unique(R.dataset,'stable');colors=lines(numel(datasets));
for k=1:numel(datasets)
  use=R.dataset==datasets(k);
  scatter(ax,R.probeOffsetDeg(use),R.correctFraction(use),12,colors(k,:), ...
    'filled','MarkerFaceAlpha',.35,'DisplayName',datasets(k));
end
xlabel(ax,'Probe offset (deg)');ylabel(ax,'Replicate correct fraction');
grid(ax,'on');box(ax,'off');legend(ax,'Location','best');
ax=nexttile(tl);hold(ax,'on');
for dataset=["IDQ","IDR"]
  use=M.dataset==dataset;
  scatter(ax,M.thresholdPC(use),M.readoutNoiseInOnePCUnits(use),20,'filled', ...
    'DisplayName',dataset);
end
xlabel(ax,'Session threshold (% coherence)');
ylabel(ax,'Preferred-readout noise SD (equivalent % coherence)');
grid(ax,'on');box(ax,'off');legend(ax,'Location','best');
end

function q=localQuantile(x,p)
x=sort(double(x(isfinite(x))));
if isempty(x),q=nan(size(p));return,end
idx=1+(numel(x)-1)*p;lo=floor(idx);hi=ceil(idx);
q=x(lo).*(hi-idx)+x(hi).*(idx-lo);same=lo==hi;q(same)=x(lo(same));
end
