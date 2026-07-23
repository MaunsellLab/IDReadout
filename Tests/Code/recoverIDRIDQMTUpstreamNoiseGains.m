function recovery=recoverIDRIDQMTUpstreamNoiseGains(context,varargin)
% recoverIDRIDQMTUpstreamNoiseGains  Refit gains to replicated choices.
%
% Uses the experimental predictors, session thresholds, pooled psychometric
% shape, and physical coherence signs. Only trial correctness is replaced.
% A small checkpoint is committed after every completed replicate.

context=validateAnalysisContext(context);
if context.Mode~="synthetic"
  error('IDReadout:SyntheticContextRequired', ...
    'Gain recovery requires a synthetic context.');
end
validateSyntheticManifest(context);

p=inputParser;
addParameter(p,'SimulationPath',defaultSimulationPath(context), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'IDQSummaryPath',defaultSummaryPath('IDQ','IDQ_AcrossSideSummary.mat'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'IDRSummaryPath',defaultSummaryPath('IDR','IDR_SideGainAnalysis.mat'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'TargetPerformance',0.75, ...
  @(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>.5&&x<1);
addParameter(p,'GainBounds',[-5 5], ...
  @(x) isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(1)<x(2));
addParameter(p,'MaxNewReplicates',inf,@validMaxNewReplicates);
addParameter(p,'Visible','off',@(x) any(strcmpi(string(x),["on","off"])));
addParameter(p,'SaveOutputs',true,@(x) islogical(x)&&isscalar(x));
parse(p,varargin{:});opts=p.Results;

if exist('fitIDQNoiseGain','file')~=2
  error('recoverIDRIDQMTUpstreamNoiseGains:MissingIDQFitter', ...
    'fitIDQNoiseGain.m is required on the MATLAB path.');
end
if exist('analyzeIDRMatchedGains','file')~=2
  error('recoverIDRIDQMTUpstreamNoiseGains:MissingIDRFitter', ...
    'analyzeIDRMatchedGains.m is required on the MATLAB path.');
end
if ~isfile(opts.SimulationPath)
  error('recoverIDRIDQMTUpstreamNoiseGains:MissingSimulation', ...
    'Synthetic-choice MAT file was not found: %s',opts.SimulationPath);
end
sourcePaths=string({opts.IDQSummaryPath;opts.IDRSummaryPath});
sourceRecords=validateSimulationSourceFiles(context,sourcePaths);
sourcePaths=string({sourceRecords.Path})';
Q=load(sourcePaths(1),'acrossSummary');R=load(sourcePaths(2),'inventory');
S=load(opts.SimulationPath,'simulation');
if ~isfield(Q,'acrossSummary')||~isfield(R,'inventory')||~isfield(S,'simulation')
  error('recoverIDRIDQMTUpstreamNoiseGains:InvalidInputFiles', ...
    'Expected acrossSummary, inventory, and simulation variables.');
end
simulation=S.simulation;
if string(simulation.context.RunID)~=context.RunID
  error('recoverIDRIDQMTUpstreamNoiseGains:RunMismatch', ...
    'Synthetic choices belong to another run.');
end

[qBase,qCorrect]=makeIDQBase(Q.acrossSummary,simulation);
[rBase,rCorrect]=makeIDRBase(R.inventory,simulation);
nReplicates=simulation.nReplicates;
if size(qCorrect,2)~=nReplicates||size(rCorrect,2)~=nReplicates
  error('recoverIDRIDQMTUpstreamNoiseGains:ReplicateCountMismatch', ...
    'Choice matrices do not match simulation.nReplicates.');
end

psych=Q.acrossSummary.psychFit;
requiredPsych={'sessionFits','alignedWeibull'};
if ~all(isfield(psych,requiredPsych))
  error('recoverIDRIDQMTUpstreamNoiseGains:MissingIDQPsychometric', ...
    'IDQ psychFit lacks sessionFits or alignedWeibull.');
end
if ~isfield(R.inventory,'alignedPsychometric')|| ...
    ~isfield(R.inventory.alignedPsychometric,'primary')
  error('recoverIDRIDQMTUpstreamNoiseGains:MissingIDRPsychometric', ...
    'IDR inventory lacks alignedPsychometric.primary.');
end
idrPsych=R.inventory.alignedPsychometric.primary;

dataFolder=analysisPath(context,'Common Code','Data','Simulation');
plotFolder=analysisPath(context,'Common Code','Plots','Simulation');
if ~isfolder(dataFolder),mkdir(dataFolder);end
if ~isfolder(plotFolder),mkdir(plotFolder);end
checkpointPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseGainRecovery_checkpoint.mat');
finalMatPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseGainRecovery.mat');
idqCsvPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseGainRecovery_IDQ.csv');
idrCsvPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseGainRecovery_IDR.csv');
idqRepCsvPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseGainReplicates_IDQ.csv');
idrRepCsvPath=fullfile(dataFolder,'IDRIDQ_MTUpstreamNoiseGainReplicates_IDR.csv');
pdfPath=fullfile(plotFolder,'IDRIDQ_MTUpstreamNoiseGainRecovery.pdf');

if isfile(finalMatPath)
  F=load(finalMatPath,'recovery');recovery=F.recovery;
  verifyRecoveryIdentity(recovery,context,simulation);
  if opts.SaveOutputs
    writeMissingFinalOutputs(context,recovery,finalMatPath,idqCsvPath, ...
      idrCsvPath,idqRepCsvPath,idrRepCsvPath,pdfPath,opts.Visible);
  end
  fprintf('Loaded completed gain recovery: %s\n',finalMatPath);
  return
end

offsets=unique(double(rBase.trialTable.probeDirDeg),'sorted');
signature=makeSignature(context,simulation,sourcePaths,offsets);
state=initializeOrResume(checkpointPath,signature,nReplicates,offsets);
if isfield(state,'observed')
  observed=state.observed;
else
  fprintf('Fitting observed-choice reference gains...\n');
  observed=fitObserved(qBase,rBase,psych,idrPsych,opts);
  state.observed=observed;state.updatedAt=datetime('now');
  saveCheckpoint(context,checkpointPath,state);
end

newlyCompleted=0;
for iRep=1:nReplicates
  if state.completed(iRep),continue,end
  if newlyCompleted>=opts.MaxNewReplicates,break,end
  fprintf('Gain recovery replicate %d/%d...\n',iRep,nReplicates);
  qT=qBase;qT.correct=logical(qCorrect(:,iRep));
  qChange=fitIDQSide(qT,'Change',psych,opts.TargetPerformance);
  qNoChange=fitIDQSide(qT,'NoChange',psych,opts.TargetPerformance);
  state.idqGain(iRep,1,:)=reshape(qChange.gain,1,1,[]);
  state.idqGain(iRep,2,:)=reshape(qNoChange.gain,1,1,[]);
  state.idqNLL(iRep,:)=[getFitNLL(qChange) getFitNLL(qNoChange)];
  state.idqExitflag(iRep,:)=[getFitExitflag(qChange) getFitExitflag(qNoChange)];

  rD=rBase;rD.trialTable.correct=logical(rCorrect(:,iRep));
  rFit=analyzeIDRMatchedGains(rD,idrPsych.betaWeibull,idrPsych.lapse, ...
    opts.TargetPerformance,opts.GainBounds);
  [tf,row]=ismember(offsets,double(rFit.summaryTable.probeDirDeg));
  if ~all(tf),error('recoverIDRIDQMTUpstreamNoiseGains:OffsetMismatch', ...
      'IDR fitter did not return every simulated offset.');end
  T=rFit.summaryTable(row,:);
  state.idrGain(iRep,:,1)=T.gCP';state.idrGain(iRep,:,2)=T.gCQ';
  state.idrGain(iRep,:,3)=T.gNP';state.idrGain(iRep,:,4)=T.gNQ';
  state.idrNLL(iRep,:)=T.fullNLL';
  state.idrDeltaNLLProbe(iRep,:)=T.deltaNLLProbeTerms';
  state.idrExitflag(iRep,:)=T.fullExitflag';
  state.completed(iRep)=true;state.updatedAt=datetime('now');
  saveCheckpoint(context,checkpointPath,state);
  newlyCompleted=newlyCompleted+1;
end

recovery=assembleRecovery(context,simulation,sourcePaths,state,observed, ...
  qBase,rBase,opts);
if ~all(state.completed)
  fprintf('Checkpointed %d/%d gain-recovery replicates. Rerun to resume.\n', ...
    sum(state.completed),nReplicates);
  return
end

if opts.SaveOutputs
  finalPaths=string({finalMatPath;idqCsvPath;idrCsvPath; ...
    idqRepCsvPath;idrRepCsvPath;pdfPath});
  recovery.outputPaths=finalPaths;
  writeMissingFinalOutputs(context,recovery,finalMatPath,idqCsvPath, ...
    idrCsvPath,idqRepCsvPath,idrRepCsvPath,pdfPath,opts.Visible);
end
end

function [T,correct]=makeIDQBase(S,simulation)
idx=find(cellfun(@(x) string(x.dataset)=="IDQ",simulation.choiceResults));
if numel(idx)~=1,error('recoverIDRIDQMTUpstreamNoiseGains:IDQConditionCount', ...
    'Expected exactly one IDQ choice condition.');end
C=simulation.choiceResults{idx};source=S.trialTable;
required={'sessionIndex','trialIndex','correct','hasStepNoise','stepCoh','dirIndex'};
requireTableVariables(source,required,'IDQ trial table');
row=matchSyntheticTrialRows(source.sessionIndex,source.trialIndex,[], ...
  C.sessionIndex,C.trialIndex,[]);
T=source(row,:);correct=C.correctReplicates;
if ~all(logical(T.hasStepNoise))|| ...
    ~isequal(logical(T.correct),logical(C.experimentalCorrect))
  error('recoverIDRIDQMTUpstreamNoiseGains:IDQTrialMismatch', ...
    'Mapped IDQ trials do not match the simulated inventory.');
end
end

function [D,correct]=makeIDRBase(S,simulation)
if ~isfield(S,'sideGainData'),error( ...
    'recoverIDRIDQMTUpstreamNoiseGains:MissingIDRSideGainData', ...
    'IDR inventory lacks sideGainData.');end
source=S.sideGainData;T=source.trialTable;
required={'sessionIndex','trialIdx','probeDirDeg','correct'};
requireTableVariables(T,required,'IDR trial table');
idx=find(cellfun(@(x) string(x.dataset)=="IDR",simulation.choiceResults));
if isempty(idx),error('recoverIDRIDQMTUpstreamNoiseGains:NoIDRConditions', ...
    'No IDR choice conditions were found.');end
session=[];trial=[];offset=[];experimental=[];correctCell=cell(numel(idx),1);
for k=1:numel(idx)
  C=simulation.choiceResults{idx(k)};
  session=[session;double(C.sessionIndex(:))]; %#ok<AGROW>
  trial=[trial;double(C.trialIndex(:))]; %#ok<AGROW>
  offset=[offset;repmat(double(C.offsetDeg),numel(C.trialIndex),1)]; %#ok<AGROW>
  experimental=[experimental;logical(C.experimentalCorrect(:))]; %#ok<AGROW>
  correctCell{k}=logical(C.correctReplicates);
end
correct=vertcat(correctCell{:});
row=matchSyntheticTrialRows(T.sessionIndex,T.trialIdx,T.probeDirDeg, ...
  session,trial,offset);
if ~isequal(logical(T.correct(row)),experimental)
  error('recoverIDRIDQMTUpstreamNoiseGains:IDRTrialMismatch', ...
    'Mapped IDR trials do not match the simulated inventory.');
end
D=struct();D.trialTable=T(row,:);D.stepFrames=source.stepFrames;
fields={'changePrefNoiseByFrameTrial','changeProbeEffectiveNoiseByFrameTrial', ...
  'noChangePrefNoiseByFrameTrial','noChangeProbeEffectiveNoiseByFrameTrial'};
for k=1:numel(fields)
  if ~isfield(source,fields{k}),error( ...
      'recoverIDRIDQMTUpstreamNoiseGains:MissingIDRPredictor', ...
      'IDR sideGainData lacks %s.',fields{k});end
  predictorMatrix=source.(fields{k});
  D.(fields{k})=predictorMatrix(:,row);
end
end

function observed=fitObserved(qBase,rBase,psych,idrPsych,opts)
observed=struct();
observed.idqChange=fitIDQSide(qBase,'Change',psych,opts.TargetPerformance);
observed.idqNoChange=fitIDQSide(qBase,'NoChange',psych,opts.TargetPerformance);
observed.idr=analyzeIDRMatchedGains(rBase,idrPsych.betaWeibull, ...
  idrPsych.lapse,opts.TargetPerformance,opts.GainBounds);
end

function fit=fitIDQSide(T,token,psych,target)
for d=1:3
  source=sprintf('rectStep%sNoisePredDir%d',token,d);
  if ~ismember(source,T.Properties.VariableNames)
    error('recoverIDRIDQMTUpstreamNoiseGains:MissingIDQPredictor', ...
      'IDQ trial table lacks %s.',source);
  end
  T.(sprintf('noisePredDir%d',d))=T.(source);
end
F=fitIDQNoiseGain(T,psych.sessionFits,psych.alignedWeibull,target);
fit=F.driftNonDrift;
end

function state=initializeOrResume(path,signature,nRep,offsets)
if isfile(path)
  C=load(path,'state');state=C.state;
  if ~isfield(state,'signature')||~isequaln(state.signature,signature)
    error('recoverIDRIDQMTUpstreamNoiseGains:CheckpointMismatch', ...
      'Existing checkpoint does not describe these choices and sources.');
  end
  fprintf('Resuming gain recovery from %d/%d completed replicates.\n', ...
    sum(state.completed),nRep);
  return
end
state=struct();state.signature=signature;state.createdAt=datetime('now');
state.updatedAt=state.createdAt;state.completed=false(nRep,1);
state.idqGain=nan(nRep,2,2);state.idqNLL=nan(nRep,2);
state.idqExitflag=nan(nRep,2);
nOffset=numel(offsets);state.idrGain=nan(nRep,nOffset,4);
state.idrNLL=nan(nRep,nOffset);state.idrDeltaNLLProbe=nan(nRep,nOffset);
state.idrExitflag=nan(nRep,nOffset);
end

function signature=makeSignature(context,simulation,sources,offsets)
signature=struct('runID',context.RunID,'choiceSeed',simulation.choiceSeed, ...
  'nReplicates',simulation.nReplicates,'sourcePaths',sources(:), ...
  'offsets',offsets(:),'fitIDQNoiseGainPath',string(which('fitIDQNoiseGain')), ...
  'analyzeIDRMatchedGainsPath',string(which('analyzeIDRMatchedGains')));
end

function verifyRecoveryIdentity(R,context,simulation)
if ~isfield(R,'context')||string(R.context.RunID)~=context.RunID|| ...
    R.choiceSeed~=simulation.choiceSeed
  error('recoverIDRIDQMTUpstreamNoiseGains:FinalIdentityMismatch', ...
    'Completed recovery does not match the supplied run.');
end
end

function saveCheckpoint(context,path,state)
temporary=path+".temporary.mat";
info=whos('state');assertSimulationRunCapacity(context,info.bytes);
save(temporary,'state','-v7');
[ok,message]=movefile(temporary,path,'f');
if ~ok,error('recoverIDRIDQMTUpstreamNoiseGains:CheckpointCommitFailed', ...
    'Could not commit checkpoint: %s',message);end
assertSimulationRunCapacity(context,0);
end

function recovery=assembleRecovery(context,simulation,sources,state,observed,qBase,rBase,opts)
recovery=struct();recovery.createdAt=datetime('now');recovery.createdBy=mfilename;
recovery.dataOrigin="synthetic_choices_from_actual_noise_trials";
recovery.context=context;recovery.choiceSeed=simulation.choiceSeed;
recovery.nReplicates=simulation.nReplicates;recovery.completed=state.completed;
recovery.sourcePaths=sources;recovery.targetPerformance=opts.TargetPerformance;
recovery.gainBounds=opts.GainBounds(:)';
recovery.modelDefinition=["Experimental rectangular-step predictors, session " + ...
  "thresholds, and pooled psychometric shapes; only correctness is replaced " + ...
  "by hard-max upstream-noise synthetic choices."];
recovery.idqTrialCount=height(qBase);recovery.idrTrialCount=height(rBase.trialTable);
[recovery.idqReplicates,recovery.idqSummary]=makeIDQTables(state,observed);
[recovery.idrReplicates,recovery.idrSummary]=makeIDRTables(state,observed);
end

function [rep,summary]=makeIDQTables(S,O)
n=size(S.idqGain,1);replicate=repelem((1:n)',2);side=repmat(["change";"noChange"],n,1);
gDrift=reshape(S.idqGain(:,:,1)',[],1);gNonDrift=reshape(S.idqGain(:,:,2)',[],1);
negLogLikelihood=reshape(S.idqNLL',[],1);exitflag=reshape(S.idqExitflag',[],1);
rep=table(replicate,side,gDrift,gNonDrift,negLogLikelihood,exitflag);
observedFits={O.idqChange,O.idqNoChange};rows=cell(2,1);
sideNames=["change","noChange"];
for k=1:2
  use=side==sideNames(k);A=rep(use,:);F=observedFits{k};
  observedSE=getFitSE(F);
  d=distributionSummary(A.gDrift);nd=distributionSummary(A.gNonDrift);
  rows{k}=struct('side',sideNames(k), ...
    'observedGDrift',F.gain(1),'observedSEGDrift',observedSE(1), ...
    'syntheticMeanGDrift',d(1), ...
    'syntheticSDGDrift',d(2),'syntheticQ05GDrift',d(3), ...
    'syntheticMedianGDrift',d(4),'syntheticQ95GDrift',d(5), ...
    'observedGNonDrift',F.gain(2),'observedSEGNonDrift',observedSE(2), ...
    'syntheticMeanGNonDrift',nd(1), ...
    'syntheticSDGNonDrift',nd(2),'syntheticQ05GNonDrift',nd(3), ...
    'syntheticMedianGNonDrift',nd(4),'syntheticQ95GNonDrift',nd(5));
end
summary=struct2table(vertcat(rows{:}));
end

function [rep,summary]=makeIDRTables(S,O)
offsets=S.signature.offsets(:);nRep=size(S.idrGain,1);nOff=numel(offsets);
replicate=repelem((1:nRep)',nOff);probeDirDeg=repmat(offsets,nRep,1);
names={'gCP','gCQ','gNP','gNQ'};values=cell(1,4);
for j=1:4,values{j}=reshape(S.idrGain(:,:,j)',[],1);end
fullNLL=reshape(S.idrNLL',[],1);deltaNLLProbeTerms=reshape(S.idrDeltaNLLProbe',[],1);
exitflag=reshape(S.idrExitflag',[],1);
rep=table(replicate,probeDirDeg,values{1},values{2},values{3},values{4}, ...
  fullNLL,deltaNLLProbeTerms,exitflag,'VariableNames', ...
  {'replicate','probeDirDeg',names{:},'fullNLL','deltaNLLProbeTerms','exitflag'});
obs=O.idr.summaryTable;[tf,row]=ismember(offsets,double(obs.probeDirDeg));
if ~all(tf),error('recoverIDRIDQMTUpstreamNoiseGains:ObservedOffsetMismatch', ...
    'Observed IDR fits lack a simulated offset.');end
rows=cell(nOff,1);
for k=1:nOff
  use=rep.probeDirDeg==offsets(k);entry=struct('probeDirDeg',offsets(k));
  for j=1:4
    replicatedValues=rep.(names{j});
    d=distributionSummary(replicatedValues(use));
    observedValues=obs.(names{j});
    entry.(['observed' names{j}])=observedValues(row(k));
    seName=['se' names{j}(2:end)];seValues=obs.(seName);
    entry.(['observedSE' names{j}])=seValues(row(k));
    entry.(['syntheticMean' names{j}])=d(1);
    entry.(['syntheticSD' names{j}])=d(2);
    entry.(['syntheticQ05' names{j}])=d(3);
    entry.(['syntheticMedian' names{j}])=d(4);
    entry.(['syntheticQ95' names{j}])=d(5);
  end
  dNLL=distributionSummary(rep.deltaNLLProbeTerms(use));
  observedDelta=obs.deltaNLLProbeTerms(row(k));
  entry.observedDeltaNLLProbeTerms=observedDelta;
  entry.syntheticMeanDeltaNLLProbeTerms=dNLL(1);
  entry.syntheticSDDeltaNLLProbeTerms=dNLL(2);
  entry.syntheticQ05DeltaNLLProbeTerms=dNLL(3);
  entry.syntheticMedianDeltaNLLProbeTerms=dNLL(4);
  entry.syntheticQ95DeltaNLLProbeTerms=dNLL(5);
  rows{k}=entry;
end
summary=struct2table(vertcat(rows{:}));
end

function d=distributionSummary(x)
x=double(x(isfinite(x)));q=localQuantile(x,[.05 .5 .95]);
d=[mean(x) std(x) q(:)'];
end

function q=localQuantile(x,p)
x=sort(double(x(:)));if isempty(x),q=nan(size(p));return,end
idx=1+(numel(x)-1).*p;lo=floor(idx);hi=ceil(idx);
q=x(lo).*(hi-idx)+x(hi).*(idx-lo);same=lo==hi;q(same)=x(lo(same));
end

function figures=makePlots(R,visible)
figures=gobjects(2,1);C=lines(2);T=R.idqSummary;
figures(1)=figure('Color','w','Visible',visible,'Units','inches','Position',[1 1 10 4]);
tl=tiledlayout(figures(1),1,2,'TileSpacing','compact','Padding','compact');
title(tl,'IDQ replicated gain recovery');
for k=1:2
  ax=nexttile(tl);hold(ax,'on');x=1:2;
  mu=[T.syntheticMeanGDrift(k) T.syntheticMeanGNonDrift(k)];
  lo=[T.syntheticQ05GDrift(k) T.syntheticQ05GNonDrift(k)];
  hi=[T.syntheticQ95GDrift(k) T.syntheticQ95GNonDrift(k)];
  obs=[T.observedGDrift(k) T.observedGNonDrift(k)];
  errorbar(ax,x,mu,mu-lo,hi-mu,'o-','Color',C(k,:),'LineWidth',1.3, ...
    'DisplayName','Synthetic 5-95%');
  obsSE=[T.observedSEGDrift(k) T.observedSEGNonDrift(k)];
  errorbar(ax,x,obs,obsSE,'ks','MarkerFaceColor','k', ...
    'DisplayName','Observed +/- SE');yline(ax,0,'k:');
  set(ax,'XTick',x,'XTickLabel',{'Drift','Non-drift'});ylabel(ax,'Gain');
  title(ax,string(T.side(k)));grid(ax,'on');box(ax,'off');legend(ax,'Location','best');
end

figures(2)=figure('Color','w','Visible',visible,'Units','inches','Position',[1 1 10 7.5]);
tl=tiledlayout(figures(2),2,2,'TileSpacing','compact','Padding','compact');
title(tl,'IDR replicated four-gain recovery by offset');S=R.idrSummary;
names={'gCP','gCQ','gNP','gNQ'};titles={'Change preferred','Change probe', ...
  'No-change preferred','No-change probe'};
for j=1:4
  ax=nexttile(tl);hold(ax,'on');x=S.probeDirDeg;
  mu=S.(['syntheticMean' names{j}]);lo=S.(['syntheticQ05' names{j}]);
  hi=S.(['syntheticQ95' names{j}]);obs=S.(['observed' names{j}]);
  fill(ax,[x;flipud(x)],[lo;flipud(hi)],C(1,:), ...
    'FaceAlpha',.18,'EdgeColor','none','DisplayName','Synthetic 5-95%');
  plot(ax,x,mu,'o-','Color',C(1,:),'LineWidth',1.3,'DisplayName','Synthetic mean');
  obsSE=S.(['observedSE' names{j}]);
  errorbar(ax,x,obs,obsSE,'ks-','MarkerFaceColor','k', ...
    'DisplayName','Observed +/- SE');yline(ax,0,'k:');
  xlabel(ax,'Probe offset (deg)');ylabel(ax,'Gain');title(ax,titles{j});
  grid(ax,'on');box(ax,'off');if j==1,legend(ax,'Location','best');end
end
end

function writeMissingFinalOutputs(context,R,matPath,idqPath,idrPath, ...
  idqRepPath,idrRepPath,pdfPath,visible)
info=whos('R');assertSimulationRunCapacity(context,info.bytes);
if ~isfile(idqPath),writeTableAtomically(R.idqSummary,idqPath);end
if ~isfile(idrPath),writeTableAtomically(R.idrSummary,idrPath);end
if ~isfile(idqRepPath),writeTableAtomically(R.idqReplicates,idqRepPath);end
if ~isfile(idrRepPath),writeTableAtomically(R.idrReplicates,idrRepPath);end
if ~isfile(pdfPath)
  figures=makePlots(R,visible);
  temporaryPDF=pdfPath+".temporary.pdf";
  exportgraphics(figures(1),temporaryPDF,'ContentType','vector');
  exportgraphics(figures(2),temporaryPDF,'ContentType','vector','Append',true);
  if strcmpi(visible,'off'),close(figures);end
  commitTemporaryFile(temporaryPDF,pdfPath);
end
if ~isfile(matPath)
  recovery=R; %#ok<NASGU>
  temporaryMAT=matPath+".temporary.mat";
  save(temporaryMAT,'recovery','-v7.3');
  commitTemporaryFile(temporaryMAT,matPath);
end
if ~isfile(matPath)
  error('recoverIDRIDQMTUpstreamNoiseGains:FinalSaveFailed', ...
    'Final recovery MAT file was not created.');
end
assertSimulationRunCapacity(context,0);
end

function writeTableAtomically(T,path)
temporary=path+".temporary.csv";
writetable(T,temporary);
commitTemporaryFile(temporary,path);
end

function commitTemporaryFile(temporary,final)
[ok,message]=movefile(temporary,final,'f');
if ~ok,error('recoverIDRIDQMTUpstreamNoiseGains:FinalCommitFailed', ...
    'Could not commit final output: %s',message);end
end

function nll=getFitNLL(fit)
for name={'negLogLikelihood','NLL','nll'}
  if isfield(fit,name{1}),nll=fit.(name{1});return,end
end
error('recoverIDRIDQMTUpstreamNoiseGains:MissingFitNLL', ...
  'IDQ gain fit has no recognized negative-log-likelihood field.');
end

function exitflag=getFitExitflag(fit)
if isfield(fit,'exitflag'),exitflag=fit.exitflag;else,exitflag=nan;end
end

function se=getFitSE(fit)
if isfield(fit,'SE'),se=fit.SE(:);return,end
if isfield(fit,'CI95')&&size(fit.CI95,2)==2
  se=(fit.CI95(:,2)-fit.CI95(:,1))./(2*1.96);return
end
se=nan(numel(fit.gain),1);
end

function requireTableVariables(T,names,label)
missing=setdiff(names,T.Properties.VariableNames);
if ~isempty(missing),error('recoverIDRIDQMTUpstreamNoiseGains:MissingVariable', ...
    '%s is missing: %s',label,strjoin(missing,', '));end
end

function path=defaultSimulationPath(context)
path=analysisPath(context,'Common Code','Data','Simulation', ...
  'IDRIDQ_MTUpstreamNoiseChoices.mat');
end

function path=defaultSummaryPath(domain,file)
path=fullfile(domainFolder(mfilename('fullpath'),domain), ...
  'Data','AcrossSessionSummaries',file);
end

function tf=validMaxNewReplicates(x)
tf=isnumeric(x)&&isscalar(x)&& ...
  ((isfinite(x)&&x>=1&&fix(x)==x)||(isinf(x)&&x>0));
end
