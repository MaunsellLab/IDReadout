function sweep=runIDRIDQMTReadoutWidthSweep(varargin)
% runIDRIDQMTReadoutWidthSweep  Resumable Gaussian readout-width sweep.

% Each width gets an independent manifest-backed run. Completed calibration,
% choice, and recovery stages are loaded rather than recomputed. A fifth
% summary run declares the member result MAT files as hashed sources.

p=inputParser;
addParameter(p,'RunIDPrefix',"",@(x) ischar(x)||isstring(x));
addParameter(p,'ReadoutWidthsDeg',[1 5 10 20],@validWidths);
addParameter(p,'NReplicates',5,@positiveInteger);
addParameter(p,'Seed',1729,@validSeed);
addParameter(p,'SigmaMTDeg',37.5,@positiveScalar);
addParameter(p,'CandidateSpacingDeg',2,@positiveScalar);
addParameter(p,'InputSpacingDeg',1,@positiveScalar);
addParameter(p,'IDQPreStepPedestalPC',7,@nonnegativeScalar);
addParameter(p,'TargetPerformance',0.75, ...
  @(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>.5&&x<1);
addParameter(p,'NMonteCarlo',2000,@positiveInteger);
addParameter(p,'RectifyCandidates',true,@(x) islogical(x)&&isscalar(x));
addParameter(p,'MaxRunBytes',2e8,@positiveScalar);
addParameter(p,'IDQSummaryPath',defaultSummaryPath('IDQ','IDQ_AcrossSideSummary.mat'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'IDRSummaryPath',defaultSummaryPath('IDR','IDR_SideGainAnalysis.mat'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'Visible','off',@(x) any(strcmpi(string(x),["on","off"])));
parse(p,varargin{:});opts=p.Results;

prefix=string(opts.RunIDPrefix);
if strlength(prefix)==0
  error('runIDRIDQMTReadoutWidthSweep:RunIDPrefixRequired', ...
    'Supply a unique RunIDPrefix.');
end
widths=unique(double(opts.ReadoutWidthsDeg(:)'),'stable');
sourcePaths=string({opts.IDQSummaryPath;opts.IDRSummaryPath});
if ~all(isfile(sourcePaths))
  error('runIDRIDQMTReadoutWidthSweep:MissingExperimentalSource', ...
    'One or both experimental summary files do not exist.');
end

nWidth=numel(widths);members=repmat(emptyMember(),nWidth,1);
for k=1:nWidth
  sigmaR=widths(k);runID=formatReadoutWidthRunID(prefix,sigmaR);
  fprintf('\n=== Readout width %.6g deg (%d/%d): %s ===\n', ...
    sigmaR,k,nWidth,runID);
  context=makeAnalysisContext("synthetic",'RunID',runID,'Seed',opts.Seed, ...
    'MaxRunBytes',opts.MaxRunBytes,'SaveLatents',false);
  modelParameters=memberModelParameters(opts,sigmaR);
  ensureMemberRun(context,modelParameters,sourcePaths);

  calibrationPath=analysisPath(context,'Common Code','Data','Simulation', ...
    'IDRIDQ_MTUpstreamNoiseCalibration.mat');
  if isfile(calibrationPath)
    C=load(calibrationPath,'calibration');calibration=C.calibration;
    verifyStageContext(calibration,context,'calibration');
    if calibration.sigmaReadoutDeg~=sigmaR
      error('runIDRIDQMTReadoutWidthSweep:CalibrationWidthMismatch', ...
        'Stored calibration width does not match run %s.',runID);
    end
    fprintf('Loaded completed calibration.\n');
  else
    calibration=calibrateIDRIDQMTUpstreamNoise(context, ...
      'IDQSummaryPath',sourcePaths(1),'IDRSummaryPath',sourcePaths(2), ...
      'SigmaMTDeg',opts.SigmaMTDeg,'SigmaReadoutDeg',sigmaR, ...
      'CandidateSpacingDeg',opts.CandidateSpacingDeg, ...
      'InputSpacingDeg',opts.InputSpacingDeg, ...
      'IDQPreStepPedestalPC',opts.IDQPreStepPedestalPC, ...
      'TargetPerformance',opts.TargetPerformance, ...
      'NMonteCarlo',opts.NMonteCarlo, ...
      'RectifyCandidates',opts.RectifyCandidates,'Visible','off');
  end

  choicePath=analysisPath(context,'Common Code','Data','Simulation', ...
    'IDRIDQ_MTUpstreamNoiseChoices.mat');
  if isfile(choicePath)
    C=load(choicePath,'simulation');simulation=C.simulation;
    verifyStageContext(simulation,context,'simulation');
    if simulation.nReplicates~=opts.NReplicates
      error('runIDRIDQMTReadoutWidthSweep:ChoiceReplicateMismatch', ...
        'Existing choices for %s contain %d rather than %d replicates.', ...
        runID,simulation.nReplicates,opts.NReplicates);
    end
    fprintf('Loaded completed replicated choices.\n');
  else
    simulation=simulateIDRIDQUpstreamNoiseChoices(context,calibration, ...
      'NReplicates',opts.NReplicates, ...
      'RectifyCandidates',opts.RectifyCandidates,'Visible','off');
  end

  recovery=recoverIDRIDQMTUpstreamNoiseGains(context, ...
    'IDQSummaryPath',sourcePaths(1),'IDRSummaryPath',sourcePaths(2), ...
    'TargetPerformance',opts.TargetPerformance,'Visible','off');
  if ~all(recovery.completed)
    error('runIDRIDQMTReadoutWidthSweep:IncompleteRecovery', ...
      'Gain recovery for %s returned before all replicates completed.',runID);
  end
  recoveryPath=analysisPath(context,'Common Code','Data','Simulation', ...
    'IDRIDQ_MTUpstreamNoiseGainRecovery.mat');
  if ~isfile(recoveryPath)
    error('runIDRIDQMTReadoutWidthSweep:MissingRecoveryOutput', ...
      'Completed recovery MAT was not saved for %s.',runID);
  end
  members(k)=struct('sigmaReadoutDeg',sigmaR,'runID',runID, ...
    'context',context,'calibrationPath',string(calibrationPath), ...
    'choicePath',string(choicePath),'recoveryPath',string(recoveryPath));
end

sweep=makeOrLoadSummaryRun(prefix,members,opts);
end

function ensureMemberRun(context,parameters,sources)
manifestPath=analysisPath(context,'Tests','manifest.mat');
if isfile(manifestPath)
  manifest=validateSyntheticManifest(context);
  if string(manifest.ModelName)~="IDRIDQMTReadoutWidthSweepMember"|| ...
      ~isequaln(manifest.ModelParameters,parameters)
    error('runIDRIDQMTReadoutWidthSweep:MemberManifestMismatch', ...
      'Existing run %s has different model parameters.',context.RunID);
  end
  validateSimulationSourceFiles(context,sources);
  fprintf('Resuming existing member run.\n');
else
  initializeSyntheticRun(context, ...
    'ModelName',"IDRIDQMTReadoutWidthSweepMember", ...
    'ModelParameters',parameters,'SourceFiles',sources);
end
end

function P=memberModelParameters(O,sigmaR)
P=struct('sigmaMTDeg',double(O.SigmaMTDeg), ...
  'sigmaReadoutDeg',double(sigmaR), ...
  'candidateSpacingDeg',double(O.CandidateSpacingDeg), ...
  'inputSpacingDeg',double(O.InputSpacingDeg), ...
  'IDQPreStepPedestalPC',double(O.IDQPreStepPedestalPC), ...
  'targetPerformance',double(O.TargetPerformance), ...
  'nMonteCarlo',double(O.NMonteCarlo), ...
  'nChoiceReplicates',double(O.NReplicates), ...
  'rectifyCandidates',logical(O.RectifyCandidates), ...
  'poolingMode',"hardMax");
end

function sweep=makeOrLoadSummaryRun(prefix,members,opts)
summaryRunID=prefix+"_summary";
if strlength(summaryRunID)>80
  error('runIDRIDQMTReadoutWidthSweep:SummaryRunIDTooLong', ...
    'Summary RunID exceeds the 80-character context limit.');
end
context=makeAnalysisContext("synthetic",'RunID',summaryRunID, ...
  'Seed',opts.Seed,'MaxRunBytes',1e8,'SaveLatents',false);
sourcePaths=strings(2*numel(members),1);
for k=1:numel(members)
  sourcePaths(2*k-1)=members(k).choicePath;
  sourcePaths(2*k)=members(k).recoveryPath;
end
parameters=struct('readoutWidthsDeg',[members.sigmaReadoutDeg], ...
  'memberRunIDs',string({members.runID}), ...
  'nChoiceReplicates',double(opts.NReplicates), ...
  'summaryDefinition',"Cross-width descriptive comparison");
manifestPath=analysisPath(context,'Tests','manifest.mat');
if isfile(manifestPath)
  manifest=validateSyntheticManifest(context);
  if string(manifest.ModelName)~="IDRIDQMTReadoutWidthSweepSummary"|| ...
      ~isequaln(manifest.ModelParameters,parameters)
    error('runIDRIDQMTReadoutWidthSweep:SummaryManifestMismatch', ...
      'Existing summary run has different member parameters.');
  end
  validateSimulationSourceFiles(context,sourcePaths);
else
  initializeSyntheticRun(context, ...
    'ModelName',"IDRIDQMTReadoutWidthSweepSummary", ...
    'ModelParameters',parameters,'SourceFiles',sourcePaths);
end
finalPath=analysisPath(context,'Common Code','Data','Simulation', ...
  'IDRIDQ_MTReadoutWidthSweep.mat');
if isfile(finalPath)
  F=load(finalPath,'sweep');sweep=F.sweep;
  sweep=saveSweepOutputs(sweep,opts.Visible);
  fprintf('Loaded completed readout-width summary: %s\n',finalPath);
  return
end
sweep=buildSweepSummary(context,members,sourcePaths,opts);
sweep=saveSweepOutputs(sweep,opts.Visible);
end

function sweep=buildSweepSummary(context,members,sourcePaths,opts)
n=numel(members);idq=cell(n,1);idr=cell(n,1);performance=cell(n,1);metrics=cell(n,1);
for k=1:n
  A=load(members(k).recoveryPath,'recovery');B=load(members(k).choicePath,'simulation');
  width=members(k).sigmaReadoutDeg;
  idq{k}=addWidth(A.recovery.idqSummary,width);
  idr{k}=addWidth(A.recovery.idrSummary,width);
  performance{k}=addWidth(B.simulation.conditionSummary,width);
  metrics{k}=addWidth(B.simulation.noiseMetricTable,width);
end
sweep=struct();sweep.createdAt=datetime('now');sweep.createdBy=mfilename;
sweep.dataOrigin="synthetic_summary_from_manifested_member_runs";
sweep.context=context;sweep.memberRuns=members;sweep.sourcePaths=sourcePaths;
sweep.readoutWidthsDeg=[members.sigmaReadoutDeg];
sweep.nReplicatesPerWidth=opts.NReplicates;
sweep.idqGainTable=vertcat(idq{:});sweep.idrGainTable=vertcat(idr{:});
sweep.performanceTable=vertcat(performance{:});sweep.noiseMetricTable=vertcat(metrics{:});
sweep.readoutPropertyTable=makeReadoutPropertyTable( ...
  sweep.readoutWidthsDeg,opts.SigmaMTDeg);
sweep.definition=["Every width is independently threshold-calibrated; actual " + ...
  "trial predictors are fixed; only internal choices vary across replicates."];
end

function T=addWidth(T,width)
T.sigmaReadoutDeg=repmat(width,height(T),1);
T=movevars(T,'sigmaReadoutDeg','Before',1);
end

function sweep=saveSweepOutputs(sweep,visible)
context=sweep.context;dataFolder=analysisPath(context,'Common Code','Data','Simulation');
plotFolder=analysisPath(context,'Common Code','Plots','Simulation');
if ~isfolder(dataFolder),mkdir(dataFolder);end
if ~isfolder(plotFolder),mkdir(plotFolder);end
matPath=fullfile(dataFolder,'IDRIDQ_MTReadoutWidthSweep.mat');
idqPath=fullfile(dataFolder,'IDRIDQ_MTReadoutWidthSweep_IDQGains.csv');
idrPath=fullfile(dataFolder,'IDRIDQ_MTReadoutWidthSweep_IDRGains.csv');
performancePath=fullfile(dataFolder,'IDRIDQ_MTReadoutWidthSweep_Performance.csv');
metricPath=fullfile(dataFolder,'IDRIDQ_MTReadoutWidthSweep_NoiseMetrics.csv');
propertyPath=fullfile(dataFolder,'IDRIDQ_MTReadoutWidthSweep_ReadoutProperties.csv');
pdfPath=fullfile(plotFolder,'IDRIDQ_MTReadoutWidthSweep.pdf');
paths=string({matPath;idqPath;idrPath;performancePath;metricPath;propertyPath;pdfPath});
sweep.outputPaths=paths;info=whos('sweep');
assertSimulationRunCapacity(context,info.bytes);
if ~isfile(idqPath),writeTableAtomic(sweep.idqGainTable,idqPath);end
if ~isfile(idrPath),writeTableAtomic(sweep.idrGainTable,idrPath);end
if ~isfile(performancePath),writeTableAtomic(sweep.performanceTable,performancePath);end
if ~isfile(metricPath),writeTableAtomic(sweep.noiseMetricTable,metricPath);end
if ~isfile(propertyPath),writeTableAtomic(sweep.readoutPropertyTable,propertyPath);end
if ~isfile(pdfPath)
  figures=makeSweepPlots(sweep,visible);temporaryPDF=pdfPath+".temporary.pdf";
  exportgraphics(figures(1),temporaryPDF,'ContentType','vector');
  exportgraphics(figures(2),temporaryPDF,'ContentType','vector','Append',true);
  exportgraphics(figures(3),temporaryPDF,'ContentType','vector','Append',true);
  if strcmpi(visible,'off'),close(figures);end
  commitTemporary(temporaryPDF,pdfPath);
end
if ~isfile(matPath)
  temporaryMAT=matPath+".temporary.mat";save(temporaryMAT,'sweep','-v7.3');
  commitTemporary(temporaryMAT,matPath);
end
assertSimulationRunCapacity(context,0);
end

function figures=makeSweepPlots(S,visible)
colors=lines(numel(S.readoutWidthsDeg));figures=gobjects(3,1);
T=S.idqGainTable;
figures(1)=figure('Color','w','Visible',visible,'Units','inches','Position',[1 1 10 7]);
tl=tiledlayout(figures(1),2,2,'TileSpacing','compact','Padding','compact');
title(tl,'IDQ gains across Gaussian readout width');
spec={"change","GDrift","Change drift";"change","GNonDrift","Change non-drift"; ...
  "noChange","GDrift","No-change drift";"noChange","GNonDrift","No-change non-drift"};
for j=1:4
  ax=nexttile(tl);hold(ax,'on');use=T.side==spec{j,1};A=sortrows(T(use,:),'sigmaReadoutDeg');
  token=string(spec{j,2});x=A.sigmaReadoutDeg;
  meanName=char("syntheticMean"+token);
  q05Name=char("syntheticQ05"+token);
  q95Name=char("syntheticQ95"+token);
  observedName=char("observed"+token);
  mu=A.(meanName);lo=A.(q05Name);hi=A.(q95Name);
  errorbar(ax,x,mu,mu-lo,hi-mu,'o-','LineWidth',1.3, ...
    'DisplayName','Synthetic 5-95%');
  observedValues=A.(observedName);
  yline(ax,observedValues(1),'k--','LineWidth',1.2, ...
    'DisplayName','Observed');yline(ax,0,'k:','HandleVisibility','off');
  xlabel(ax,'Readout sigma (deg)');ylabel(ax,'Gain');title(ax,spec{j,3});
  grid(ax,'on');box(ax,'off');if j==1,legend(ax,'Location','best');end
end

T=S.idrGainTable;names={'gCP','gCQ','gNP','gNQ'};
titles={'Change preferred','Change probe','No-change preferred','No-change probe'};
figures(2)=figure('Color','w','Visible',visible,'Units','inches','Position',[1 1 10 7]);
tl=tiledlayout(figures(2),2,2,'TileSpacing','compact','Padding','compact');
title(tl,'IDR gains across Gaussian readout width');
for j=1:4
  ax=nexttile(tl);hold(ax,'on');name=string(names{j});
  syntheticName=char("syntheticMean"+name);
  observedName=char("observed"+name);
  observedSEName=char("observedSE"+name);
  for k=1:numel(S.readoutWidthsDeg)
    use=T.sigmaReadoutDeg==S.readoutWidthsDeg(k);A=sortrows(T(use,:),'probeDirDeg');
    plot(ax,A.probeDirDeg,A.(syntheticName),'o-', ...
      'Color',colors(k,:),'LineWidth',1.2, ...
      'DisplayName',sprintf('sigma_R = %g deg',S.readoutWidthsDeg(k)));
  end
  A=sortrows(T(T.sigmaReadoutDeg==S.readoutWidthsDeg(1),:),'probeDirDeg');
  errorbar(ax,A.probeDirDeg,A.(observedName),A.(observedSEName), ...
    'ks-','MarkerFaceColor','k','DisplayName','Observed');
  yline(ax,0,'k:','HandleVisibility','off');xlabel(ax,'Probe offset (deg)');
  ylabel(ax,'Gain');title(ax,titles{j});grid(ax,'on');box(ax,'off');
  if j==1,legend(ax,'Location','best');end
end

T=S.performanceTable;figures(3)=figure('Color','w','Visible',visible, ...
  'Units','inches','Position',[1 1 12 4]);
tl=tiledlayout(figures(3),1,3,'TileSpacing','compact','Padding','compact');
title(tl,'Performance across Gaussian readout width');
ax=nexttile(tl);hold(ax,'on');A=sortrows(T(T.dataset=="IDQ",:),'sigmaReadoutDeg');
plot(ax,A.sigmaReadoutDeg,A.syntheticMeanCorrectFraction,'o-','LineWidth',1.3, ...
  'DisplayName','Synthetic');yline(ax,A.experimentalCorrectFraction(1),'k--', ...
  'DisplayName','Observed');xlabel(ax,'Readout sigma (deg)');ylabel(ax,'P(correct)');
title(ax,'IDQ');grid(ax,'on');box(ax,'off');legend(ax,'Location','best');
ax=nexttile(tl);hold(ax,'on');
for k=1:numel(S.readoutWidthsDeg)
  use=T.dataset=="IDR"&T.sigmaReadoutDeg==S.readoutWidthsDeg(k);
  A=sortrows(T(use,:),'probeOffsetDeg');plot(ax,A.probeOffsetDeg, ...
    A.syntheticMeanCorrectFraction,'o-','Color',colors(k,:),'LineWidth',1.2, ...
    'DisplayName',sprintf('sigma_R = %g deg',S.readoutWidthsDeg(k)));
end
A=sortrows(T(T.dataset=="IDR"&T.sigmaReadoutDeg==S.readoutWidthsDeg(1),:), ...
  'probeOffsetDeg');plot(ax,A.probeOffsetDeg,A.experimentalCorrectFraction, ...
  'ks-','MarkerFaceColor','k','DisplayName','Observed');
xlabel(ax,'Probe offset (deg)');ylabel(ax,'P(correct)');title(ax,'IDR');
grid(ax,'on');box(ax,'off');legend(ax,'Location','best');
ax=nexttile(tl);A=S.readoutPropertyTable;
yyaxis(ax,'left');plot(ax,A.sigmaReadoutDeg,A.approxEffectiveSigmaDeg, ...
  'o-','LineWidth',1.3);ylabel(ax,'Approx. effective sigma (deg)');
yyaxis(ax,'right');plot(ax,A.sigmaReadoutDeg,A.preferredNullRatio, ...
  's-','LineWidth',1.3);ylabel(ax,'Preferred:null sensitivity');
xlabel(ax,'Readout sigma (deg)');title(ax,'Linear readout properties');
grid(ax,'on');box(ax,'off');
end

function T=makeReadoutPropertyTable(widths,sigmaMT)
mt=makeMTReadoutForwardModel('sigmaMTDeg',sigmaMT);n=numel(widths);
sigmaReadoutDeg=double(widths(:));sigmaMTDeg=repmat(double(sigmaMT),n,1);
preferredSensitivity=nan(n,1);nullSensitivity=nan(n,1);
for k=1:n
  bank=makeGaussianReadoutBank(0,sigmaReadoutDeg(k),mt);
  weights=double(bank.weightsPhi(1,:))';
  preferredSensitivity(k)=mtPopulationTemplate(0,mt)*weights;
  nullSensitivity(k)=mtPopulationTemplate(180,mt)*weights;
end
preferredNullRatio=preferredSensitivity./abs(nullSensitivity);
nullPreferredMagnitude=abs(nullSensitivity)./preferredSensitivity;
approxEffectiveSigmaDeg=sqrt(sigmaMTDeg.^2+sigmaReadoutDeg.^2);
T=table(sigmaReadoutDeg,sigmaMTDeg,approxEffectiveSigmaDeg, ...
  preferredSensitivity,nullSensitivity,preferredNullRatio,nullPreferredMagnitude);
end

function writeTableAtomic(T,path)
temporary=path+".temporary.csv";writetable(T,temporary);
commitTemporary(temporary,path);
end

function commitTemporary(temporary,final)
[ok,message]=movefile(temporary,final,'f');
if ~ok,error('runIDRIDQMTReadoutWidthSweep:OutputCommitFailed', ...
    'Could not commit summary output: %s',message);end
end

function verifyStageContext(S,context,label)
if ~isstruct(S)||~isfield(S,'context')||string(S.context.RunID)~=context.RunID
  error('runIDRIDQMTReadoutWidthSweep:StageContextMismatch', ...
    'Stored %s does not belong to run %s.',label,context.RunID);
end
end

function M=emptyMember()
M=struct('sigmaReadoutDeg',nan,'runID',"",'context',struct(), ...
  'calibrationPath',"",'choicePath',"",'recoveryPath',"");
end

function path=defaultSummaryPath(domain,file)
path=fullfile(domainFolder(mfilename('fullpath'),domain), ...
  'Data','AcrossSessionSummaries',file);
end

function tf=validWidths(x)
tf=isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>0);
end
function tf=positiveInteger(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&fix(x)==x;
end
function tf=validSeed(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=2^32-1&&fix(x)==x;
end
function tf=positiveScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0;
end
function tf=nonnegativeScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0;
end
