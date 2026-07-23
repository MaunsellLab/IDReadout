function summary = analyzeIDRIDQMTTrialCandidates(context,varargin)
% analyzeIDRIDQMTTrialCandidates  Deterministic actual-trial MT candidates.
%
% This diagnostic adds no internal variability and does not generate
% choices. It evaluates signed and candidate-rectified Gaussian readouts,
% with a hard maximum within each patch, over actual rectangular-step IDQ
% and IDR coherence distributions.

context = validateAnalysisContext(context);
if context.Mode ~= "synthetic"
  error('IDReadout:SyntheticContextRequired', ...
    'MT trial-candidate analysis requires a synthetic context.');
end
validateSyntheticManifest(context);

p=inputParser;
addParameter(p,'IDQSummaryPath',defaultSummaryPath('IDQ','IDQ_AcrossSideSummary.mat'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'IDRSummaryPath',defaultSummaryPath('IDR','IDR_SideGainAnalysis.mat'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'SigmaMTDeg',37.5,@isPositiveScalar);
addParameter(p,'IDQPedestalPC',nan, ...
  @(x) isnumeric(x)&&isvector(x)&&~isempty(x));
addParameter(p,'SigmaReadoutDeg',[15 20 25], ...
  @(x) isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0));
addParameter(p,'LatentSigmaReadoutDeg',20,@isPositiveScalar);
addParameter(p,'Visible','off',@(x) any(strcmpi(string(x),["on","off"])));
addParameter(p,'SaveOutputs',true,@(x) islogical(x)&&isscalar(x));
parse(p,varargin{:}); opts=p.Results;
if any(~isfinite(opts.IDQPedestalPC)) || any(opts.IDQPedestalPC<0)
  error('analyzeIDRIDQMTTrialCandidates:IDQPedestalRequired', ...
    ['IDQPedestalPC is required and must give the positive coherence ' ...
     'pedestal amplitude in percent coherence, either as one common value ' ...
     'or one value indexed by IDQ sessionIndex.']);
end

sourcePaths = string({opts.IDQSummaryPath;opts.IDRSummaryPath});
sourceRecords = validateSimulationSourceFiles(context,sourcePaths);
conditions = makeIDRIDQMTTrialInputs( ...
  sourcePaths(1),sourcePaths(2),opts.IDQPedestalPC);
mtModel = makeMTReadoutForwardModel('sigmaMTDeg',opts.SigmaMTDeg);
sigmaR = unique(double(opts.SigmaReadoutDeg(:)),'stable');
if context.SaveLatents && ~any(abs(sigmaR-opts.LatentSigmaReadoutDeg)<1e-9)
  error('analyzeIDRIDQMTTrialCandidates:LatentWidthNotSwept', ...
    'LatentSigmaReadoutDeg must be included in SigmaReadoutDeg.');
end

rows = cell(0,1); winnerRows=cell(0,1); signRows=cell(0,1);
latentResults = cell(0,1);
for iSigma=1:numel(sigmaR)
  for iRect=1:2
    rectify = iRect==2;
    for iCondition=1:numel(conditions)
      R = evaluateMTTrialCandidateCondition(conditions(iCondition), ...
        mtModel,sigmaR(iSigma),rectify);
      rows{end+1,1}=summarizeResult(R); %#ok<AGROW>
      winnerRows{end+1,1}=summarizeWinners(R); %#ok<AGROW>
      signRows{end+1,1}=summarizeNoiseSigns(R,conditions(iCondition)); %#ok<AGROW>
      if context.SaveLatents && abs(sigmaR(iSigma)-opts.LatentSigmaReadoutDeg)<1e-9
        latentResults{end+1,1}=R; %#ok<AGROW>
      end
    end
  end
end
diagnosticTable = struct2table(vertcat(rows{:}));
winnerRoleTable = vertcat(winnerRows{:});
signedNoiseAssociationTable = vertcat(signRows{:});

summary=struct();
summary.createdAt=datetime('now');
summary.createdBy=mfilename;
summary.dataOrigin="synthetic_from_manifested_experimental_predictors";
summary.context=context;
summary.sourceRecords=sourceRecords;
summary.sigmaMTDeg=opts.SigmaMTDeg;
summary.sigmaReadoutDeg=sigmaR';
summary.idqPedestalPC=double(opts.IDQPedestalPC(:)');
summary.latentSigmaReadoutDeg=opts.LatentSigmaReadoutDeg;
summary.diagnosticTable=diagnosticTable;
summary.winnerRoleTable=winnerRoleTable;
summary.signedNoiseAssociationTable=signedNoiseAssociationTable;
summary.modelDefinition=[ ...
  "Actual step-rectangular physical coherence including the positive IDQ " + ...
  "three-direction noise pedestals; signed mean-subtracted MT " + ...
  "population; equal unit-peak Gaussian direction readouts; physical plus " + ...
  "opponent candidates with coincident mechanisms collapsed; signed or " + ...
  "candidate-rectified activation; hard within-patch maximum; no internal " + ...
  "variability and no generated choices"];

fig=makePlot(summary,opts.LatentSigmaReadoutDeg,opts.Visible);
if opts.SaveOutputs
  dataFolder=analysisPath(context,'Common Code','Data','Simulation');
  plotFolder=analysisPath(context,'Common Code','Plots','Simulation');
  if ~isfolder(dataFolder),mkdir(dataFolder);end
  if ~isfolder(plotFolder),mkdir(plotFolder);end
  matPath=fullfile(dataFolder,'IDRIDQ_MTTrialCandidateSummary.mat');
  csvPath=fullfile(dataFolder,'IDRIDQ_MTTrialCandidateSummary.csv');
  winnerCsvPath=fullfile(dataFolder,'IDRIDQ_MTTrialCandidateWinners.csv');
  signCsvPath=fullfile(dataFolder,'IDRIDQ_MTTrialCandidateSignedNoise.csv');
  pdfPath=fullfile(plotFolder,'IDRIDQ_MTTrialCandidateSummary.pdf');
  paths=string({matPath;csvPath;winnerCsvPath;signCsvPath;pdfPath});
  if context.SaveLatents
    latentPath=fullfile(dataFolder,'IDRIDQ_MTTrialCandidateLatents.mat');
    paths(end+1,1)=latentPath;
  end
  summary.outputPaths=paths;
  bytes=whos('summary','latentResults');
  assertSimulationRunCapacity(context,sum([bytes.bytes]));
  save(matPath,'summary','-v7.3');
  writetable(diagnosticTable,csvPath);
  writetable(winnerRoleTable,winnerCsvPath);
  writetable(signedNoiseAssociationTable,signCsvPath);
  exportgraphics(fig,pdfPath,'ContentType','vector');
  if context.SaveLatents
    latentDefinition=summary.modelDefinition; %#ok<NASGU>
    save(latentPath,'latentResults','latentDefinition','-v7.3');
  end
  assertSimulationRunCapacity(context,0);
end
if strcmpi(opts.Visible,'off'),close(fig);end
end

function row=summarizeResult(R)
n=numel(R.patchMargin);
primary=find(abs(R.candidateDirectionsDeg)<1e-9,1,'first');
null=find(abs(abs(R.candidateDirectionsDeg)-180)<1e-9,1,'first');
row=struct();
row.dataset=string(R.dataset); row.probeOffsetDeg=R.offsetDeg;
row.sigmaMTDeg=R.sigmaMTDeg; row.sigmaReadoutDeg=R.sigmaReadoutDeg;
row.rectified=R.rectifyCandidates; row.nTrials=n;
row.nExcludedNonfinite=R.nExcludedNonfinite;
row.nCandidates=numel(R.candidateDirectionsDeg);
row.changePrimaryWinnerFraction=mean(R.change.winner==primary);
row.noChangePrimaryWinnerFraction=mean(R.noChange.winner==primary);
row.changeOpponentOnlyWinnerFraction=mean(R.isOpponentOnlyCandidate(R.change.winner));
row.noChangeOpponentOnlyWinnerFraction=mean(R.isOpponentOnlyCandidate(R.noChange.winner));
nonPrimaryPhysical=R.isPhysicalCandidate;
nonPrimaryPhysical(primary)=false;
row.changePhysicalNonPrimaryWinnerFraction=mean(nonPrimaryPhysical(R.change.winner));
row.noChangePhysicalNonPrimaryWinnerFraction=mean(nonPrimaryPhysical(R.noChange.winner));
row.changeWinnerChangedByNoise=mean(R.change.winnerChangedByNoise,'omitnan');
row.noChangeWinnerChangedByNoise=mean(R.noChange.winnerChangedByNoise,'omitnan');
row.changePatchSelectedFraction=mean(R.changePatchSelected);
row.meanPatchMargin=mean(R.patchMargin,'omitnan');
row.sdPatchMargin=std(R.patchMargin,'omitnan');
q=localQuantile(R.patchMargin,[.05 .5 .95]);
row.q05PatchMargin=q(1); row.medianPatchMargin=q(2); row.q95PatchMargin=q(3);
row.meanChangePrimaryEvidence=mean(R.change.evidence(:,primary),'omitnan');
row.meanNoChangePrimaryEvidence=mean(R.noChange.evidence(:,primary),'omitnan');
if isempty(null)
  row.meanChangeOpponentEvidence=nan; row.meanNoChangeOpponentEvidence=nan;
else
  row.meanChangeOpponentEvidence=mean(R.change.evidence(:,null),'omitnan');
  row.meanNoChangeOpponentEvidence=mean(R.noChange.evidence(:,null),'omitnan');
end
end

function T=summarizeWinners(R)
nCandidates=numel(R.candidateDirectionsDeg);
dataset=repmat(string(R.dataset),2*nCandidates,1);
probeOffsetDeg=repmat(R.offsetDeg,2*nCandidates,1);
sigmaMTDeg=repmat(R.sigmaMTDeg,2*nCandidates,1);
sigmaReadoutDeg=repmat(R.sigmaReadoutDeg,2*nCandidates,1);
rectified=repmat(R.rectifyCandidates,2*nCandidates,1);
patch=[repmat("change",nCandidates,1);repmat("noChange",nCandidates,1)];
candidateDirectionDeg=repmat(R.candidateDirectionsDeg,2,1);
roleLabel=repmat(R.roleLabels,2,1);
winnerFraction=nan(2*nCandidates,1);
signalWinnerFraction=nan(2*nCandidates,1);
for k=1:nCandidates
  winnerFraction(k)=mean(R.change.winner==k);
  winnerFraction(nCandidates+k)=mean(R.noChange.winner==k);
  signalWinnerFraction(k)=conditionalSignalWinner(R.change,k);
  signalWinnerFraction(nCandidates+k)=conditionalSignalWinner(R.noChange,k);
end
T=table(dataset,probeOffsetDeg,sigmaMTDeg,sigmaReadoutDeg,rectified,patch, ...
  candidateDirectionDeg,roleLabel,winnerFraction,signalWinnerFraction);
end

function f=conditionalSignalWinner(P,k)
use=P.signalWinnerUnique;
if any(use),f=mean(P.signalWinner(use)==k);else,f=nan;end
end

function T=summarizeNoiseSigns(R,C)
rows=cell(0,1);
for side=["change","noChange"]
  noise=C.(side+"Noise"); P=R.(side);
  for component=1:numel(C.componentDirectionsDeg)
    target=findDirection(R.candidateDirectionsDeg,C.componentDirectionsDeg(component));
    opponent=findDirection(R.candidateDirectionsDeg,C.componentDirectionsDeg(component)+180);
    for signValue=[-1 1]
      use=signValue.*noise(:,component)>0;
      row=struct('dataset',string(R.dataset),'probeOffsetDeg',R.offsetDeg, ...
        'sigmaMTDeg',R.sigmaMTDeg,'sigmaReadoutDeg',R.sigmaReadoutDeg, ...
        'rectified',R.rectifyCandidates,'patch',side, ...
        'componentRole',string(C.componentRoles(component)), ...
        'componentDirectionDeg',C.componentDirectionsDeg(component), ...
        'noiseSign',signValue,'nObservations',sum(use), ...
        'meanNoise',safeMean(noise(use,component)), ...
        'meanDeltaTargetEvidence',safeMean(P.evidence(use,target)-P.signalEvidence(use,target)), ...
        'meanDeltaOpponentEvidence',safeMean(P.evidence(use,opponent)-P.signalEvidence(use,opponent)), ...
        'meanDeltaTargetActivation',safeMean(P.activation(use,target)-P.signalActivation(use,target)), ...
        'meanDeltaOpponentActivation',safeMean(P.activation(use,opponent)-P.signalActivation(use,opponent)), ...
        'targetWinnerFraction',safeMean(P.winner(use)==target), ...
        'opponentWinnerFraction',safeMean(P.winner(use)==opponent));
      rows{end+1,1}=row; %#ok<AGROW>
    end
  end
end
T=struct2table(vertcat(rows{:}));
end

function row=findDirection(candidateDirections,direction)
d=mod(double(candidateDirections(:))-double(direction)+180,360)-180;
row=find(abs(d)<1e-9,1,'first');
if isempty(row)
  error('analyzeIDRIDQMTTrialCandidates:MissingCandidateDirection', ...
    'Candidate bank lacks direction %.6g.',direction);
end
end

function y=safeMean(x)
if isempty(x),y=nan;else,y=mean(double(x),'omitnan');end
end

function q=localQuantile(x,p)
x=sort(double(x(isfinite(x))));
if isempty(x),q=nan(size(p));return,end
idx=1+(numel(x)-1)*p; lo=floor(idx); hi=ceil(idx);
q=x(lo).*(hi-idx)+x(hi).*(idx-lo);
same=lo==hi; q(same)=x(lo(same));
end

function fig=makePlot(summary,selectedSigma,visible)
T=summary.diagnosticTable;
use=abs(T.sigmaReadoutDeg-selectedSigma)<1e-9;
T=T(use,:);
fig=figure('Color','w','Visible',visible,'Units','inches','Position',[1 1 10 7]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('Actual-trial deterministic MT candidates (sigma_M_T %.1f, sigma_R %.1f deg)', ...
  summary.sigmaMTDeg,selectedSigma));
metrics={'changePrimaryWinnerFraction','changeOpponentOnlyWinnerFraction', ...
  'changeWinnerChangedByNoise','changePatchSelectedFraction'};
titles={'Changed patch: primary wins','Changed patch: opponent-only wins', ...
  'Changed-patch winner altered by noise','Changed patch has larger maximum'};
for k=1:4
  ax=nexttile(tl);hold(ax,'on');
  for rect=[false true]
    R=T(T.rectified==rect,:); idr=R.dataset=="IDR"; idq=R.dataset=="IDQ";
    style='--'; marker='o'; if rect,style='-';marker='s';end
    plot(ax,R.probeOffsetDeg(idr),R.(metrics{k})(idr), ...
      'LineStyle',style,'Marker',marker,'LineWidth',1.2, ...
      'DisplayName',signedLabel(rect));
    plot(ax,R.probeOffsetDeg(idq),R.(metrics{k})(idq),'d', ...
      'MarkerSize',7,'LineWidth',1.2,'HandleVisibility','off');
  end
  xlabel(ax,'Offset (deg); diamond = IDQ 120 deg');ylabel(ax,'Fraction');
  title(ax,titles{k});ylim(ax,[0 1]);xlim(ax,[0 185]);grid(ax,'on');box(ax,'off');
  if k==1,legend(ax,'Location','best');end
end
end

function label=signedLabel(rect)
if rect,label='Rectified candidates';else,label='Signed candidates';end
end

function path=defaultSummaryPath(domain,file)
path=fullfile(domainFolder(mfilename('fullpath'),domain), ...
  'Data','AcrossSessionSummaries',file);
end

function tf=isPositiveScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0;
end
