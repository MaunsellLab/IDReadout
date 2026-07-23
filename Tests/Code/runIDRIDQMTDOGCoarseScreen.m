function screen=runIDRIDQMTDOGCoarseScreen(varargin)
% runIDRIDQMTDOGCoarseScreen  Resumable expected-choice DOG grid screen.
%
% Only IDR gCQ(theta) enters model ranking. Other IDR gains and performance
% are descriptive held-out consequences. Full replicated IDQ/IDR recovery is
% deferred until a small number of finalists has been selected.

p=inputParser;
addParameter(p,'RunIDPrefix',"",@(x)ischar(x)||isstring(x));
addParameter(p,'SigmaCenterDeg',5,@positiveScalar);
addParameter(p,'SigmaSurroundDeg',[20 40 80 160],@positiveVector);
addParameter(p,'SurroundMagnitude',[.05 .15 .30 .60],@nonnegativeVector);
addParameter(p,'NChoiceDraws',20,@positiveInteger);
addParameter(p,'Seed',1729,@validSeed);
addParameter(p,'SigmaMTDeg',37.5,@positiveScalar);
addParameter(p,'CandidateSpacingDeg',2,@positiveScalar);
addParameter(p,'InputSpacingDeg',1,@positiveScalar);
addParameter(p,'IDQPreStepPedestalPC',7,@nonnegativeScalar);
addParameter(p,'TargetPerformance',.75,@validTarget);
addParameter(p,'NMonteCarlo',2000,@positiveInteger);
addParameter(p,'RectifyCandidates',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'MaxRunBytes',2e8,@positiveScalar);
addParameter(p,'IDQSummaryPath',defaultSummaryPath('IDQ','IDQ_AcrossSideSummary.mat'), ...
  @(x)ischar(x)||isstring(x));
addParameter(p,'IDRSummaryPath',defaultSummaryPath('IDR','IDR_SideGainAnalysis.mat'), ...
  @(x)ischar(x)||isstring(x));
addParameter(p,'Visible','off',@(x)any(strcmpi(string(x),["on","off"])));
parse(p,varargin{:});opts=p.Results;

prefix=string(opts.RunIDPrefix);
if ~isscalar(prefix)||strlength(prefix)==0
  error('runIDRIDQMTDOGCoarseScreen:RunIDPrefixRequired', ...
    'Supply a unique RunIDPrefix.');
end
sigmaS=unique(double(opts.SigmaSurroundDeg(:)'),'stable');
strength=unique(double(opts.SurroundMagnitude(:)'),'stable');
if any(sigmaS<=opts.SigmaCenterDeg)
  error('runIDRIDQMTDOGCoarseScreen:SurroundNotBroader', ...
    'Every surround sigma must exceed SigmaCenterDeg.');
end
strength=strength(strength>0);
sourcePaths=string({opts.IDQSummaryPath;opts.IDRSummaryPath});
if ~all(isfile(sourcePaths))
  error('runIDRIDQMTDOGCoarseScreen:MissingSource','Source summary file missing.');
end

grid=makeGrid(opts.SigmaCenterDeg,sigmaS,strength);
members=repmat(emptyMember(),height(grid),1);
for k=1:height(grid)
  sigmaSurround=grid.sigmaSurroundDeg(k);a=grid.surroundMagnitude(k);
  runID=formatDOGRunID(prefix,opts.SigmaCenterDeg,sigmaSurround,a);
  fprintf('\n=== DOG member %d/%d: %s ===\n',k,height(grid),runID);
  context=makeAnalysisContext("synthetic",'RunID',runID,'Seed',opts.Seed, ...
    'MaxRunBytes',opts.MaxRunBytes,'SaveLatents',false);
  parameters=memberParameters(opts,sigmaSurround,a);
  ensureMemberRun(context,parameters,sourcePaths);

  calibrationPath=analysisPath(context,'Common Code','Data','Simulation', ...
    'IDRIDQ_MTUpstreamNoiseCalibration.mat');
  if isfile(calibrationPath)
    A=load(calibrationPath,'calibration');calibration=A.calibration;
    verifyStage(calibration,context,a,sigmaSurround,'calibration');
    fprintf('Loaded completed calibration.\n');
  else
    calibration=calibrateIDRIDQMTUpstreamNoise(context, ...
      'IDQSummaryPath',sourcePaths(1),'IDRSummaryPath',sourcePaths(2), ...
      'SigmaMTDeg',opts.SigmaMTDeg,'SigmaReadoutDeg',opts.SigmaCenterDeg, ...
      'SurroundMagnitude',a,'SigmaSurroundDeg',sigmaSurround, ...
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
    A=load(choicePath,'simulation');simulation=A.simulation;
    verifyStage(simulation,context,a,sigmaSurround,'simulation');
    if simulation.nReplicates~=opts.NChoiceDraws
      error('runIDRIDQMTDOGCoarseScreen:ChoiceDrawMismatch', ...
        'Existing member has a different number of choice draws.');
    end
    fprintf('Loaded completed choices.\n');
  else
    simulation=simulateIDRIDQUpstreamNoiseChoices(context,calibration, ...
      'NReplicates',opts.NChoiceDraws, ...
      'RectifyCandidates',opts.RectifyCandidates,'Visible','off');
  end

  expected=fitIDRExpectedChoiceDOGGains(context, ...
    'IDRSummaryPath',sourcePaths(2), ...
    'TargetPerformance',opts.TargetPerformance);
  expectedPath=analysisPath(context,'Common Code','Data','Simulation', ...
    'IDRIDQ_MTDOGExpectedGains.mat');
  members(k)=struct('sigmaCenterDeg',opts.SigmaCenterDeg, ...
    'sigmaSurroundDeg',sigmaSurround,'surroundMagnitude',a, ...
    'runID',runID,'context',context,'choicePath',string(choicePath), ...
    'expectedPath',string(expectedPath),'gCQChiSquare',expected.gCQChiSquare);
end
screen=makeOrLoadSummary(prefix,members,sourcePaths,opts);
end

function grid=makeGrid(sigmaCenter,sigmaSurround,strength)
sigmaCenterDeg=sigmaCenter;sigmaSurroundDeg=nan;surroundMagnitude=0;
grid=table(sigmaCenterDeg,sigmaSurroundDeg,surroundMagnitude);
[S,A]=ndgrid(sigmaSurround,strength);
extra=table(repmat(sigmaCenter,numel(S),1),S(:),A(:), ...
  'VariableNames',grid.Properties.VariableNames);
grid=[grid;extra];
end

function P=memberParameters(O,sigmaSurround,a)
P=struct('sigmaMTDeg',double(O.SigmaMTDeg), ...
  'sigmaCenterDeg',double(O.SigmaCenterDeg), ...
  'sigmaSurroundDeg',double(sigmaSurround), ...
  'surroundMagnitude',double(a), ...
  'candidateSpacingDeg',double(O.CandidateSpacingDeg), ...
  'inputSpacingDeg',double(O.InputSpacingDeg), ...
  'IDQPreStepPedestalPC',double(O.IDQPreStepPedestalPC), ...
  'targetPerformance',double(O.TargetPerformance), ...
  'nMonteCarlo',double(O.NMonteCarlo), ...
  'nChoiceDraws',double(O.NChoiceDraws), ...
  'rectifyCandidates',logical(O.RectifyCandidates),'poolingMode',"hardMax", ...
  'selectionTarget',"IDR gCQ(theta) only");
end

function ensureMemberRun(context,parameters,sources)
manifestPath=analysisPath(context,'Tests','manifest.mat');
if isfile(manifestPath)
  M=validateSyntheticManifest(context);
  if string(M.ModelName)~="IDRIDQMTDOGCoarseScreenMember"|| ...
      ~isequaln(M.ModelParameters,parameters)
    error('runIDRIDQMTDOGCoarseScreen:MemberManifestMismatch', ...
      'Existing member manifest has different parameters.');
  end
  validateSimulationSourceFiles(context,sources);
  fprintf('Resuming existing member run.\n');
else
  initializeSyntheticRun(context,'ModelName',"IDRIDQMTDOGCoarseScreenMember", ...
    'ModelParameters',parameters,'SourceFiles',sources);
end
end

function verifyStage(S,context,a,sigmaSurround,label)
if ~isstruct(S)||~isfield(S,'context')||string(S.context.RunID)~=context.RunID
  error('runIDRIDQMTDOGCoarseScreen:StageContextMismatch', ...
    'Stored %s belongs to another run.',label);
end
if ~isfield(S,'surroundMagnitude')||S.surroundMagnitude~=a|| ...
    ~isfield(S,'sigmaSurroundDeg')||~isequaln(S.sigmaSurroundDeg,sigmaSurround)
  error('runIDRIDQMTDOGCoarseScreen:StageParameterMismatch', ...
    'Stored %s has different DOG parameters.',label);
end
end

function screen=makeOrLoadSummary(prefix,members,experimentalSources,opts)
runID=prefix+"_summary";
if strlength(runID)>80,error('runIDRIDQMTDOGCoarseScreen:SummaryIDTooLong', ...
    'Summary RunID exceeds 80 characters.');end
context=makeAnalysisContext("synthetic",'RunID',runID,'Seed',opts.Seed, ...
  'MaxRunBytes',1e8,'SaveLatents',false);
sources=strings(2*numel(members),1);
for k=1:numel(members)
  sources(2*k-1)=members(k).choicePath;sources(2*k)=members(k).expectedPath;
end
parameters=struct('memberRunIDs',string({members.runID}), ...
  'sigmaCenterDeg',double(opts.SigmaCenterDeg), ...
  'sigmaSurroundDeg',double(opts.SigmaSurroundDeg(:)'), ...
  'surroundMagnitude',double(opts.SurroundMagnitude(:)'), ...
  'nChoiceDraws',double(opts.NChoiceDraws), ...
  'selectionTarget',"IDR gCQ(theta) only");
manifestPath=analysisPath(context,'Tests','manifest.mat');
if isfile(manifestPath)
  M=validateSyntheticManifest(context);
  if string(M.ModelName)~="IDRIDQMTDOGCoarseScreenSummary"|| ...
      ~isequaln(M.ModelParameters,parameters)
    error('runIDRIDQMTDOGCoarseScreen:SummaryManifestMismatch', ...
      'Existing summary manifest has different parameters.');
  end
  validateSimulationSourceFiles(context,sources);
else
  initializeSyntheticRun(context,'ModelName',"IDRIDQMTDOGCoarseScreenSummary", ...
    'ModelParameters',parameters,'SourceFiles',sources);
end
finalPath=analysisPath(context,'Common Code','Data','Simulation', ...
  'IDRIDQ_MTDOGCoarseScreen.mat');
if isfile(finalPath)
  A=load(finalPath,'screen');screen=A.screen;
  screen=saveOutputs(screen,opts.Visible);
  fprintf('Loaded completed DOG coarse screen: %s\n',finalPath);return
end
screen=buildSummary(context,members,sources,experimentalSources,opts);
screen=saveOutputs(screen,opts.Visible);
end

function screen=buildSummary(context,members,sources,experimentalSources,opts)
n=numel(members);gain=cell(n,1);performance=cell(n,1);metrics=cell(n,1);
sigmaCenterDeg=nan(n,1);sigmaSurroundDeg=nan(n,1);surroundMagnitude=nan(n,1);
gCQChiSquare=nan(n,1);runID=strings(n,1);
for k=1:n
  E=load(members(k).expectedPath,'expected');C=load(members(k).choicePath,'simulation');
  sigmaCenterDeg(k)=members(k).sigmaCenterDeg;
  sigmaSurroundDeg(k)=members(k).sigmaSurroundDeg;
  surroundMagnitude(k)=members(k).surroundMagnitude;
  gCQChiSquare(k)=E.expected.gCQChiSquare;runID(k)=members(k).runID;
  gain{k}=addParameters(E.expected.summaryTable,members(k));
  performance{k}=addParameters(C.simulation.conditionSummary,members(k));
  metrics{k}=addParameters(C.simulation.noiseMetricTable,members(k));
end
[bestValue,bestIndex]=min(gCQChiSquare);deltaChiSquare=gCQChiSquare-bestValue;
isBest=(1:n)'==bestIndex;
rankingTable=table(sigmaCenterDeg,sigmaSurroundDeg,surroundMagnitude, ...
  gCQChiSquare,deltaChiSquare,isBest,runID);
rankingTable=sortrows(rankingTable,'gCQChiSquare');
screen=struct();screen.createdAt=datetime('now');screen.createdBy=mfilename;
screen.dataOrigin="synthetic_expected_choice_DOG_screen";screen.context=context;
screen.memberRuns=members;screen.sourcePaths=sources;
screen.experimentalSourcePaths=experimentalSources;
screen.nChoiceDraws=opts.NChoiceDraws;screen.fitTarget="IDR gCQ(theta) only";
screen.rankingTable=rankingTable;screen.bestMember=members(bestIndex);
screen.idrGainTable=vertcat(gain{:});screen.performanceTable=vertcat(performance{:});
screen.noiseMetricTable=vertcat(metrics{:});
screen.heldOutDefinition=["gCP, gNP, gNQ, performance, and all IDQ gains " + ...
  "do not enter coarse model selection. IDQ gains are deferred to full " + ...
  "replicated recovery for selected finalists."];
end

function T=addParameters(T,M)
T.sigmaCenterDeg=repmat(M.sigmaCenterDeg,height(T),1);
T.sigmaSurroundDeg=repmat(M.sigmaSurroundDeg,height(T),1);
T.surroundMagnitude=repmat(M.surroundMagnitude,height(T),1);
T=movevars(T,{'sigmaCenterDeg','sigmaSurroundDeg','surroundMagnitude'},'Before',1);
end

function screen=saveOutputs(screen,visible)
context=screen.context;dataFolder=analysisPath(context,'Common Code','Data','Simulation');
plotFolder=analysisPath(context,'Common Code','Plots','Simulation');
if ~isfolder(dataFolder),mkdir(dataFolder);end
if ~isfolder(plotFolder),mkdir(plotFolder);end
matPath=fullfile(dataFolder,'IDRIDQ_MTDOGCoarseScreen.mat');
rankingPath=fullfile(dataFolder,'IDRIDQ_MTDOGCoarseScreen_Ranking.csv');
gainPath=fullfile(dataFolder,'IDRIDQ_MTDOGCoarseScreen_IDRGains.csv');
performancePath=fullfile(dataFolder,'IDRIDQ_MTDOGCoarseScreen_Performance.csv');
metricPath=fullfile(dataFolder,'IDRIDQ_MTDOGCoarseScreen_NoiseMetrics.csv');
pdfPath=fullfile(plotFolder,'IDRIDQ_MTDOGCoarseScreen.pdf');
screen.outputPaths=string({matPath;rankingPath;gainPath;performancePath;metricPath;pdfPath});
info=whos('screen');assertSimulationRunCapacity(context,info.bytes);
if ~isfile(rankingPath),writeTableAtomic(screen.rankingTable,rankingPath);end
if ~isfile(gainPath),writeTableAtomic(screen.idrGainTable,gainPath);end
if ~isfile(performancePath),writeTableAtomic(screen.performanceTable,performancePath);end
if ~isfile(metricPath),writeTableAtomic(screen.noiseMetricTable,metricPath);end
if ~isfile(pdfPath)
  figs=makePlots(screen,visible);temporary=pdfPath+".temporary.pdf";
  exportgraphics(figs(1),temporary,'ContentType','vector');
  exportgraphics(figs(2),temporary,'ContentType','vector','Append',true);
  exportgraphics(figs(3),temporary,'ContentType','vector','Append',true);
  if strcmpi(visible,'off'),close(figs);end;commitTemporary(temporary,pdfPath);
end
if ~isfile(matPath)
  temporary=matPath+".temporary.mat";save(temporary,'screen','-v7.3');
  commitTemporary(temporary,matPath);
end
assertSimulationRunCapacity(context,0);
end

function figs=makePlots(S,visible)
R=S.rankingTable;G=S.idrGainTable;figs=gobjects(3,1);
figs(1)=figure('Color','w','Visible',visible,'Units','inches','Position',[1 1 10 5]);
tl=tiledlayout(figs(1),1,2,'TileSpacing','compact','Padding','compact');
title(tl,'DOG coarse screen: selection uses IDR gCQ only');
ax=nexttile(tl);hold(ax,'on');widths=unique(R.sigmaSurroundDeg(isfinite(R.sigmaSurroundDeg)));
colors=lines(numel(widths));
for k=1:numel(widths)
  A=sortrows(R(R.sigmaSurroundDeg==widths(k),:),'surroundMagnitude');
  plot(ax,A.surroundMagnitude,A.gCQChiSquare,'o-','Color',colors(k,:), ...
    'LineWidth',1.2,'DisplayName',sprintf('sigma_S = %g deg',widths(k)));
end
baseline=R(R.surroundMagnitude==0,:);yline(ax,baseline.gCQChiSquare,'k--', ...
  'DisplayName','Gaussian baseline');xlabel(ax,'Surround magnitude');
ylabel(ax,'gCQ weighted squared error');grid(ax,'on');box(ax,'off');legend(ax,'Location','best');
ax=nexttile(tl);axis(ax,'off');best=R(1,:);
text(ax,.02,.92,sprintf([ ...
  'Best model\ncenter sigma: %.3g deg\nsurround sigma: %.3g deg\n' ...
  'surround magnitude: %.3g\ngCQ criterion: %.3f\nbaseline criterion: %.3f\n' ...
  'Delta baseline-best: %.3f'],best.sigmaCenterDeg,best.sigmaSurroundDeg, ...
  best.surroundMagnitude,best.gCQChiSquare,baseline.gCQChiSquare, ...
  baseline.gCQChiSquare-best.gCQChiSquare),'VerticalAlignment','top','FontSize',11);

figs(2)=figure('Color','w','Visible',visible,'Units','inches','Position',[1 1 10 7]);
tl=tiledlayout(figs(2),2,2,'TileSpacing','compact','Padding','compact');
title(tl,'Best DOG: fitted target and held-out IDR gains');
names={'gCP','gCQ','gNP','gNQ'};titles={'Change preferred','Change probe (fit target)', ...
  'No-change preferred (held out)','No-change probe (held out)'};
bestUse=G.surroundMagnitude==best.surroundMagnitude& ...
  sameNaN(G.sigmaSurroundDeg,best.sigmaSurroundDeg);
baseUse=G.surroundMagnitude==0;B=sortrows(G(bestUse,:),'probeDirDeg');
Z=sortrows(G(baseUse,:),'probeDirDeg');
for j=1:4
  ax=nexttile(tl);hold(ax,'on');name=names{j};
  plot(ax,Z.probeDirDeg,Z.(['synthetic' name]),'o--','Color',[.5 .5 .5], ...
    'DisplayName','Gaussian baseline');
  plot(ax,B.probeDirDeg,B.(['synthetic' name]),'o-','LineWidth',1.4, ...
    'DisplayName','Best DOG');
  errorbar(ax,B.probeDirDeg,B.(['observed' name]),B.(['observedSE' name]), ...
    'ks-','MarkerFaceColor','k','DisplayName','Observed');yline(ax,0,'k:');
  xlabel(ax,'Probe offset (deg)');ylabel(ax,'Gain');title(ax,titles{j});
  grid(ax,'on');box(ax,'off');if j==1,legend(ax,'Location','best');end
end

figs(3)=figure('Color','w','Visible',visible,'Units','inches','Position',[1 1 11 6]);
tl=tiledlayout(figs(3),2,2,'TileSpacing','compact','Padding','compact');
title(tl,'DOG gCQ families (the only selection target)');
for k=1:numel(widths)
  ax=nexttile(tl);hold(ax,'on');A=G(G.sigmaSurroundDeg==widths(k),:);
  strengths=unique(A.surroundMagnitude);c=lines(numel(strengths));
  for j=1:numel(strengths)
    Q=sortrows(A(A.surroundMagnitude==strengths(j),:),'probeDirDeg');
    plot(ax,Q.probeDirDeg,Q.syntheticgCQ,'o-','Color',c(j,:), ...
      'DisplayName',sprintf('a = %.3g',strengths(j)));
  end
  Q=sortrows(A(A.surroundMagnitude==strengths(1),:),'probeDirDeg');
  errorbar(ax,Q.probeDirDeg,Q.observedgCQ,Q.observedSEgCQ,'ks-', ...
    'MarkerFaceColor','k','DisplayName','Observed');yline(ax,0,'k:');
  title(ax,sprintf('Surround sigma = %g deg',widths(k)));xlabel(ax,'Offset (deg)');
  ylabel(ax,'gCQ');grid(ax,'on');box(ax,'off');if k==1,legend(ax,'Location','best');end
end
end

function tf=sameNaN(x,y)
tf=(x==y)|(isnan(x)&isnan(y));
end
function writeTableAtomic(T,path)
temporary=path+".temporary.csv";writetable(T,temporary);commitTemporary(temporary,path);
end
function commitTemporary(temporary,final)
[ok,message]=movefile(temporary,final,'f');if ~ok,error( ...
  'runIDRIDQMTDOGCoarseScreen:OutputCommitFailed','Could not commit output: %s',message);end
end
function M=emptyMember()
M=struct('sigmaCenterDeg',nan,'sigmaSurroundDeg',nan,'surroundMagnitude',nan, ...
  'runID',"",'context',struct(),'choicePath',"",'expectedPath',"", ...
  'gCQChiSquare',nan);
end
function path=defaultSummaryPath(domain,file)
path=fullfile(domainFolder(mfilename('fullpath'),domain), ...
  'Data','AcrossSessionSummaries',file);
end
function tf=positiveScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0;
end
function tf=nonnegativeScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0;
end
function tf=positiveVector(x)
tf=isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>0);
end
function tf=nonnegativeVector(x)
tf=isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0);
end
function tf=positiveInteger(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&fix(x)==x;
end
function tf=validSeed(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=2^32-1&&fix(x)==x;
end
function tf=validTarget(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>.5&&x<1;
end
