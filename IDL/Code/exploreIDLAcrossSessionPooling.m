function [result, fig] = exploreIDLAcrossSessionPooling(varargin)
% exploreIDLAcrossSessionPooling  Meetz-only, across-session sigma validation.
% TwoAllocationMode=false selects the original style; true selects 100:0/50:50.
%
% result = exploreIDLAcrossSessionPooling()
% result = exploreIDLAcrossSessionPooling('SessionFolder',folder)
%
% For each held-out session: estimate the 100:0 correct-minus-error temporal
% kernel on all other Meetz sessions, integrate the framewise stimulus with
% that kernel, fit a separate L1/L2/Max sigma on ALL allocations in those
% other sessions, and score every trial in the held-out session. No choice
% from the test session contributes to either the kernel or sigma.
%
% Per-session all-allocation sigmas and psychometric thresholds are included
% ONLY as descriptive diagnostics. They are not used in held-out predictions.
% The divided-only scores permit comparison with the existing 100:0 transfer
% analysis; the all-allocation score answers the new prediction question.
% Nothing here changes or overwrites ModelRecovery.mat.

p = inputParser;
p.KeepUnmatched = true;
p.PartialMatching = false;
p.addParameter('TwoAllocationMode',false,@(x) islogical(x)&&isscalar(x));
p.addParameter('SessionFolder','',@(x) ischar(x) || isstring(x));
p.addParameter('SaveResult',true,@(x) islogical(x) && isscalar(x));
p.addParameter('ExportPDF',true,@(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
root = domainFolder(mfilename('fullpath'));
folder = char(string(p.Results.SessionFolder));
if isempty(folder), folder = fullfile(root,'Data','SessionAnalysis'); end
assert(isfolder(folder),'exploreIDLAcrossSessionPooling:MissingFolder', ...
  'SessionAnalysis folder not found: %s',folder);
selection=p.Unmatched;
if isfield(selection,'Animal')
  assert(strcmpi(string(selection.Animal),'Meetz'), ...
    'exploreIDLAcrossSessionPooling:Animal','This exploration is Meetz-only.');
  selection=rmfield(selection,'Animal');
end
names = fieldnames(selection);
args = cell(1,2*numel(names));
for k=1:numel(names)
  args{2*k-1}=names{k}; args{2*k}=selection.(names{k});
end
paths = selectAnalysisFiles(folder,'Animal','Meetz', ...
  'TwoAllocationMode',p.Results.TwoAllocationMode,args{:});
if ischar(paths), paths={paths}; end
if isstring(paths), paths=cellstr(paths(:)); end
paths=paths(:);
nSessions=numel(paths);
assert(nSessions>=3,'exploreIDLAcrossSessionPooling:TooFewSessions', ...
  'At least three Meetz sessions are needed for this analysis.');
A=cell(nSessions,1);
for s=1:nSessions
  S=load(paths{s},'sessionAnalysis');
  assert(isfield(S,'sessionAnalysis'), ...
    'exploreIDLAcrossSessionPooling:MissingAnalysis', ...
    '%s lacks sessionAnalysis.',paths{s});
  A{s}=S.sessionAnalysis;
  idlAssertAnalysisMode(A{s},p.Results.TwoAllocationMode);
  validateSession(A{s});
end
nFrames=numel(A{1}.stepTimeMS);
assert(all(cellfun(@(x) isequal(double(x.stepTimeMS(:)), ...
  double(A{1}.stepTimeMS(:))),A)), ...
  'exploreIDLAcrossSessionPooling:TimeGrid', ...
  'All sessions must share a video-frame time grid.');

modelNames=["L1" "L2" "Max"];
ids=strings(nSessions,1);
nTrials=zeros(nSessions,1);
nDivided=zeros(nSessions,1);
sigmaTraining=nan(nSessions,3);
sigmaSessionDiagnostic=nan(nSessions,3);
threshold100Diagnostic=nan(nSessions,1);
llAll=nan(nSessions,3);
llDivided=nan(nSessions,3);
llReference=nan(nSessions,3);
looKernel=nan(nFrames,nSessions);
predictions=cell(nSessions,1);

for held=1:nSessions
  test=A{held};
  ids(held)=string(test.sourceSessionID);
  keep=setdiff(1:nSessions,held);
  w=trainingKernel(A(keep));
  looKernel(:,held)=w;
  trainEvidence=cell(1,3);
  trainChoice=cell(numel(keep),1);
  for j=1:numel(keep)
    trainChoice{j}=double(A{keep(j)}.trialTable.correct(:));
  end
  yTrain=vertcat(trainChoice{:});
  yTest=double(test.trialTable.correct(:));
  divided=abs(double(test.trialTable.unsignedAllocPC)-100)>1e-6;
  nTrials(held)=numel(yTest);
  nDivided(held)=sum(divided);
  probs=nan(nTrials(held),3);
  for m=1:3
    X=cell(numel(keep),1);
    for j=1:numel(keep)
      X{j}=poolEvidence(A{keep(j)},w,modelNames(m));
    end
    trainEvidence{m}=vertcat(X{:});
    sigmaTraining(held,m)=fitSigma(trainEvidence{m},yTrain);
    xTest=poolEvidence(test,w,modelNames(m));
    probs(:,m)=normalCDF(xTest/sigmaTraining(held,m));
    llAll(held,m)=bernoulliLL(probs(:,m),yTest);
    llDivided(held,m)=bernoulliLL(probs(divided,m),yTest(divided));
    llReference(held,m)=bernoulliLL(probs(~divided,m),yTest(~divided));
    % Descriptive fit to held-out outcomes; NEVER used to score them.
    sigmaSessionDiagnostic(held,m)=fitSigma(xTest,yTest);
  end
  T=test.trialTable;
  fit=fitIDLSessionPsychometrics(T);
  threshold100Diagnostic(held)=fit.thresholdTable.threshold(1);
  predictions{held}=table(repmat(ids(held),nTrials(held),1), ...
    double(T.sourceTrialNumber(:)),logical(T.correct(:)), ...
    double(T.unsignedAllocPC(:)),double(T.totalStepCohPC(:)), ...
    probs(:,1),probs(:,2),probs(:,3), ...
    'VariableNames',{'sourceSessionID','sourceTrialNumber','correct', ...
    'unsignedAllocPC','totalStepCohPC','pL1','pL2','pMax'});
end

sessionTable=table(ids,nTrials,nDivided,threshold100Diagnostic, ...
  sigmaTraining(:,1),sigmaTraining(:,2),sigmaTraining(:,3), ...
  sigmaSessionDiagnostic(:,1),sigmaSessionDiagnostic(:,2), ...
  sigmaSessionDiagnostic(:,3), ...
  llAll(:,3)-llAll(:,2),llAll(:,3)-llAll(:,1), ...
  llDivided(:,3)-llDivided(:,2),llDivided(:,3)-llDivided(:,1), ...
  'VariableNames',{'sourceSessionID','nTrials','nDividedTrials', ...
  'threshold100Diagnostic','trainSigmaL1','trainSigmaL2','trainSigmaMax', ...
  'sessionSigmaL1Diagnostic','sessionSigmaL2Diagnostic', ...
  'sessionSigmaMaxDiagnostic','allMaxMinusL2','allMaxMinusL1', ...
  'dividedMaxMinusL2','dividedMaxMinusL1'});
summaryTable=makeSummary(sessionTable);
result=struct();
result.twoAllocationMode=p.Results.TwoAllocationMode;
tag=idlAllocationTag(p.Results.TwoAllocationMode);
result.createdAt=datetime('now');
result.animal="Meetz";
result.method=['Held-out whole sessions; one animal-wide sigma per model ' ...
  'fit on other sessions/all allocations; kernel from other sessions/100:0'];
result.selection=selection;
result.sourceFiles=string(paths);
result.modelNames=modelNames;
result.looKernel=looKernel;
result.stepTimeMS=double(A{1}.stepTimeMS(:));
result.sessionTable=sessionTable;
result.summaryTable=summaryTable;
result.predictionTable=vertcat(predictions{:});
result.heldOutLLAll=llAll;
result.heldOutLLDivided=llDivided;
result.heldOutLL100=llReference;
result.trainSigma=sigmaTraining;
result.sessionSigmaDiagnostic=sigmaSessionDiagnostic;
result.threshold100Diagnostic=threshold100Diagnostic;
result.note=['The session-specific sigma and threshold diagnostics use the ' ...
  'held-out outcomes and must not be used to interpret the held-out LL.'];
if p.Results.SaveResult
  dest=validFolder(fullfile(root,'Data','AcrossSessions'));
  save(fullfile(dest,['AcrossSessionPooling_Meetz' tag '.mat']),'result','-v7.3');
end
fig=plotModelSummary(result);
diagnosticFig=plotDiagnostics(result);
if p.Results.ExportPDF
  plotFolder=validFolder(fullfile(root,'Plots','AcrossSessions'));
  exportgraphics(fig,fullfile(plotFolder, ...
    ['AcrossSessionPooling_Meetz' tag '.pdf']),'ContentType','vector');
  exportgraphics(diagnosticFig,fullfile(plotFolder, ...
    ['AcrossSessionPoolingDiagnostics_Meetz' tag '.pdf']),'ContentType','vector');
end
disp(summaryTable)
end

function fig=plotModelSummary(R)
% Match the 3-by-3 layout of fitIDLModelRecoveryAcrossSessions. The top
% comparisons use the divided allocations, while the summary also reports
% the all-allocation predictive comparison.
colors=[.10 .45 .80; .90 .45 .10; .20 .60 .35];
T=R.sessionTable;
fig=figure('Name','IDL model recovery | Meetz | shared sigma', ...
  'NumberTitle','off','Color','w','Position',[100 100 1600 980]);
fig.WindowStyle='docked';
tl=tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
sgtitle(tl,'IDL model recovery: Meetz (across-session sigma; LOO kernel)', ...
  'FontWeight','bold');

ax=nexttile(tl,1); axis(ax,'off');
info=sprintf(['Meetz\n\nSessions: %d\nAnalysis trials: %d\n' ...
  '100:0 trials: %d\nDivided trials: %d\n\nPrepared: %s\n\n' ...
  'Kernel: other sessions, 100:0\n' ...
  'Sigma: other sessions, all allocations\n' ...
  'Plots/tests: held-out divided trials'], ...
  height(T),sum(T.nTrials),sum(T.nTrials-T.nDividedTrials), ...
  sum(T.nDividedTrials),datestr(R.createdAt,'yyyy-mm-dd HH:MM'));
text(ax,.02,.98,info,'Units','normalized', ...
  'VerticalAlignment','top','Interpreter','none','FontSize',10);

ax=nexttile(tl,2);
d=[T.dividedMaxMinusL1 T.dividedMaxMinusL2];
plotAdvantages(ax,T.sourceSessionID,d,false);
ax=nexttile(tl,3);
plotAdvantages(ax,T.sourceSessionID,d,true);

ax=nexttile(tl,6);
plot(ax,R.stepTimeMS,R.looKernel,'-', ...
  'Color',[.75 .75 .75],'LineWidth',.8); hold(ax,'on');
plot(ax,R.stepTimeMS,mean(R.looKernel,2),'k-', ...
  'LineWidth',2);
xlabel(ax,'Time from step onset (ms)');
ylabel(ax,'Normalized weight');
title(ax,'Temporal kernels (other sessions, 100:0)');
grid(ax,'on'); box(ax,'off');

allocations=[100 75 62.5 50];
tiles=[4 5 7 8];
P=R.predictionTable;
stepMax=max(5,5*ceil(max(P.totalStepCohPC)/5));
for a=1:4
  ax=nexttile(tl,tiles(a)); hold(ax,'on');
  select=abs(P.unsignedAllocPC-allocations(a))<1e-6;
  stepValues=unique(P.totalStepCohPC(select));
  x=zeros(numel(stepValues),1);
  observed=x; predicted=zeros(numel(x),3);
  for j=1:numel(stepValues)
    rows=select & abs(P.totalStepCohPC-stepValues(j))<1e-6;
    x(j)=mean(P.totalStepCohPC(rows));
    observed(j)=mean(P.correct(rows));
    predicted(j,:)=mean(P{rows,{'pL1','pL2','pMax'}},1);
  end
  plot(ax,x,observed,'ko','MarkerFaceColor','k', ...
    'DisplayName','Observed');
  for m=1:3
    plot(ax,x,predicted(:,m),'-','Color',colors(m,:), ...
      'LineWidth',1.5+.5*(m==3),'DisplayName',char(R.modelNames(m)));
  end
  yline(ax,.75,':','HandleVisibility','off');
  xlabel(ax,'Total coherence step (%)'); ylabel(ax,'P(correct)');
  xlim(ax,[0 stepMax]); xticks(ax,0:5:stepMax);
  ylim(ax,[.45 1.02]); grid(ax,'on'); box(ax,'off');
  title(ax,sprintf('%g:%g (held out)', ...
    allocations(a),100-allocations(a)));
  if a==1, legend(ax,'Location','southeast'); end
end

ax=nexttile(tl,9); axis(ax,'off');
G=R.summaryTable;
textLines=sprintf('Divided allocations (primary comparison)\n\n');
for k=3:4
  textLines=[textLines sprintf(['%s: %d/%d positive\n' ...
    'Mean %.2f LL/1000; 95%% CI [%.2f, %.2f]\n' ...
    'One-sided p = %.3g; total %+.1f LL\n\n'], ...
    char(G.comparison(k)),G.nPositive(k),G.nSessions(k), ...
    G.meanPer1000(k),G.lower95(k),G.upper95(k), ...
    G.oneSidedP(k),G.totalDeltaLL(k))]; %#ok<AGROW>
end
textLines=[textLines sprintf(['All allocations (including 100:0)\n' ...
  'Max - L2: %+.1f LL; Max - L1: %+.1f LL\n\n' ...
  'All choices held out by session.'], ...
  G.totalDeltaLL(1),G.totalDeltaLL(2))]; %#ok<AGROW>
text(ax,.02,.98,textLines,'Units','normalized', ...
  'VerticalAlignment','top','Interpreter','none','FontSize',9);
end

function plotAdvantages(ax,ids,d,cumulative)
n=size(d,1); x=1:n;
if cumulative
  values=cumsum(d,1);
  ylabel(ax,'Cumulative held-out LL advantage');
  title(ax,sprintf('Cumulative: %+.1f over L1, %+.1f over L2', ...
    values(end,1),values(end,2)));
else
  values=d;
  ylabel(ax,'Session held-out LL advantage');
  title(ax,'Individual sessions (divided allocations)');
end
plot(ax,x,values(:,1),'o-','LineWidth',1.2, ...
  'DisplayName','Max - L1'); hold(ax,'on');
plot(ax,x,values(:,2),'s-','LineWidth',1.2, ...
  'DisplayName','Max - L2');
yline(ax,0,':','HandleVisibility','off');
labels=string(ids);
for k=1:n
  date=regexp(char(labels(k)),'[0-9]{8}$','match','once');
  if ~isempty(date)
    labels(k)=string(date(5:6))+"/"+string(date(7:8));
  end
end
ticks=1:max(1,ceil(n/12)):n;
set(ax,'XTick',ticks,'XTickLabel',labels(ticks));
xtickangle(ax,45); xlabel(ax,'Session date');
legend(ax,'Location','best'); grid(ax,'on'); box(ax,'off');
end

function fig=plotDiagnostics(R)
T=R.sessionTable;
fig=figure('Name','IDL across-session pooling | Meetz', ...
  'NumberTitle','off','Color','w');
fig.WindowStyle='docked';
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
x=1:height(T);
ax=nexttile(tl,1);
plot(ax,x,1000*T.allMaxMinusL2./T.nTrials,'o-', ...
  'LineWidth',1.4,'DisplayName','Max - L2'); hold(ax,'on');
plot(ax,x,1000*T.allMaxMinusL1./T.nTrials,'s-', ...
  'LineWidth',1.4,'DisplayName','Max - L1');
yline(ax,0,':'); grid(ax,'on');
ylabel(ax,'Held-out ΔLL / 1000 trials');
title(ax,'All allocations'); legend(ax,'Location','best');
ax=nexttile(tl,2);
plot(ax,x,cumsum(T.allMaxMinusL2),'o-', ...
  'LineWidth',1.4,'DisplayName','Max - L2'); hold(ax,'on');
plot(ax,x,cumsum(T.allMaxMinusL1),'s-', ...
  'LineWidth',1.4,'DisplayName','Max - L1');
yline(ax,0,':'); grid(ax,'on');
ylabel(ax,'Cumulative held-out ΔLL');
title(ax,'Session order in selected file list');
ax=nexttile(tl,3);
plot(ax,x,T.trainSigmaMax,'o-','LineWidth',1.4, ...
  'DisplayName','Other-session sigma'); hold(ax,'on');
plot(ax,x,T.sessionSigmaMaxDiagnostic,'s-', ...
  'LineWidth',1.4,'DisplayName','Own-session sigma (diagnostic)');
grid(ax,'on'); ylabel(ax,'Max sigma (% coherence)');
xlabel(ax,'Held-out session'); legend(ax,'Location','best');
ax=nexttile(tl,4);
scatter(ax,T.threshold100Diagnostic, ...
  T.sessionSigmaMaxDiagnostic./T.trainSigmaMax,40,'filled');
yline(ax,1,':'); grid(ax,'on');
xlabel(ax,'Own-session 100:0 threshold (%)');
ylabel(ax,'Own / other-session Max sigma');
title(ax,'Sensitivity variation (descriptive)');
sgtitle(tl,sprintf('Meetz | %d sessions | shared sigma trained on other sessions', ...
  height(T)));
end

function validateSession(A)
required={'sessionHeader','sourceSessionID','trialTable','stepTimeMS', ...
  'changeNoiseByFrameTrial','noChangeNoiseByFrameTrial'};
assert(all(isfield(A,required)), ...
  'exploreIDLAcrossSessionPooling:MissingFields', ...
  'Session analysis lacks fields needed for model evidence.');
assert(strcmpi(string(A.sessionHeader.animal),'Meetz'), ...
  'exploreIDLAcrossSessionPooling:Animal','Only Meetz sessions are supported.');
n=height(A.trialTable);
f=numel(A.stepTimeMS);
assert(isequal(size(A.changeNoiseByFrameTrial),[f n 2]) && ...
  isequal(size(A.noChangeNoiseByFrameTrial),[f n 2]), ...
  'exploreIDLAcrossSessionPooling:NoiseShape', ...
  'Step noise arrays do not match the trial table.');
end

function w=trainingKernel(train)
f=numel(train{1}.stepTimeMS);
correctSum=zeros(f,1); errorSum=zeros(f,1);
nCorrect=0; nError=0;
for s=1:numel(train)
  A=train{s}; T=A.trialTable;
  rows=find(abs(double(T.unsignedAllocPC)-100)<1e-6);
  X=zeros(f,numel(rows));
  for j=1:numel(rows)
    stream=double(T.highStreamIndex(rows(j)));
    assert(ismember(stream,[1 2]), ...
      'exploreIDLAcrossSessionPooling:ReferenceStream', ...
      'Invalid highStreamIndex on a 100:0 trial.');
    X(:,j)=double(A.changeNoiseByFrameTrial(:,rows(j),stream));
  end
  good=logical(T.correct(rows));
  correctSum=correctSum+sum(X(:,good),2);
  errorSum=errorSum+sum(X(:,~good),2);
  nCorrect=nCorrect+sum(good); nError=nError+sum(~good);
end
assert(nCorrect>0 && nError>0, ...
  'exploreIDLAcrossSessionPooling:KernelOutcomes', ...
  'Training sessions need correct and error 100:0 trials.');
w=correctSum/nCorrect-errorSum/nError;
if sum(w)<0, w=-w; end
assert(sum(abs(w))>0, ...
  'exploreIDLAcrossSessionPooling:ZeroKernel','Temporal kernel has zero mass.');
w=w/sum(abs(w));
end

function evidence=poolEvidence(A,w,name)
T=A.trialTable;
n=height(T);
allocation=[double(T.alloc0(:)) double(T.alloc1(:))];
step=allocation.*double(T.totalStepCohPC(:));
pedestal=double(T.prestepCohPC(:));
change=zeros(n,2); noChange=zeros(n,2);
for j=1:2
  C=double(A.changeNoiseByFrameTrial(:,:,j))+ ...
    reshape(pedestal+step(:,j),1,n);
  N=double(A.noChangeNoiseByFrameTrial(:,:,j))+reshape(pedestal,1,n);
  change(:,j)=C'*w;
  noChange(:,j)=N'*w;
end
switch name
  case "L1"
    C=sum(abs(change),2); N=sum(abs(noChange),2);
  case "L2"
    C=sqrt(sum(change.^2,2)); N=sqrt(sum(noChange.^2,2));
  case "Max"
    C=max(abs(change),[],2); N=max(abs(noChange),[],2);
end
evidence=C-N;
end

function sigma=fitSigma(x,y)
scale=max(std(x),1);
objective=@(z) -bernoulliLL(normalCDF(x/exp(z)),y);
bounds=log([scale*1e-3 scale*1e3]);
z=fminbnd(objective,bounds(1),bounds(2), ...
  optimset('Display','off','TolX',1e-8));
sigma=exp(z);
end

function p=normalCDF(x)
p=.5*erfc(-double(x)/sqrt(2));
p=min(max(p,1e-12),1-1e-12);
end

function ll=bernoulliLL(p,y)
p=min(max(double(p(:)),1e-12),1-1e-12);
y=double(y(:));
ll=sum(y.*log(p)+(1-y).*log(1-p));
end

function Q=makeSummary(T)
comparison=["Max - L2";"Max - L1";"Max - L2";"Max - L1"];
scoreSet=["All";"All";"Divided";"Divided"];
nSessions=repmat(height(T),4,1);
nScoredTrials=zeros(4,1); nPositive=zeros(4,1);
meanPer1000=zeros(4,1); sdPer1000=zeros(4,1);
lower95=zeros(4,1); upper95=zeros(4,1); oneSidedP=zeros(4,1);
totalDeltaLL=zeros(4,1); pooledPer1000=zeros(4,1);
fields={'allMaxMinusL2','allMaxMinusL1', ...
  'dividedMaxMinusL2','dividedMaxMinusL1'};
for k=1:4
  d=double(T.(fields{k}));
  if k<=2, n=double(T.nTrials); else, n=double(T.nDividedTrials); end
  rate=1000*d./n;
  nScoredTrials(k)=sum(n); nPositive(k)=sum(d>0);
  meanPer1000(k)=mean(rate);
  sdPer1000(k)=std(rate);
  se=sdPer1000(k)/sqrt(height(T));
  df=height(T)-1;
  if se>0
    critical=fzero(@(z) tUpper(z,df)-.025,[0 1000]);
    lower95(k)=meanPer1000(k)-critical*se;
    upper95(k)=meanPer1000(k)+critical*se;
    oneSidedP(k)=tUpper(meanPer1000(k)/se,df);
  else
    lower95(k)=meanPer1000(k); upper95(k)=meanPer1000(k);
    oneSidedP(k)=double(meanPer1000(k)<0)+.5*double(meanPer1000(k)==0);
  end
  totalDeltaLL(k)=sum(d);
  pooledPer1000(k)=1000*sum(d)/sum(n);
end
Q=table(scoreSet,comparison,nSessions,nScoredTrials,nPositive, ...
  meanPer1000,sdPer1000,lower95,upper95,oneSidedP, ...
  totalDeltaLL,pooledPer1000);
end

function p=tUpper(t,df)
q=.5*betainc(df/(df+t.^2),df/2,.5);
if t>=0, p=q; else, p=1-q; end
end
