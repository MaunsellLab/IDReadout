function M = plotIDRIDQMatchedSummary(varargin)
% plotIDRIDQMatchedSummary  Multipage matched IDQ-IDR descriptive report.
%
% Pages:
%   1. Head-to-head empirical dashboard
%   2. IDQ task-specific detail
%   3. IDR task-specific detail, using the same visual grammar
%   4. Side-by-side stimulus and analysis specification
%
% This function loads saved fits through makeIDRIDQMatchedSummary. It does
% not refit behavioral models.
%
% Name-value arguments:
%   IDQSummaryFolder     Folder containing IDQ across-session MAT files
%   IDRSummaryFolder     Folder containing IDR across-session MAT files
%   CommonSummaryFolder  Folder for the compact matched-summary MAT file
%   OutputPath     Output PDF path
%   Visible        Figure visibility (default 'off')


p=inputParser;
addParameter(p,'IDQSummaryFolder',defaultDomainSummaryFolder('IDQ'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'IDRSummaryFolder',defaultDomainSummaryFolder('IDR'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'CommonSummaryFolder',defaultDomainSummaryFolder('Common Code'), ...
  @(x) ischar(x)||isstring(x));
addParameter(p,'OutputPath','',@(x) ischar(x)||isstring(x));
addParameter(p,'Visible','off',@(x) any(strcmpi(string(x),["on","off"])));
parse(p,varargin{:});
opts=p.Results;

idqSummaryFolder=char(opts.IDQSummaryFolder);
idrSummaryFolder=char(opts.IDRSummaryFolder);
commonSummaryFolder=char(opts.CommonSummaryFolder);
if strlength(string(opts.OutputPath))==0
  outputPath=fullfile(defaultPlotFolder(),'IDRIDQ_MatchedSummary.pdf');
else
  outputPath=char(opts.OutputPath);
end

M=makeIDRIDQMatchedSummary('IDQSummaryFolder',idqSummaryFolder, ...
  'IDRSummaryFolder',idrSummaryFolder, ...
  'CommonSummaryFolder',commonSummaryFolder, ...
  'SavePath',fullfile(commonSummaryFolder,'IDRIDQ_MatchedSummary.mat'), ...
  'SaveSummary',true);

if isfile(outputPath), delete(outputPath); end
figures=gobjects(4,1);
figures(1)=plotDashboard(M,opts.Visible);
figures(2)=plotIDQDetail(M,opts.Visible);
figures(3)=plotIDRDetail(M,opts.Visible);
figures(4)=plotMethodsPage(M,opts.Visible);
for k=1:numel(figures)
  exportgraphics(figures(k),outputPath,'ContentType','vector', ...
    'Append',k>1);
  if strcmpi(opts.Visible,'off'), close(figures(k)); end
end
M.reportPath=string(outputPath);
fprintf('Saved matched IDQ-IDR report: %s\n',outputPath);
end

%% ------------------------------------------------------------------------
function folder=defaultDomainSummaryFolder(domainName)
try
  root=domainFolder(mfilename('fullpath'),domainName);
  folder=fullfile(root,'Data','AcrossSessionSummaries');
catch
  folder=pwd;
end
end

function folder=defaultPlotFolder()
try
  root=domainFolder(mfilename('fullpath'),'Common Code');
  folder=validFolder(fullfile(root,'Plots','AcrossSessionSummaries'));
catch
  folder=pwd;
end
end

function fig=newPage(visible)
fig=figure('Color','w','Visible',visible,'Units','inches', ...
  'Position',[0.5 0.5 11 8.5],'PaperPositionMode','auto');
end

function C=reportColors()
C=struct();
C.preferred=[0.48 0.25 0.72];
C.probe=[0.00 0.48 0.72];
C.positive=[0.82 0.35 0.25];
C.negative=[0.25 0.55 0.78];
C.neutral=[0.40 0.40 0.40];
C.sum=[0.35 0.60 0.82];
C.max=[0.82 0.45 0.35];
end

%% ------------------------------------------------------------------------
function fig=plotDashboard(M,visible)
C=reportColors(); fig=newPage(visible);
tl=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
title(tl,'Matched IDQ-IDR Empirical Summary','FontWeight','bold','FontSize',15);

psychAx=gobjects(2,1);
psychAx(1)=nexttile(tl,1,[1 2]); plotIDQPsychometric(psychAx(1),M.idq.psychometric);
title(psychAx(1),sprintf('IDQ psychometric (%d sessions)',M.idq.nSessions));
psychAx(2)=nexttile(tl,3,[1 2]); plotIDRPsychometric(psychAx(2),M.idr.psychometric);
title(psychAx(2),sprintf('IDR psychometric (%d sessions)',M.idr.nSessions));
xl=cell2mat(get(psychAx,'XLim')); yl=cell2mat(get(psychAx,'YLim'));
set(psychAx,'XLim',[min(xl(:,1)) max(xl(:,2))], ...
  'YLim',[min(yl(:,1)) max(yl(:,2))]);

gainAx=gobjects(4,1);
gainAx(1)=nexttile(tl,5); plotIDQLinearGains(gainAx(1),M.idq.linear.change,'IDQ change side',C);
gainAx(2)=nexttile(tl,6); plotIDQNoChangeSigned(gainAx(2),M.idq.signed.noChange,C);
gainAx(3)=nexttile(tl,7); plotIDRMatchedChange(gainAx(3),M.idr.matchedGains,C);
gainAx(4)=nexttile(tl,8); plotIDRNoChangeSigned(gainAx(4),M.idr.signedGains,C);
set(gainAx,'YLim',[-2 2]);

ax=nexttile(tl,9,[1 2]); plotIDQProfile(ax,M.idq.pNorm.controlled,C);
title(ax,'IDQ controlled no-change pooling profile');
ax=nexttile(tl,11,[1 2]); plotIDRProfile(ax,M.idr.pooling,C);
title(ax,'IDR controlled no-change pooling profile');
end

%% ------------------------------------------------------------------------
function fig=plotIDQDetail(M,visible)
C=reportColors(); Q=M.idq; fig=newPage(visible);
tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('IDQ Detail: Three Direction Mechanisms (%d sessions)',Q.nSessions), ...
  'FontWeight','bold','FontSize',14);

ax=nexttile(tl,1,[1 2]); plotIDQPsychometric(ax,Q.psychometric);
title(ax,'Aligned pooled psychometric');
kernelAx=gobjects(2,1);
kernelAx(1)=nexttile(tl,3); plotKernel(kernelAx(1),Q.kernel.change,'Change-side aggregate',C.preferred);
kernelAx(2)=nexttile(tl,4); plotKernel(kernelAx(2),Q.kernel.noChange,'No-change aggregate',C.probe);
yl=cell2mat(get(kernelAx,'YLim'));
set(kernelAx,'YLim',[min(yl(:,1)) max(yl(:,2))]);

ax=nexttile(tl,5); plotIDQLinearGains(ax,Q.linear.change,'Change linear gains',C);
ax=nexttile(tl,6); plotIDQLinearGains(ax,Q.linear.noChange,'No-change linear gains',C);
ax=nexttile(tl,7); plotIDQChangeSigned(ax,Q.signed.change,C);
ax=nexttile(tl,8); plotIDQNoChangeSigned(ax,Q.signed.noChange,C);

ax=nexttile(tl,9); plotIDQOverlap(ax,Q.overlap.change,'Change candidate margin',C);
ax=nexttile(tl,10); plotIDQOverlap(ax,Q.overlap.noChange,'No-change candidate margin',C);
ax=nexttile(tl,11); plotIDQSumMax(ax,Q.pooling,C);
ax=nexttile(tl,12); plotIDQProfile(ax,Q.pNorm.controlled,C);

ax=nexttile(tl,13); plotIDQAbsoluteGains(ax,Q.linear.change,C);
ax=nexttile(tl,14); plotIDQInteractions(ax,Q.interaction,C);
ax=nexttile(tl,15); plotIDQWinnerFractions(ax,Q.pooling.winnerStats,C);
ax=nexttile(tl,16); plotIDQText(ax,Q);
end

%% ------------------------------------------------------------------------
function fig=plotIDRDetail(M,visible)
C=reportColors(); D=M.idr; fig=newPage(visible);
tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('IDR Detail: Preferred and Probe Mechanisms (%d sessions)',D.nSessions), ...
  'FontWeight','bold','FontSize',14);

ax=nexttile(tl,1,[1 2]); plotIDRPsychometric(ax,D.psychometric);
title(ax,'Aligned pooled psychometric');
ax=nexttile(tl,3); plotIDRCommonKernel(ax,D.commonKernel,C.preferred);
ax=nexttile(tl,4); plotIDRPredictorDiagnostics(ax,D.predictorDiagnostics,C);

ax=nexttile(tl,5); plotIDRMatchedChange(ax,D.matchedGains,C);
ax=nexttile(tl,6); plotIDRMatchedNoChange(ax,D.matchedGains,C);
ax=nexttile(tl,7); plotIDRChangeSigned(ax,D.signedGains,C);
ax=nexttile(tl,8); plotIDRNoChangeSigned(ax,D.signedGains,C);

ax=nexttile(tl,9); plotIDROverlapQuantiles(ax,D.candidateOverlap,C);
ax=nexttile(tl,10); plotIDRCrossings(ax,D.candidateOverlap,C);
ax=nexttile(tl,11); plotIDRSumMaxByOffset(ax,D.pooling.offset,C);
ax=nexttile(tl,12); plotIDRProfile(ax,D.pooling,C);

ax=nexttile(tl,13); plotIDRProbeEvidence(ax,D.matchedGains,C);
ax=nexttile(tl,14); plotIDRPoolingWinners(ax,D.pooling.offset,C);
ax=nexttile(tl,15); plotIDRSideGains(ax,D.linearSide,C);
ax=nexttile(tl,16); plotIDRText(ax,D);
end

%% ------------------------------------------------------------------------
function plotIDQPsychometric(ax,P)
B=P.binned; F=P.fit; hold(ax,'on');
use=B.nTrials>0 & isfinite(B.alignedCoh) & isfinite(B.pCorrect);
plot(ax,B.alignedCoh(use),B.pCorrect(use),'ko','MarkerFaceColor','k','MarkerSize',4);
xMax=max([B.alignedCoh(use);F.threshold])*1.08;
x=linspace(0,xMax,300);
y=1-F.lapse-(0.5-F.lapse).*exp(-(max(x,0)./F.alpha).^F.betaWeibull);
plot(ax,x,y,'LineWidth',1.4); xline(ax,1,'k--');
yline(ax,F.thresholdPerformance,'k:');
xlabel(ax,'Coherence / session c_{75}'); ylabel(ax,'P(correct)');
ylim(ax,[0.48 1.01]); xlim(ax,[0 xMax]); grid(ax,'on'); box(ax,'off');
text(ax,.98,.03,sprintf('%d trials\n\\beta %.2f; lapse %.3f', ...
  P.nTrials,F.betaWeibull,F.lapse),'Units','normalized', ...
  'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',7);
end

function plotIDRPsychometric(ax,F)
B=F.binned; hold(ax,'on');
plot(ax,B.alignedCoh,B.pCorrect,'ko','MarkerFaceColor','k','MarkerSize',4);
xMax=max([B.alignedCoh;F.threshold])*1.08; x=linspace(0,xMax,300);
y=1-F.lapse-(0.5-F.lapse).*exp(-(max(x,0)./F.alpha).^F.betaWeibull);
plot(ax,x,y,'LineWidth',1.4); xline(ax,1,'k--');
yline(ax,F.thresholdPerformance,'k:');
xlabel(ax,'Coherence / session c_{75}'); ylabel(ax,'P(correct)');
ylim(ax,[0.48 1.01]); xlim(ax,[0 xMax]); grid(ax,'on'); box(ax,'off');
text(ax,.98,.03,sprintf('%d trials\n\\beta %.2f; lapse %.3f', ...
  F.nTrials,F.betaWeibull,F.lapse),'Units','normalized', ...
  'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',7);
end

function plotKernel(ax,K,titleText,color)
hold(ax,'on'); plot(ax,K.tMS,K.meanDiff,'Color',color,'LineWidth',1.2);
plot(ax,K.tMS,K.rectReference,'k-','LineWidth',1.0); yline(ax,0,'k:');
xlabel(ax,'Trial time (ms)'); ylabel(ax,'Kernel'); title(ax,titleText);
grid(ax,'on'); box(ax,'off');
end

function plotIDRCommonKernel(ax,K,color)
hold(ax,'on'); plot(ax,K.tMS,K.kernel,'Color',color,'LineWidth',1.2);
yline(ax,0,'k:'); xlabel(ax,'Trial time (ms)'); ylabel(ax,'Kernel');
title(ax,'Common change-minus-no-change kernel'); grid(ax,'on'); box(ax,'off');
end

%% ------------------------------------------------------------------------
function plotIDQLinearGains(ax,G,titleText,C)
F=G.driftNonDrift; values=F.gain(:); se=F.SE(:);
bar(ax,1:2,values,'FaceColor',C.preferred,'EdgeColor','k'); hold(ax,'on');
errorbar(ax,1:2,values,se,'k.','LineWidth',1.1); yline(ax,0,'k:');
set(ax,'XTick',1:2,'XTickLabel',{'Drift','Non-drift'});
ylabel(ax,'Absolute gain'); title(ax,titleText); grid(ax,'on'); box(ax,'off');
end

function plotIDQAbsoluteGains(ax,G,C)
F=G.absolute; v=F.gain(:); se=F.SE(:);
bar(ax,1:3,v,'FaceColor',C.probe,'EdgeColor','k'); hold(ax,'on');
errorbar(ax,1:3,v,se,'k.'); yline(ax,0,'k:');
set(ax,'XTick',1:3,'XTickLabel',{'0-119','120-239','240-359'});
xtickangle(ax,25); ylabel(ax,'Absolute gain'); title(ax,'Absolute direction gains');
grid(ax,'on'); box(ax,'off');
end

function plotIDQChangeSigned(ax,F,C)
D=F.driftSignedSlopes; N=F.nonDriftSignedSlopes;
values=[D.gain(:),N.gain(:)]; b=bar(ax,1:2,values,'grouped');
b(1).FaceColor=C.preferred; b(2).FaceColor=C.probe; hold(ax,'on'); yline(ax,0,'k:');
set(ax,'XTick',1:2,'XTickLabel',D.signName); ylabel(ax,'Signed gain');
title(ax,'Change signed gains'); legend(ax,{'Drift','Non-drift'},'FontSize',6);
grid(ax,'on'); box(ax,'off');
end

function plotIDQNoChangeSigned(ax,F,C)
S=F.signedSlopes; v=S.gain(:); se=0.5*(S.CI95High-S.CI95Low);
bar(ax,1:2,v,'FaceColor',C.probe,'EdgeColor','k'); hold(ax,'on');
errorbar(ax,1:2,v,se,'k.'); yline(ax,0,'k:');
set(ax,'XTick',1:2,'XTickLabel',S.signName); ylabel(ax,'Signed gain');
title(ax,'No-change signed gains'); grid(ax,'on'); box(ax,'off');
end

%% ------------------------------------------------------------------------
function plotIDRMatchedChange(ax,S,C)
hold(ax,'on');
errorbar(ax,S.probeDirDeg,S.gCP,S.seCP,'o-','Color',C.preferred, ...
  'MarkerFaceColor',C.preferred,'DisplayName','Preferred');
errorbar(ax,S.probeDirDeg,S.gCQ,S.seCQ,'o-','Color',C.probe, ...
  'MarkerFaceColor',C.probe,'DisplayName','Probe');
yline(ax,0,'k:','HandleVisibility','off'); xlabel(ax,'Probe offset (deg)'); ylabel(ax,'Absolute gain');
title(ax,'Change-side gains'); legend(ax,'FontSize',6,'Location','best');
grid(ax,'on'); box(ax,'off');
end

function plotIDRMatchedNoChange(ax,S,C)
hold(ax,'on');
errorbar(ax,S.probeDirDeg,S.gNP,S.seNP,'o-','Color',C.preferred, ...
  'MarkerFaceColor',C.preferred,'DisplayName','Preferred');
errorbar(ax,S.probeDirDeg,S.gNQ,S.seNQ,'o-','Color',C.probe, ...
  'MarkerFaceColor',C.probe,'DisplayName','Probe');
yline(ax,0,'k:','HandleVisibility','off'); xlabel(ax,'Probe offset (deg)'); ylabel(ax,'Absolute gain');
title(ax,'No-change-side gains'); legend(ax,'FontSize',6,'Location','best');
grid(ax,'on'); box(ax,'off');
end

function plotIDRChangeSigned(ax,S,C)
hold(ax,'on');
plot(ax,S.probeDirDeg,S.allGCPpos,'o-','Color',C.positive,'DisplayName','Pref +');
plot(ax,S.probeDirDeg,S.allGCPneg,'o--','Color',C.negative,'DisplayName','Pref -');
plot(ax,S.probeDirDeg,S.allGCQpos,'s-','Color',C.positive,'DisplayName','Probe +');
plot(ax,S.probeDirDeg,S.allGCQneg,'s--','Color',C.negative,'DisplayName','Probe -');
yline(ax,0,'k:','HandleVisibility','off'); xlabel(ax,'Probe offset (deg)'); ylabel(ax,'Signed gain');
title(ax,'Change signed gains'); legend(ax,'FontSize',5,'Location','best');
grid(ax,'on'); box(ax,'off');
end

function plotIDRNoChangeSigned(ax,S,C)
hold(ax,'on');
plot(ax,S.probeDirDeg,S.gNPpos,'o-','Color',C.positive,'DisplayName','Pref +');
plot(ax,S.probeDirDeg,S.gNPneg,'o--','Color',C.negative,'DisplayName','Pref -');
plot(ax,S.probeDirDeg,S.gNQpos,'s-','Color',C.positive,'DisplayName','Probe +');
plot(ax,S.probeDirDeg,S.gNQneg,'s--','Color',C.negative,'DisplayName','Probe -');
yline(ax,0,'k:','HandleVisibility','off'); xlabel(ax,'Probe offset (deg)'); ylabel(ax,'Signed gain');
title(ax,'No-change signed gains'); legend(ax,'FontSize',5,'Location','best');
grid(ax,'on'); box(ax,'off');
end

%% ------------------------------------------------------------------------
function plotIDQOverlap(ax,O,titleText,C)
H=O.histogram; bar(ax,H.centers,H.probability,1,'FaceColor',C.neutral, ...
  'EdgeColor','none'); hold(ax,'on'); xline(ax,0,'k--');
xlabel(ax,'Candidate margin (% coherence)'); ylabel(ax,'Probability');
title(ax,titleText); grid(ax,'on'); box(ax,'off');
text(ax,.03,.97,sprintf('cross %.4f\nq_{99} %.2f', ...
  O.crossingProbability,O.q99),'Units','normalized','VerticalAlignment','top', ...
  'FontSize',7);
end

function plotIDROverlapQuantiles(ax,S,C)
hold(ax,'on');
plot(ax,S.probeDirDeg,S.q95ChangeMargin,'o-','Color',C.negative,'DisplayName','q95');
plot(ax,S.probeDirDeg,S.q99ChangeMargin,'s-','Color',C.positive,'DisplayName','q99');
plot(ax,S.probeDirDeg,S.maxChangeMargin,'^-','Color',C.neutral,'DisplayName','max');
yline(ax,0,'k:','HandleVisibility','off'); xlabel(ax,'Probe offset (deg)'); ylabel(ax,'Margin (% coherence)');
title(ax,'Change candidate margin'); legend(ax,'FontSize',6,'Location','best');
grid(ax,'on'); box(ax,'off');
end

function plotIDRCrossings(ax,S,C)
semilogy(ax,S.probeDirDeg,max(S.pChangeCrossed,1e-6),'o-', ...
  'Color',C.positive,'MarkerFaceColor',C.positive); hold(ax,'on');
semilogy(ax,S.probeDirDeg,max(S.pNoChangeCrossed,1e-6),'s-', ...
  'Color',C.negative,'MarkerFaceColor',C.negative);
xlabel(ax,'Probe offset (deg)'); ylabel(ax,'Crossing probability');
title(ax,'Physical candidate crossings'); legend(ax,{'Change','No-change'},'FontSize',6);
grid(ax,'on'); box(ax,'off');
end

%% ------------------------------------------------------------------------
function plotIDQSumMax(ax,F,C)
v=[F.controlledSum.negLogLikelihood,F.controlledMax.negLogLikelihood];
v=v-min(v); b=bar(ax,1:2,v); b.FaceColor='flat'; b.CData=[C.sum;C.max];
set(ax,'XTick',1:2,'XTickLabel',{'Sum','Hard max'}); ylabel(ax,'NLL - best');
title(ax,sprintf('Controlled pooling; max over sum %.2f', ...
  F.controlledContrast.deltaNLLMaxOverSum)); grid(ax,'on'); box(ax,'off');
end

function plotIDQProfile(ax,F,C)
P=F.profile; hold(ax,'on');
use=isfinite(P.p) & P.p>=1 & isfinite(P.deltaNLL);
semilogx(ax,P.p(use),P.deltaNLL(use),'o-','Color',C.max,'MarkerFaceColor',C.max);
set(ax,'XScale','log');
if any(use), xlim(ax,[1 max(P.p(use))]); end
yline(ax,0.5*3.84145882069412,'k:'); xlabel(ax,'Fixed p'); ylabel(ax,'Profile \DeltaNLL');
grid(ax,'on'); box(ax,'off');
text(ax,.97,.95,sprintf('95%% lower %.2f',P.profileCI95Low), ...
  'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
end

function plotIDRSumMaxByOffset(ax,S,C)
x=1:height(S);
bar(ax,x,S.deltaNLLSumMinusMax,'FaceColor',C.max,'EdgeColor','k');
set(ax,'XTick',x,'XTickLabel',compose('%g',S.probeDirDeg)); xtickangle(ax,45);
yline(ax,0,'k:'); xlabel(ax,'Probe offset (deg)'); ylabel(ax,'NLL(sum)-NLL(max)');
title(ax,'Hard max versus sum'); grid(ax,'on'); box(ax,'off');
end

function plotIDRProfile(ax,F,C)
P=F.profile; use=isfinite(P.p) & P.p>=1 & isfinite(P.deltaNLL);
semilogx(ax,P.p(use),P.deltaNLL(use),'o-','Color',C.max, ...
  'MarkerFaceColor',C.max); hold(ax,'on');
set(ax,'XScale','log');
if any(use), xlim(ax,[1 max(P.p(use))]); end
yline(ax,1.920729,'k:'); xlabel(ax,'Fixed shared p'); ylabel(ax,'Profile \DeltaNLL');
grid(ax,'on'); box(ax,'off');
text(ax,.97,.95,sprintf('1-sided 95%% lower %.2f',F.profileLowerOneSided95), ...
  'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
end

%% ------------------------------------------------------------------------
function plotIDQInteractions(ax,F,C)
v=[F.driftNonDriftInteraction.deltaNLL,F.nonDriftInteraction.deltaNLL, ...
  F.bothInteractions.deltaNLL];
bar(ax,1:3,v,'FaceColor',C.neutral,'EdgeColor','k');
set(ax,'XTick',1:3,'XTickLabel',{'D x N','N x N','Both'});
ylabel(ax,'\DeltaNLL'); title(ax,'Change-side interactions');
grid(ax,'on'); box(ax,'off');
end

function plotIDQWinnerFractions(ax,W,C)
v=[W.changeDriftMatchedFraction,W.changeNonDriftFraction; ...
  W.noChangePreferredFraction,W.noChangeOpponentFraction];
b=bar(ax,v,'stacked'); b(1).FaceColor=C.preferred; b(2).FaceColor=C.probe;
set(ax,'XTick',1:2,'XTickLabel',{'Change','No-change'}); ylim(ax,[0 1]);
ylabel(ax,'Fraction'); title(ax,'Model-implied max candidate');
legend(ax,{'Preferred','Other'},'FontSize',6); grid(ax,'on'); box(ax,'off');
end

function plotIDRPredictorDiagnostics(ax,T,C)
row=strcmp(string(T.predictorType),'rectStep');
if ~any(row), axis(ax,'off'); title(ax,'Predictor diagnostics unavailable'); return; end
v=[T.sdChange(row),T.sdNoChange(row)];
bar(ax,1:2,v,'FaceColor',C.neutral,'EdgeColor','k');
set(ax,'XTick',1:2,'XTickLabel',{'Change','No-change'});
ylabel(ax,'Predictor SD (% coherence)'); title(ax,'Preferred predictor scale');
grid(ax,'on'); box(ax,'off');
end

function plotIDRProbeEvidence(ax,S,C)
x=1:height(S);
bar(ax,x,S.deltaNLLProbeTerms,'FaceColor',C.probe,'EdgeColor','k');
set(ax,'XTick',x,'XTickLabel',compose('%g',S.probeDirDeg)); xtickangle(ax,45);
xlabel(ax,'Probe offset (deg)'); ylabel(ax,'\DeltaNLL');
title(ax,'Evidence for probe terms'); grid(ax,'on'); box(ax,'off');
end

function plotIDRPoolingWinners(ax,S,C)
plot(ax,S.probeDirDeg,S.maxProbeWinnerFraction,'o-','Color',C.probe, ...
  'MarkerFaceColor',C.probe); ylim(ax,[0 1]);
xlabel(ax,'Probe offset (deg)'); ylabel(ax,'Probe winner fraction');
title(ax,'Model-implied hard-max winner'); grid(ax,'on'); box(ax,'off');
end

function plotIDRSideGains(ax,F,C)
U=F.unrestricted; v=U.gain(:); se=U.SE(:);
bar(ax,1:2,v,'FaceColor',C.preferred,'EdgeColor','k'); hold(ax,'on');
errorbar(ax,1:2,v,se,'k.'); yline(ax,0,'k:');
set(ax,'XTick',1:2,'XTickLabel',{'Change','No-change'});
ylabel(ax,'Gain'); title(ax,'Preferred-stream side gains'); grid(ax,'on'); box(ax,'off');
end

%% ------------------------------------------------------------------------
function plotIDQText(ax,Q)
axis(ax,'off'); C=Q.pooling.controlledContrast;
txt={sprintf('%d sessions; %d noisy trials',Q.nSessions,Q.nStepNoiseTrials); ...
  sprintf('Predictor: %s',Q.predictorType); ...
  sprintf('Mean step %.2f%%',Q.meanStepCoh); ...
  sprintf('Stream SD mean %.2f%%',mean(Q.stepNoiseSD)); ...
  sprintf('g_D %.3f; g_N %.3f',Q.linear.change.driftNonDrift.gain); ...
  sprintf('Sum-max \\DeltaNLL %.3f',C.deltaNLLMaxOverSum); ...
  'Winner identity is model-implied'};
text(ax,0,1,txt,'Units','normalized','VerticalAlignment','top', ...
  'FontName','Menlo','FontSize',7);
end

function plotIDRText(ax,D)
axis(ax,'off'); O=D.pooling.overall;
txt={sprintf('%d sessions; %d noisy trials',D.nSessions,D.nStepNoiseTrials); ...
  sprintf('Predictor: %s',D.predictorType); ...
  sprintf('Excluded offset: %g deg',D.excludedProbeDirectionsDeg); ...
  sprintf('Shared p %.2f',D.pooling.sharedP); ...
  sprintf('Sum-max \\DeltaNLL %.3f',O.deltaNLLSumMinusMax); ...
  sprintf('1-sided lower p %.2f',D.pooling.profileLowerOneSided95); ...
  'Winner identity is model-implied'};
text(ax,0,1,txt,'Units','normalized','VerticalAlignment','top', ...
  'FontName','Menlo','FontSize',7);
end

%% ------------------------------------------------------------------------
function fig=plotMethodsPage(M,visible)
fig=newPage(visible);
tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
title(tl,'IDQ-IDR Stimulus and Analysis Contract','FontWeight','bold','FontSize',15);

ax=nexttile(tl,1,[2 1]); plotMethodColumn(ax,'IDQ',M.analysisContract.idq, ...
  [0.96 0.94 1.00]);
ax=nexttile(tl,2,[2 1]); plotMethodColumn(ax,'IDR',M.analysisContract.idr, ...
  [0.93 0.97 1.00]);
ax=nexttile(tl,5,[1 2]); axis(ax,'off');
common=[{'COMMON ANALYSIS CONTRACT'};M.analysisContract.common;{''; ...
  'Interpretive boundary: physical-coherence overlap is not neural candidate overlap.'; ...
  'A shared forward model must add MT tuning, divisive normalization, internal variability,'; ...
  'within-patch selection, and comparison of patch outputs.'}];
text(ax,.02,.95,common,'Units','normalized','VerticalAlignment','top', ...
  'FontName','Menlo','FontSize',9,'Interpreter','none');
end

function plotMethodColumn(ax,heading,items,bg)
axis(ax,'off'); rectangle(ax,'Position',[0 0 1 1],'FaceColor',bg, ...
  'EdgeColor',[.75 .75 .75]);
lines=cell(numel(items)*2+1,1); lines{1}=heading; j=2;
for k=1:numel(items)
  lines{j}=sprintf('%02d  %s',k,items{k}); j=j+1;
  lines{j}=''; j=j+1;
end
text(ax,.04,.96,lines,'Units','normalized','VerticalAlignment','top', ...
  'FontName','Menlo','FontSize',8.5,'Interpreter','none');
end
