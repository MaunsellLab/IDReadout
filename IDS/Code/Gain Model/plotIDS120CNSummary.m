function summaries = plotIDS120CNSummary(varargin)
% plotIDS120CNSummary  Analyze C/N step-noise gains in 120-degree IDQ/IDS.
%
% The three nonoverlapping datasets are IDS, randomized-orientation IDQ,
% and fixed-orientation IDQ.  File exclusions remain under the control of
% selectAnalysisFiles.  Kernels span prestep and step for display, whereas
% regression gains use step-period noise only.

p = inputParser;
p.addParameter('Animal','All',@(x) ischar(x) || isstring(x));
p.addParameter('Dataset','All',@(x) ischar(x) || isstring(x));
p.addParameter('NBootstrap',1000,@(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('ExportPDF',true,@(x) islogical(x) && isscalar(x));
p.addParameter('ReportExcluded',false,@(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;

root = domainFolder(mfilename('fullpath'));
analysisFolder = fullfile(root,'Data','Probe120','SessionAnalysis');
if ~isfolder(analysisFolder)
  error('plotIDS120CNSummary:MissingAnalysisFolder', ...
    'Session-analysis folder not found: %s',analysisFolder);
end

specs = datasetSpecifications(opts.Dataset);
summaries = cell(numel(specs),1);
for k = 1:numel(specs)
  args = {'FilePattern','*_sessionAnalysis.mat','Animal',opts.Animal, ...
    'IncludeIDQ',true,'SourceDomain',specs(k).sourceDomain, ...
    'ReportExcluded',opts.ReportExcluded};
  if ~isempty(specs(k).randomBaseDir)
    args = [args {'RandomBaseDir',specs(k).randomBaseDir}]; %#ok<AGROW>
  end
  dataFilePaths = selectAnalysisFiles(analysisFolder,args{:});
  if isempty(dataFilePaths)
    warning('plotIDS120CNSummary:NoFiles', ...
      'No files selected for %s (%s).',specs(k).label,char(opts.Animal));
    summaries{k} = [];
    continue
  end
  analyses = loadAnalyses(dataFilePaths);
  summaries{k} = summarizeDataset(root,specs(k),analyses,char(opts.Animal), ...
    round(opts.NBootstrap),opts.ExportPDF);
end
summaries = summaries(~cellfun(@isempty,summaries));
if isscalar(summaries)
  summaries = summaries{1};
end
end

function specs = datasetSpecifications(requested)
allSpecs = struct( ...
  'name',{'IDS','IDQRandom','IDQFixed'}, ...
  'label',{'IDS','IDQ randomized','IDQ fixed'}, ...
  'sourceDomain',{'IDS','IDQ','IDQ'}, ...
  'randomBaseDir',{[],true,false});
requested = char(requested);
if strcmpi(requested,'All')
  specs = allSpecs;
  return
end
match = find(strcmpi({allSpecs.name},requested) | strcmpi({allSpecs.label},requested));
if numel(match) ~= 1
  error('plotIDS120CNSummary:BadDataset', ...
    'Dataset must be All, IDS, IDQRandom, or IDQFixed.');
end
specs = allSpecs(match);
end

function analyses = loadAnalyses(paths)
analyses = cell(numel(paths),1);
for k = 1:numel(paths)
  S = load(paths{k},'sessionAnalysis');
  if ~isfield(S,'sessionAnalysis')
    error('plotIDS120CNSummary:MissingAnalysis', ...
      '%s lacks sessionAnalysis.',paths{k});
  end
  validateAnalysis(S.sessionAnalysis,paths{k});
  analyses{k} = S.sessionAnalysis;
end
end

function validateAnalysis(A,path)
required = {'trialTable','tMS','stepFrames','sessionHeader','sessionProbeHeader'};
missing = required(~isfield(A,required));
if ~isempty(missing)
  error('plotIDS120CNSummary:IncompleteAnalysis', ...
    '%s lacks required field(s): %s.',path,strjoin(missing,', '));
end
changeNoise = getStreamNoise(A,'change');
noChangeNoise = getStreamNoise(A,'noChange');
if size(changeNoise,3) ~= 3 || size(noChangeNoise,3) ~= 3
  error('plotIDS120CNSummary:BadStreamCount', ...
    '%s does not contain three streams per RDK.',path);
end
if A.sessionProbeHeader.probeOffsetDeg ~= 120
  error('plotIDS120CNSummary:WrongOffset','%s is not a 120-degree analysis.',path);
end
end

function S = summarizeDataset(root,spec,analyses,animal,nBootstrap,exportPDF)
T = concatenateTables(analyses);
[psych,T] = fitPsychometric(T);
P = computeLOOStepPredictors(analyses);
predictorNames = {'noisePredCDrift','noisePredCPlus','noisePredCMinus', ...
  'noisePredNDrift','noisePredNPlus','noisePredNMinus'};
for k = 1:numel(predictorNames)
  T.(predictorNames{k}) = P(:,k);
end
gain = fitIDS120CNGain(T,psych.sessionFits,psych.alignedWeibull,0.75, ...
  'NBootstrap',nBootstrap,'BootstrapSeed',3100+datasetSeed(spec.name));
if nBootstrap > 0 && gain.sideRole.bootstrapNSuccess < 0.9*nBootstrap
  warning('plotIDS120CNSummary:GainBootstrapFailures', ...
    'Only %d of %d gain bootstrap refits succeeded for %s.', ...
    gain.sideRole.bootstrapNSuccess,nBootstrap,spec.label);
end

S = struct();
S.dataset = spec.name;
S.datasetLabel = spec.label;
S.sourceDomain = spec.sourceDomain;
S.randomBaseDir = spec.randomBaseDir;
S.animal = animal;
S.probeOffsetDeg = 120;
S.createdBy = mfilename;
S.createdAt = datetime('now');
S.nSessions = numel(analyses);
S.nTrials = height(T);
S.nStepNoiseTrials = sum(T.hasStepNoise);
S.nBootstrap = nBootstrap;
S.noiseAmplitudesPC = unique(round(T.cohNoisePC(isfinite(T.cohNoisePC)),6))';
S.trialTable = T;
S.psychFit = psych;
S.kernels = computeKernelSet(analyses,nBootstrap);
S.noiseGain = gain;
S.gainTable = makeGainSummaryTable(gain);
S.noChangeProbeInteractionTable = makeInteractionSummaryTable(gain);
S.sessionRecords = makeSessionRecords(analyses,psych.sessionFits);

summaryFolder = validFolder(fullfile(root,'Data','Probe120','AcrossSessionSummaries'));
plotFolder = validFolder(fullfile(root,'Plots','Probe120'));
tag = sprintf('IDS120CN_%s_%s',spec.name,matlab.lang.makeValidName(animal));
acrossSummary = S; %#ok<NASGU>
save(fullfile(summaryFolder,[tag '_AcrossSessionSummary.mat']),'acrossSummary','-v7.3');
writetable(S.gainTable,fullfile(summaryFolder,[tag '_GainSummary.csv']));
writetable(S.noChangeProbeInteractionTable, ...
  fullfile(summaryFolder,[tag '_NoChangeProbeInteraction.csv']));
plotSummary(S);
if exportPDF
  exportFig = figure('Color','w','WindowStyle','normal','Visible','off', ...
    'Units','inches','Position',[1 1 8 10.5]);
  cleanup = onCleanup(@() delete(exportFig)); %#ok<NASGU>
  plotSummary(S,exportFig);
  drawnow;
  exportgraphics(exportFig,fullfile(plotFolder,[tag '_Summary.pdf']), ...
    'ContentType','vector','Units','inches','Width',8,'Height',10.5, ...
    'Padding',0.1,'PreserveAspectRatio','on');
end
end

function seed = datasetSeed(name)
switch char(name)
  case 'IDS'
    seed = 1;
  case 'IDQRandom'
    seed = 2;
  case 'IDQFixed'
    seed = 3;
  otherwise
    seed = 0;
end
end

function T = concatenateTables(analyses)
parts = cell(numel(analyses),1);
for k = 1:numel(analyses)
  q = analyses{k}.trialTable;
  q.sessionIndex = repmat(k,height(q),1);
  q.noisyStepCoh = repmat(analyses{k}.noisyStepCoh,height(q),1);
  parts{k} = q;
end
T = vertcat(parts{:});
T = movevars(T,'sessionIndex','Before',1);
end

function [psych,T] = fitPsychometric(T)
target = 0.75;
initial = fitIDQInitialSessionThresholds(T,2,0.02,target);
aligned = fitIDQAcrossAlignedWeibull(initial.alignedCoh,T.correct,target,[0 0.05]);
T.alignedCoh = initial.alignedCoh;
T.noisyStepAlignedCoh = initial.noisyStepAlignedCoh;
psych = struct('title','All Directions','targetPerformance',target, ...
  'sessionFits',initial.sessionFits,'alignedWeibull',aligned);
end

function P = computeLOOStepPredictors(analyses)
nSessions = numel(analyses);
if nSessions < 2
  error('plotIDS120CNSummary:TooFewSessions', ...
    'Leave-one-session-out weighting requires at least two sessions.');
end
parts = cell(nSessions,1);
for held = 1:nSessions
  A = analyses{held};
  correctNoise = [];
  errorNoise = [];
  for source = 1:nSessions
    if source == held
      continue
    end
    B = analyses{source};
    validateTimeBase(A,B);
    use = logical(B.trialTable.hasStepNoise);
    changeSum = sum(getStreamNoise(B,'change'),3,'omitnan');
    correctNoise = [correctNoise changeSum(:,use & B.trialTable.correct)]; %#ok<AGROW>
    errorNoise = [errorNoise changeSum(:,use & ~B.trialTable.correct)]; %#ok<AGROW>
  end
  kernel = mean(correctNoise,2,'omitnan')-mean(errorNoise,2,'omitnan');
  weights = kernel(A.stepFrames);
  weightSum = sum(weights,'omitnan');
  if ~isfinite(weightSum) || abs(weightSum) < 100*eps(max(1,sum(abs(weights),'omitnan')))
    error('plotIDS120CNSummary:BadLOOKernel', ...
      'Step-kernel weights sum to zero for held-out session %s.',A.fileName);
  end
  weights = weights/weightSum;
  C = getStreamNoise(A,'change');
  N = getStreamNoise(A,'noChange');
  part = nan(A.nTrials,6);
  for stream = 1:3
    part(:,stream) = sum(C(A.stepFrames,:,stream).*weights,1,'omitnan')';
    part(:,stream+3) = sum(N(A.stepFrames,:,stream).*weights,1,'omitnan')';
  end
  parts{held} = part;
end
P = vertcat(parts{:});
end

function noise = getStreamNoise(A,side)
switch lower(side)
  case 'change'
    candidates = {'changeStreamNoiseByFrameTrial','changeDirNoiseByFrameTrial'};
  case 'nochange'
    candidates = {'noChangeStreamNoiseByFrameTrial','noChangeDirNoiseByFrameTrial'};
  otherwise
    error('plotIDS120CNSummary:BadSide','Unknown side %s.',side);
end
for k = 1:numel(candidates)
  if isfield(A,candidates{k})
    noise = A.(candidates{k});
    return
  end
end
error('plotIDS120CNSummary:MissingNoise','Analysis lacks %s-side stream noise.',side);
end

function validateTimeBase(A,B)
if numel(A.tMS) ~= numel(B.tMS) || any(A.tMS(:) ~= B.tMS(:)) || ...
    numel(A.stepFrames) ~= numel(B.stepFrames) || ...
    any(A.stepFrames(:) ~= B.stepFrames(:))
  error('plotIDS120CNSummary:TimeVectorMismatch', ...
    'Sessions %s and %s have incompatible time bases.',A.fileName,B.fileName);
end
end

function kernels = computeKernelSet(analyses,nBootstrap)
definitions = struct( ...
  'name',{'cDrift','cProbe','nDrift','nProbe'}, ...
  'title',{'Change: Drift','Change: Probes','NoChange: Drift','NoChange: Probes'}, ...
  'side',{'change','change','noChange','noChange'}, ...
  'streams',{1,[2 3],1,[2 3]});
classKernels = struct();
for k = 1:numel(definitions)
  stats = collectKernelStats(analyses,definitions(k).side,definitions(k).streams);
  classKernels.(definitions(k).name) = packageKernel(stats,definitions(k).title, ...
    nBootstrap,1000+k);
end

individualNames = {'cDrift','cPlus','cMinus','nDrift','nPlus','nMinus'};
individualTitles = {'C Drift','C +Probe','C -Probe','N Drift','N +Probe','N -Probe'};
individualSides = {'change','change','change','noChange','noChange','noChange'};
individualStreams = [1 2 3 1 2 3];
individual = struct();
for k = 1:numel(individualNames)
  stats = collectKernelStats(analyses,individualSides{k},individualStreams(k));
  individual.(individualNames{k}) = packageKernel(stats,individualTitles{k}, ...
    nBootstrap,2000+k);
end
kernels = struct();
kernels.class = classKernels;
kernels.individual = individual;
kernels.classOrder = {definitions.name};
kernels.individualOrder = individualNames;
end

function stats = collectKernelStats(analyses,side,streams)
nSessions = numel(analyses);
tMS = analyses{1}.tMS(:);
nFrames = numel(tMS);
sumCorrect = zeros(nFrames,nSessions);
sumError = zeros(nFrames,nSessions);
nCorrect = zeros(1,nSessions);
nError = zeros(1,nSessions);
for s = 1:nSessions
  A = analyses{s};
  validateTimeBase(analyses{1},A);
  noise = getStreamNoise(A,side);
  x = reshape(noise(:,:,streams),nFrames,[]);
  useTrial = logical(A.trialTable.hasStepNoise(:));
  correct = logical(A.trialTable.correct(:));
  use = repmat(useTrial,numel(streams),1);
  choices = repmat(correct,numel(streams),1);
  sumCorrect(:,s) = sum(x(:,use & choices),2,'omitnan');
  sumError(:,s) = sum(x(:,use & ~choices),2,'omitnan');
  nCorrect(s) = sum(use & choices);
  nError(s) = sum(use & ~choices);
end
stats = struct('tMS',tMS,'stepFrames',analyses{1}.stepFrames(:), ...
  'sumCorrect',sumCorrect,'sumError',sumError, ...
  'nCorrect',nCorrect,'nError',nError);
end

function K = packageKernel(stats,titleText,nBootstrap,seed)
K = kernelFromSessionSample(stats,1:numel(stats.nCorrect));
K.title = titleText;
K.nSessions = numel(stats.nCorrect);
K.preStepIntegralCI95 = [nan nan];
K.stepIntegralCI95 = [nan nan];
if nBootstrap > 0 && K.nSessions > 1
  oldState = rng;
  cleanup = onCleanup(@() rng(oldState)); %#ok<NASGU>
  rng(seed,'twister');
  bootPre = nan(nBootstrap,1);
  bootStep = nan(nBootstrap,1);
  for b = 1:nBootstrap
    sample = randi(K.nSessions,K.nSessions,1);
    B = kernelFromSessionSample(stats,sample);
    bootPre(b) = B.preStepIntegral;
    bootStep(b) = B.stepIntegral;
  end
  K.preStepIntegralCI95 = empiricalCI(bootPre);
  K.stepIntegralCI95 = empiricalCI(bootStep);
end
end

function K = kernelFromSessionSample(stats,sample)
nC = sum(stats.nCorrect(sample));
nE = sum(stats.nError(sample));
if nC == 0 || nE == 0
  error('plotIDS120CNSummary:EmptyChoiceClass', ...
    'A kernel sample lacks correct or error trials.');
end
K = struct();
K.tMS = stats.tMS;
K.stepFrames = stats.stepFrames;
K.meanCorrect = sum(stats.sumCorrect(:,sample),2)/nC;
K.meanError = sum(stats.sumError(:,sample),2)/nE;
K.meanDiff = K.meanCorrect-K.meanError;
K.nCorrect = nC;
K.nError = nE;
K.nTrials = nC+nE;
preFrames = find(K.tMS < K.tMS(K.stepFrames(1)));
K.preFrames = preFrames;
K.preStepMean = mean(K.meanDiff(preFrames),'omitnan');
K.stepMean = mean(K.meanDiff(K.stepFrames),'omitnan');
K.preStepIntegral = sum(K.meanDiff(preFrames),'omitnan');
K.stepIntegral = sum(K.meanDiff(K.stepFrames),'omitnan');
K.meanReference = nan(size(K.meanDiff));
K.meanReference(preFrames) = K.preStepMean;
K.meanReference(K.stepFrames) = K.stepMean;
end

function ci = empiricalCI(values)
values = sort(values(isfinite(values)));
if isempty(values)
  ci = [nan nan];
  return
end
n = numel(values);
indices = max(1,min(n,round([0.025 0.975]*(n-1)+1)));
ci = values(indices)';
end

function T = makeGainSummaryTable(G)
fit = G.sideRole;
[gainCI,gainCIMethod] = preferredGainCI(fit);
contrast = G.contrasts.changeProbeMinusNoChangeProbe;
[contrastCI,contrastP,contrastCIMethod] = preferredContrastInference(contrast);

term = [string(fit.predictorNames(:)); string(contrast.label)];
estimate = [fit.gain(:); contrast.estimate];
CI95Lower = [gainCI(:,1); contrastCI(1)];
CI95Upper = [gainCI(:,2); contrastCI(2)];
PValue = [nan(numel(fit.gain),1); contrastP];
CIMethod = [repmat(string(gainCIMethod),numel(fit.gain),1); ...
  string(contrastCIMethod)];
isContrast = [false(numel(fit.gain),1); true];
nEffectiveCohClipped = repmat(fit.nEffectiveCohClipped,numel(term),1);
nTrials = repmat(fit.nTrials,numel(term),1);
T = table(term,estimate,CI95Lower,CI95Upper,PValue,CIMethod,isContrast, ...
  nEffectiveCohClipped,nTrials);
end

function T = makeInteractionSummaryTable(G)
fit = G.noChangeProbeInteraction;
k = numel(fit.gain);
[ci,ciMethod] = preferredGainCI(fit);
C = fit.comparison;
scaling = fit.interactionScaling;
T = table(fit.gain(k),ci(k,1),ci(k,2),string(ciMethod), ...
  C.deltaNegLogLikelihood,C.likelihoodRatio,C.df,C.pValue, ...
  scaling.plusPredictorSD,scaling.minusPredictorSD, ...
  scaling.centeredPredictorCorrelation,scaling.interactionSD, ...
  fit.nEffectiveCohClipped,fit.nTrials, ...
  'VariableNames',{'interactionGain','CI95Lower','CI95Upper','CIMethod', ...
  'deltaNegLogLikelihood','likelihoodRatio','likelihoodRatioDF', ...
  'likelihoodRatioP','nPlusPredictorSD','nMinusPredictorSD', ...
  'centeredPredictorCorrelation','interactionPredictorSD', ...
  'nEffectiveCohClipped','nTrials'});
end

function [ci,label] = preferredGainCI(fit)
useBootstrap = isfield(fit,'bootstrapNSuccess') && fit.bootstrapNSuccess > 0 && ...
  isfield(fit,'bootstrapCI95') && all(isfinite(fit.bootstrapCI95),'all');
if useBootstrap
  ci = fit.bootstrapCI95;
  label = sprintf('Session bootstrap (%d/%d successful)', ...
    fit.bootstrapNSuccess,fit.bootstrapNRequested);
else
  ci = fit.CI95;
  label = fit.covarianceMethod;
end
end

function [ci,pValue,label] = preferredContrastInference(contrast)
if all(isfinite(contrast.bootstrapCI95)) && isfinite(contrast.bootstrapPValue)
  ci = contrast.bootstrapCI95;
  pValue = contrast.bootstrapPValue;
  label = 'Session bootstrap';
else
  ci = contrast.CI95;
  pValue = contrast.pValue;
  label = 'Covariance contrast';
end
end

function R = makeSessionRecords(analyses,fits)
n = numel(analyses);
fileName = strings(n,1);
sourceSessionID = strings(n,1);
sourceDomain = strings(n,1);
randomBaseDir = nan(n,1);
nTrials = nan(n,1);
nStepNoiseTrials = nan(n,1);
threshold75 = nan(n,1);
for k = 1:n
  A = analyses{k};
  fileName(k) = A.fileName;
  sourceSessionID(k) = A.sourceSessionID;
  sourceDomain(k) = string(A.sessionProbeHeader.sourceDomain);
  if isfield(A.sessionProbeHeader,'randomBaseDir')
    randomBaseDir(k) = double(A.sessionProbeHeader.randomBaseDir);
  end
  nTrials(k) = A.nTrials;
  nStepNoiseTrials(k) = A.nStepNoiseTrials;
  threshold75(k) = fits.threshold(fits.sessionIndex == k);
end
R = table(fileName,sourceSessionID,sourceDomain,randomBaseDir,nTrials, ...
  nStepNoiseTrials,threshold75);
end

function fig = plotSummary(S,fig)
if nargin < 2
  fig = figure('Color','w');
end
tl = tiledlayout(fig,4,3,'TileSpacing','loose','Padding','compact');
title(tl,sprintf('120° C/N Summary: %s, %s (%d sessions)', ...
  S.datasetLabel,S.animal,S.nSessions),'FontWeight','bold');
plotPsychometric(nexttile(tl,1),S.trialTable,S.psychFit);
plotDataText(nexttile(tl,2),S);
plotModelText(nexttile(tl,3),S);

colors = struct('cDrift',[0 .45 .8],'cProbe',[.25 .7 .85], ...
  'nDrift',[.85 .35 .2],'nProbe',[.9 .65 .25]);
names = S.kernels.classOrder;
kernelAxes = gobjects(4,1);
patches = gobjects(4,1);
for k = 1:4
  kernelAxes(k) = nexttile(tl,3+k);
  K = S.kernels.class.(names{k});
  patches(k) = plotKernel(kernelAxes(k),K,colors.(names{k}));
end
limits = cell2mat(get(kernelAxes,'YLim'));
shared = [min(limits(:,1)) max(limits(:,2))];
if diff(shared) == 0
  shared = shared + [-1 1];
end
set(kernelAxes,'YLim',shared);
for k = 1:4
  patches(k).YData = shared([1 1 2 2]);
end

gainLim = gainLimits(S.noiseGain);
plotGain(nexttile(tl,8),S.noiseGain.sideRole,'C/N Drift-Probe Gains',gainLim,[.45 .45 .75]);
plotGain(nexttile(tl,9),S.noiseGain.full,'Six Unconstrained Gains',gainLim,[.6 .6 .6]);
plotGain(nexttile(tl,10),S.noiseGain.cnp,'C, N, Probe Gains',gainLim,[.45 .65 .45]);
plotGain(nexttile(tl,11),S.noiseGain.side,'RDK-Side Gains',gainLim,[.65 .45 .65]);
plotKernelIntegrals(nexttile(tl,12),S.kernels.class,names,colors);
end

function plotPsychometric(ax,T,psych)
hold(ax,'on');
fit = psych.alignedWeibull;
x = T.alignedCoh;
y = T.correct;
edges = linspace(min(x),max(x),10);
centers = (edges(1:end-1)+edges(2:end))/2;
for k = 1:numel(centers)
  use = x >= edges(k) & x < edges(k+1);
  if k == numel(centers)
    use = x >= edges(k) & x <= edges(k+1);
  end
  if any(use)
    plot(ax,centers(k),mean(y(use)),'ko','MarkerFaceColor','k','MarkerSize',4);
  end
end
gridX = linspace(0,max(x)*1.05,300);
plot(ax,gridX,idqWeibullP(gridX,fit.alpha,fit.betaWeibull,fit.lapse), ...
  '-','LineWidth',1.5);
xline(ax,fit.threshold,':');
yline(ax,fit.thresholdPerformance,':');
xlabel(ax,'Norm. Coh. Step'); ylabel(ax,'P(hit)');
title(ax,'All Directions');
ylim(ax,[.48 1.02]); xlim(ax,[0 max(gridX)]); grid(ax,'on'); box(ax,'off');
normCoh = median(T.stepCoh./T.alignedCoh,'omitnan');
txt = sprintf('%d trials\n%.1f%% correct\nthresh. = %.1f%%\n\\beta=%.2f\n\\lambda=%.3f', ...
  height(T),100*mean(T.correct),fit.threshold*normCoh,fit.betaWeibull,fit.lapse);
text(ax,.98,.02,txt,'Units','normalized','HorizontalAlignment','right', ...
  'VerticalAlignment','bottom','FontSize',8);
end

function plotDataText(ax,S)
axis(ax,'off');
noiseText = sprintf('%.4g, ',S.noiseAmplitudesPC);
if isempty(noiseText)
  noiseText = 'unknown';
else
  noiseText = noiseText(1:end-2);
end
txt = {sprintf('Created: %s',string(S.createdAt,'dd/MM/yyyy HH:mm')); ...
  sprintf('Dataset: %s',S.datasetLabel); sprintf('Animal: %s',S.animal); ...
  sprintf('Sessions: %d',S.nSessions); sprintf('Trials: %d',S.nTrials); ...
  sprintf('Noise trials: %d',S.nStepNoiseTrials); ...
  sprintf('Noise amplitude(s): %s%%',noiseText); ...
  'Gain window: step only'; 'Kernel display: prestep + step'; ...
  sprintf('Kernel bootstrap: %d session resamples',S.nBootstrap); ...
  sprintf('Gain bootstrap: %d/%d successful', ...
    S.noiseGain.sideRole.bootstrapNSuccess, ...
    S.noiseGain.sideRole.bootstrapNRequested)};
text(ax,0,1,txt,'Units','normalized','VerticalAlignment','top', ...
  'FontName','Menlo','FontSize',7.5);
end

function plotModelText(ax,S)
axis(ax,'off');
Q = S.noiseGain.symmetry;
[ci,ciLabel] = preferredGainCI(S.noiseGain.sideRole);
g = S.noiseGain.sideRole.gain;
C = S.noiseGain.contrasts.changeProbeMinusNoChangeProbe;
[contrastCI,contrastP] = preferredContrastInference(C);
nClipped = S.noiseGain.sideRole.nEffectiveCohClipped;
nGainTrials = S.noiseGain.sideRole.nTrials;
interactionFit = S.noiseGain.noChangeProbeInteraction;
[interactionCI,interactionCILabel] = preferredGainCI(interactionFit);
interactionIndex = numel(interactionFit.gain);
interaction = interactionFit.gain(interactionIndex);
interactionComparison = interactionFit.comparison;
predictorCorrelation = interactionFit.interactionScaling.centeredPredictorCorrelation;
txt = {'Primary model:'; ...
  'c_{eff}=c_{step}+g_{CD}n_{CD}+g_{CP}(n_{C+}+n_{C-})'; ...
  '       +g_{ND}n_{ND}+g_{NP}(n_{N+}+n_{N-})'; ...
  sprintf('Gain CI: %s',ciLabel); ...
  sprintf('CD %.3f [%.3f, %.3f]; CP %.3f [%.3f, %.3f]', ...
    g(1),ci(1,:),g(2),ci(2,:)); ...
  sprintf('ND %.3f [%.3f, %.3f]; NP %.3f [%.3f, %.3f]', ...
    g(3),ci(3,:),g(4),ci(4,:)); ...
  sprintf('CP-NP %.3f [%.3f, %.3f], p=%.3g', ...
    C.estimate,contrastCI,contrastP); ...
  sprintf('N probe interaction %.3f [%.3f, %.3f]', ...
    interaction,interactionCI(interactionIndex,:)); ...
  sprintf('Interaction CI: %s',interactionCILabel); ...
  sprintf('Interaction dNLL %.2f, LRT p=%.3g; predictor r=%.3f', ...
    interactionComparison.deltaNegLogLikelihood, ...
    interactionComparison.pValue,predictorCorrelation); ...
  sprintf('Effective coherence clipped: %d/%d (%.2f%%)', ...
    nClipped,nGainTrials,100*nClipped/nGainTrials); ...
  sprintf('Symmetry p: C probes %.3g; N probes %.3g; N all %.3g', ...
    Q.changeProbe.pValue,Q.noChangeProbe.pValue,Q.noChangeAll.pValue)};
text(ax,0,1,txt,'Units','normalized','VerticalAlignment','top', ...
  'FontName','Menlo','FontSize',5.7);
end

function patchHandle = plotKernel(ax,K,color)
hold(ax,'on');
plot(ax,K.tMS,K.meanDiff,'Color',color,'LineWidth',1.3);
plot(ax,K.tMS,K.meanReference,'k-','LineWidth',1);
yline(ax,0,'k:');
yl = ylim(ax);
x1 = K.tMS(K.stepFrames(1));
if numel(K.tMS) > K.stepFrames(end)
  x2 = K.tMS(K.stepFrames(end)+1);
else
  x2 = K.tMS(K.stepFrames(end)) + median(diff(K.tMS));
end
patchHandle = patch(ax,[x1 x2 x2 x1],yl([1 1 2 2]),[.7 .7 .7], ...
  'EdgeColor','none','FaceAlpha',.25,'HandleVisibility','off');
plot(ax,K.tMS,K.meanDiff,'Color',color,'LineWidth',1.3);
plot(ax,K.tMS,K.meanReference,'k-','LineWidth',1);
title(ax,K.title); xlabel(ax,'Trial Time (ms)'); ylabel(ax,'Kernel (% coh.)');
grid(ax,'on'); box(ax,'off');
txt = sprintf(['Pre: mean %.3f, \\Sigma %.2f [%.2f, %.2f]\n' ...
  'Step: mean %.3f, \\Sigma %.2f [%.2f, %.2f]'], ...
  K.preStepMean,K.preStepIntegral,K.preStepIntegralCI95, ...
  K.stepMean,K.stepIntegral,K.stepIntegralCI95);
text(ax,.02,.98,txt,'Units','normalized','VerticalAlignment','top','FontSize',7.5);
end

function plotGain(ax,fit,titleText,limits,color)
g = fit.gain(:);
[ci,~] = preferredGainCI(fit);
x = (1:numel(g))';
bar(ax,x,g,.55,'FaceColor',color,'EdgeColor','k'); hold(ax,'on');
for k = 1:numel(g)
  plot(ax,[x(k) x(k)],ci(k,:),'k-','LineWidth',1.3);
end
plot(ax,x,g,'ko','MarkerFaceColor','k');
yline(ax,0,'k:');
set(ax,'XTick',x,'XTickLabel',fit.predictorNames);
xtickangle(ax,25);
xlim(ax,[.25 numel(g)+.75]); ylim(ax,limits);
ylabel(ax,'Noise gain'); title(ax,titleText); grid(ax,'on'); box(ax,'off');
end

function plotKernelIntegrals(ax,kernels,names,colors)
values = nan(4,1); ci = nan(4,2); labels = cell(4,1); C = nan(4,3);
for k = 1:4
  K = kernels.(names{k});
  values(k) = K.stepIntegral;
  ci(k,:) = K.stepIntegralCI95;
  labels{k} = K.title;
  C(k,:) = colors.(names{k});
end
b = bar(ax,1:4,values,.65,'FaceColor','flat','EdgeColor','k');
b.CData = C; hold(ax,'on');
for k = 1:4
  plot(ax,[k k],ci(k,:),'k-','LineWidth',1.3);
end
yline(ax,0,'k:');
set(ax,'XTick',1:4,'XTickLabel',labels); xtickangle(ax,25);
ylabel(ax,'Step \Sigma kernel (% coh.-frames)');
title(ax,'Step Kernel Integrals'); grid(ax,'on'); box(ax,'off');
end

function limits = gainLimits(G)
fits = {G.full,G.sideRole,G.cnp,G.side};
values = [];
for k = 1:numel(fits)
  [ci,~] = preferredGainCI(fits{k});
  values = [values; fits{k}.gain(:); ci(:)]; %#ok<AGROW>
end
values = values(isfinite(values));
if isempty(values)
  limits = [-1 1];
  return
end
span = max(1,max(values)-min(values));
limits = [min(0,min(values)-.1*span) max(0,max(values)+.1*span)];
if diff(limits) == 0
  limits = limits + [-1 1];
end
end
