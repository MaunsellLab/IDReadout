function summaries = plotIDSSummary(varargin)
% plotIDSSummary  Create one IDQ-style across-session summary per IDS offset.

p = inputParser;
p.addParameter('Animal', 'All', @(x) ischar(x) || isstring(x));
p.addParameter('IncludeIDQ', true, @(x) islogical(x) && isscalar(x));
p.addParameter('NoisePredictorType', 'looStep', @(x) ischar(x) || isstring(x));
p.addParameter('ProbeOffsets', [], @isnumeric);
p.parse(varargin{:});
opts = p.Results;

root = domainFolder(mfilename('fullpath'));
dirs = dir(fullfile(root, 'Data', 'Probe*'));
dirs = dirs([dirs.isdir]);
summaries = cell(0,1);
for k = 1:numel(dirs)
  token = regexp(dirs(k).name, '^Probe(\d+)$', 'tokens', 'once');
  if isempty(token), continue; end
  offset = str2double(token{1});
  if ~isempty(opts.ProbeOffsets) && ~ismember(offset, opts.ProbeOffsets)
    continue; 
  end
  analysisFolder = fullfile(dirs(k).folder, dirs(k).name, 'SessionAnalysis');
  if ~isfolder(analysisFolder)
    continue;
  end
  dataFilePaths = selectAnalysisFiles(analysisFolder, ...
    'FilePattern', '*_sessionAnalysis.mat', 'Animal', opts.Animal, ...
    'IncludeIDQ', opts.IncludeIDQ);
  if isempty(dataFilePaths), continue; end
  analyses = cell(numel(dataFilePaths),1);
  for j = 1:numel(dataFilePaths)
    S = load(dataFilePaths{j}, 'sessionAnalysis');
    analyses{j} = S.sessionAnalysis;
  end
  summaries{end+1,1} = summarizeOffset(root, dirs(k).name, offset, analyses, ...
    char(opts.NoisePredictorType), char(opts.Animal)); %#ok<AGROW>
end
if isscalar(summaries), summaries = summaries{1}; end
end

function S = summarizeOffset(root, probeTag, offset, analyses, predictorType, animal)
for k = 1:numel(analyses)
  assert(analyses{k}.sessionProbeHeader.probeOffsetDeg == offset, ...
    'plotIDSSummary:OffsetMismatch', 'Mixed offsets in %s.', probeTag);
end
T = concatenateTables(analyses);
[psych, T] = fitPsychometric(T);
T = selectPredictors(T, analyses, predictorType);
gain = fitIDSNoiseGain(T, psych.sessionFits, psych.alignedWeibull, 0.75);

S = struct();
S.probeTag = probeTag;
S.probeOffsetDeg = offset;
S.createdBy = mfilename;
S.createdAt = datetime('now');
S.animal = animal;
S.noisePredictorType = predictorType;
S.nSessions = numel(analyses);
S.nTrials = height(T);
S.nStepNoiseTrials = sum(T.hasStepNoise);
S.noiseAmplitudesPC = unique(round(T.cohNoisePC(isfinite(T.cohNoisePC)), 6))';
S.trialTable = T;
S.psychFit = psych;
S.absKernels = computeAbsoluteKernels(analyses);
S.relKernels = computeRelativeKernels(analyses, offset);
S.noiseGain = gain;
S.sessionRecords = makeSessionRecords(analyses, psych.sessionFits);

summaryFolder = validFolder(fullfile(root, 'Data', probeTag, 'AcrossSessionSummaries'));
plotFolder = validFolder(fullfile(root, 'Plots', probeTag));
acrossSummary = S; 
save(fullfile(summaryFolder, sprintf('IDS_%s_AcrossSessionSummary.mat', probeTag)), ...
  'acrossSummary', '-v7.3');
% Create the usual docked figure for interactive viewing.
plotSummary(S);
% Independently render the same summary into a temporary portrait figure.
exportFig = figure('Color', 'w', 'WindowStyle', 'normal', 'Visible', 'off', 'Units', 'inches', 'Position', [1 1 8 10.5]);
% Ensure that the temporary figure is deleted even if export fails.
cleanup = onCleanup(@() delete(exportFig));
plotSummary(S, exportFig);
drawnow;
pdfPath = fullfile(plotFolder, sprintf('IDS%sSessionSummary.pdf', probeTag));
exportgraphics(exportFig, pdfPath, 'ContentType', 'vector', 'Units', 'inches', ...
  'Width', 8, 'Height', 10.5, 'Padding', 0.1, 'PreserveAspectRatio', 'on');
end

function T = concatenateTables(analyses)
parts = cell(numel(analyses),1);
for k = 1:numel(analyses)
  q = analyses{k}.trialTable;
  q.sessionIndex = repmat(k, height(q), 1);
  q.noisyStepCoh = repmat(analyses{k}.noisyStepCoh, height(q), 1);
  parts{k} = q;
end
T = vertcat(parts{:});
T = movevars(T, 'sessionIndex', 'Before', 1);
end

function [psych, T] = fitPsychometric(T)
target = 0.75;
initial = fitIDQInitialSessionThresholds(T, 2, 0.02, target);
aligned = fitIDQAcrossAlignedWeibull(initial.alignedCoh, T.correct, target, [0 0.05]);
T.alignedCoh = initial.alignedCoh;
T.noisyStepAlignedCoh = initial.noisyStepAlignedCoh;
psych = struct('title', 'All Directions', 'targetPerformance', target, ...
  'sessionFits', initial.sessionFits, 'alignedWeibull', aligned);
end

function T = selectPredictors(T, analyses, predictorType)
switch predictorType
  case {'rectStep','rectFull'}
    stem = predictorType;
    T.noisePredDrift = T.([stem 'ChangeNoisePredDrift']);
    T.noisePredPlus = T.([stem 'ChangeNoisePredPlus']);
    T.noisePredMinus = T.([stem 'ChangeNoisePredMinus']);
  case {'looStep','looFull'}
    P = computeLOOPredictors(analyses, lower(erase(predictorType, 'loo')));
    T.noisePredDrift = P(:,1); T.noisePredPlus = P(:,2); T.noisePredMinus = P(:,3);
  otherwise
    error('plotIDSSummary:BadPredictorType', 'Unknown predictor type %s.', predictorType);
end
end

%%-------------------------------------------------
function P = computeLOOPredictors(analyses, window)
nSessions = numel(analyses);
parts = cell(nSessions,1);
for held = 1:nSessions
  A = analyses{held};
  useFrames = chooseFrames(A, window);
  correctNoise = []; errorNoise = [];
  for source = 1:nSessions
    if nSessions > 1 && source == held, continue; end
    B = analyses{source};
    validateTimeBase(A, B);
    use = B.trialTable.hasStepNoise;
    correctNoise = [correctNoise B.changeSumNoiseByFrameTrial(:,use & B.trialTable.correct)]; %#ok<AGROW>
    errorNoise = [errorNoise B.changeSumNoiseByFrameTrial(:,use & ~B.trialTable.correct)]; %#ok<AGROW>
  end
  kernel = mean(correctNoise,2,'omitnan') - mean(errorNoise,2,'omitnan');
  weights = kernel(useFrames);
  weightSum = sum(weights, 'omitnan');
  if ~isfinite(weightSum) || weightSum == 0
    error('plotIDSSummary:BadLOOKernel', 'Kernel weights sum to zero for %s.', A.fileName);
  end
  weights = weights / weightSum;
  part = nan(A.nTrials,3);
  for stream = 1:3
    x = A.changeStreamNoiseByFrameTrial(useFrames,:,stream);
    part(:,stream) = sum(x .* weights, 1, 'omitnan')';
  end
  parts{held} = part;
end
P = vertcat(parts{:});
end

function frames = chooseFrames(A, window)
switch window
  case 'step', frames = A.stepFrames;
  case 'full', frames = (1:numel(A.tMS))';
  otherwise, error('plotIDSSummary:BadWindow', 'Unknown temporal window %s.', window);
end
end

function validateTimeBase(A, B)
if numel(A.tMS) ~= numel(B.tMS) || any(A.tMS(:) ~= B.tMS(:)) || ...
    numel(A.stepFrames) ~= numel(B.stepFrames) || any(A.stepFrames(:) ~= B.stepFrames(:))
  error('plotIDSSummary:TimeVectorMismatch', ...
    'Sessions %s and %s have incompatible time bases.', A.fileName, B.fileName);
end
end

function kernels = computeRelativeKernels(analyses, offset)
titles = {'Drift Direction', sprintf('Drift + %d°', offset), sprintf('Drift - %d°', offset)};
kernels = cell(1,3);
for stream = 1:3
  [correctNoise, errorNoise, tMS, stepFrames] = collectRelative(analyses, stream);
  kernels{stream} = packageKernel(tMS, stepFrames, correctNoise, errorNoise, titles{stream});
end
end

function [C,E,tMS,stepFrames] = collectRelative(analyses, stream)
tMS = analyses{1}.tMS; stepFrames = analyses{1}.stepFrames; C = []; E = [];
for k = 1:numel(analyses)
  A = analyses{k}; validateTimeBase(analyses{1}, A);
  use = A.trialTable.hasStepNoise;
  x = A.changeStreamNoiseByFrameTrial(:,:,stream);
  C = [C x(:,use & A.trialTable.correct)]; %#ok<AGROW>
  E = [E x(:,use & ~A.trialTable.correct)]; %#ok<AGROW>
end
end

function kernels = computeAbsoluteKernels(analyses)
edges = [0 120 240 360];
titles = {'0°-119°','120°-239°','240°-359°'};
kernels = cell(1,3);
for sector = 1:3
  C = []; E = []; tMS = analyses{1}.tMS; stepFrames = analyses{1}.stepFrames;
  for k = 1:numel(analyses)
    A = analyses{k}; validateTimeBase(analyses{1}, A);
    T = A.trialTable;
    directions = [T.driftDirDeg T.plusProbeDirDeg T.minusProbeDirDeg];
    for stream = 1:3
      inSector = directions(:,stream) >= edges(sector) & directions(:,stream) < edges(sector+1);
      use = T.hasStepNoise & inSector;
      x = A.changeStreamNoiseByFrameTrial(:,:,stream);
      C = [C x(:,use & T.correct)]; %#ok<AGROW>
      E = [E x(:,use & ~T.correct)]; %#ok<AGROW>
    end
  end
  kernels{sector} = packageKernel(tMS, stepFrames, C, E, titles{sector});
end
end

function K = packageKernel(tMS, stepFrames, C, E, titleText)
K = struct(); K.title = titleText; K.tMS = tMS; K.stepFrames = stepFrames;
K.meanCorrect = mean(C,2,'omitnan'); K.meanError = mean(E,2,'omitnan');
K.meanDiff = K.meanCorrect - K.meanError;
K.nCorrect = size(C,2); K.nError = size(E,2); K.nTrials = K.nCorrect + K.nError;
preFrames = find(tMS < tMS(stepFrames(1)));
K.preStepMean = mean(K.meanDiff(preFrames),'omitnan');
K.stepMean = mean(K.meanDiff(stepFrames),'omitnan');
K.rectReference = repmat(K.preStepMean, size(K.meanDiff));
K.rectReference(stepFrames) = K.stepMean;
end

function R = makeSessionRecords(analyses, fits)
n = numel(analyses); fileName = strings(n,1); sourceSessionID = strings(n,1);
sourceDomain = strings(n,1);
nTrials = nan(n,1); nStepNoiseTrials = nan(n,1); threshold75 = nan(n,1);
for k = 1:n
  fileName(k) = analyses{k}.fileName; sourceSessionID(k) = analyses{k}.sourceSessionID;
  if isfield(analyses{k}.sessionProbeHeader, 'sourceDomain')
    sourceDomain(k) = string(analyses{k}.sessionProbeHeader.sourceDomain);
  end
  nTrials(k) = analyses{k}.nTrials; nStepNoiseTrials(k) = analyses{k}.nStepNoiseTrials;
  threshold75(k) = fits.threshold(fits.sessionIndex == k);
end
R = table(fileName, sourceSessionID, sourceDomain, nTrials, nStepNoiseTrials, threshold75);
end

%%====================================
function fig = plotSummary(S, fig)
if nargin < 2
  fig = figure('Color','w');   % Uses your docked-figure default
end
tl = tiledlayout(fig, 4, 3, 'TileSpacing','loose', 'Padding','compact');
title(tl, sprintf('IDS Summary: ±%d° Probe (%d sessions)', ...
  S.probeOffsetDeg, S.nSessions), 'FontWeight','bold');
plotPsychometric(nexttile(tl,1), S.trialTable, S.psychFit);
plotDataText(nexttile(tl,2), S);
plotModelText(nexttile(tl,3), S);
absColor = [0 .5 .8]; relColor = [.8 .5 1]; allColor = [.5 .5 .5];
kernelAxes = gobjects(6,1); patches = gobjects(6,1);
for k = 1:3
  kernelAxes(k) = nexttile(tl,3+k);
  patches(k) = plotKernel(kernelAxes(k),S.absKernels{k},absColor, false);
  kernelAxes(k+3) = nexttile(tl,6+k);
  patches(k+3) = plotKernel(kernelAxes(k+3),S.relKernels{k},relColor, k == 3);
end
limits = cell2mat(get(kernelAxes,'YLim')); shared = [min(limits(:,1)) max(limits(:,2))];
set(kernelAxes,'YLim',shared);
for k = 1:6, patches(k).YData = shared([1 1 2 2]); end
gainLim = gainLimits(S.noiseGain);
plotGain(nexttile(tl,10),S.noiseGain.driftRelative,'Drift-Rel. Gains',gainLim,relColor);
plotGain(nexttile(tl,11),S.noiseGain.combined,'Combined Gain',gainLim,allColor);
plotGain(nexttile(tl,12),S.noiseGain.driftProbe,'Drift v. Probe Gains',gainLim,relColor);
end

function plotPsychometric(ax,T,psych)
hold(ax,'on'); fit = psych.alignedWeibull; x = T.alignedCoh; y = T.correct;
edges = linspace(min(x),max(x),10); centers = (edges(1:end-1)+edges(2:end))/2;
for k = 1:numel(centers)
  use = x>=edges(k) & x<(edges(k+1)); if k==numel(centers), use=x>=edges(k)&x<=edges(k+1); end
  if any(use)
    plot(ax,centers(k),mean(y(use)),'ko','MarkerFaceColor','k','MarkerSize',4);
    text(ax,centers(k),mean(y(use)),sprintf(' %d',sum(use)),'FontSize',7,'VerticalAlignment','bottom');
  end
end
gridX = linspace(0,max(x)*1.05,300);
plot(ax,gridX,idqWeibullP(gridX,fit.alpha,fit.betaWeibull,fit.lapse),'-','LineWidth',1.5);
xline(ax,fit.threshold,':'); yline(ax,fit.thresholdPerformance,':');
xlabel(ax,'Norm. Coh. Step'); ylabel(ax,'P(hit)'); title(ax,'All Directions');
ylim(ax,[.48 1.02]); xlim(ax,[0 max(gridX)]); grid(ax,'on'); box(ax,'off');
normCoh = median(T.stepCoh./T.alignedCoh,'omitnan');
txt = sprintf('%d trials\n%.1f%% correct\n1 = %.1f%% coh.\nthresh. = %.1f%%\n\\beta=%.2f\n\\lambda=%.3f', ...
  height(T),100*mean(T.correct),normCoh,fit.threshold*normCoh,fit.betaWeibull,fit.lapse);
text(ax,.98,.02,txt,'Units','normalized','HorizontalAlignment','right', ...
  'VerticalAlignment','bottom','FontSize',8);
end

function plotDataText(ax,S)
axis(ax,'off');
noiseText = sprintf('%.4g, ', S.noiseAmplitudesPC);
if isempty(noiseText)
  noiseText = 'unknown';
else
  noiseText = noiseText(1:end - 2);
end
txt = {sprintf('Created: %s',string(S.createdAt,'dd/MM/yyyy HH:mm')); ...
  sprintf('Offset: ±%d°',S.probeOffsetDeg); sprintf('Sessions: %d',S.nSessions); ...
  sprintf('Trials: %d',S.nTrials); sprintf('Noise trials: %d',S.nStepNoiseTrials); ...
  sprintf('Noise amplitude(s): %s%%', noiseText); 'Streams: independent'; ...
  sprintf('Predictor: %s',S.noisePredictorType); ...
  'Probe0 = drift + offset'; 'Probe1 = drift - offset'};
text(ax,0,1,txt,'Units','normalized','VerticalAlignment','top','FontName','Menlo','FontSize',7.5);
end

function plotModelText(ax,S)
axis(ax,'off'); N=S.noiseGain.normalized;
txt = {'Drift-relative:'; ' c_eff = c_step + g_D n_D + g_+ n_+ + g_- n_-'; ''; ...
  'Combined:'; ' c_eff = c_step + g_all(n_D + n_+ + n_-)'; ''; ...
  'Drift v. Probe:'; ' c_eff = c_step + g_D n_D + g_P(n_+ + n_-)'; ''; ...
  sprintf('g_P/g_D = %.3f  [%.3f, %.3f]',N.pooledProbeOverDrift,N.pooledCI95); ...
  sprintf('g_+/g_D = %.3f',N.plusProbeOverDrift); ...
  sprintf('g_-/g_D = %.3f',N.minusProbeOverDrift)};
text(ax,0,1,txt,'Units','normalized','VerticalAlignment','top','FontName','Menlo','FontSize',7);
end

function patchHandle = plotKernel(ax, K, color, legendOn)
hold(ax,'on');
plot(ax,K.tMS,K.meanDiff,'Color',color,'LineWidth', 1.2);
plot(ax,K.tMS,K.rectReference,'k-','LineWidth',1);
yl=ylim(ax);
x1=K.tMS(K.stepFrames(1));
x2=K.tMS(K.stepFrames(end));
patchHandle = patch(ax,[x1 x2 x2 x1],yl([1 1 2 2]),[.7 .7 .7], ...
  'EdgeColor','none','FaceAlpha',.35,'HandleVisibility','off');
plot(ax,K.tMS,K.meanDiff,'Color',color,'LineWidth',1.2,'DisplayName','Kernel');
plot(ax,K.tMS,K.rectReference,'k-','LineWidth',1,'DisplayName','Mean Pre/Post Step');
title(ax,K.title); xlabel(ax,'Trial Time (ms)'); ylabel(ax,'Change Side Kernel'); grid(ax,'on'); box(ax,'off');
text(ax,.02,.98,sprintf('PreStep Mean %.2f%%\nStep Mean %.2f%%',K.preStepMean,K.stepMean), ...
  'Units','normalized','VerticalAlignment','top','FontSize',8);
if legendOn
  legend(ax, {'', '', 'Kernel', 'Mean Pre/Post Step'}, 'Location', 'southwest', 'FontSize', 6);
end
end

function plotGain(ax, fit, titleText, limits, color)

g = fit.gain(:);
x = (1:numel(g))';
bar(ax,x,g,.55,'FaceColor',color,'EdgeColor','k'); hold(ax,'on');
for k=1:numel(g)
  plot(ax,[x(k) x(k)],fit.CI95(k,:), 'k-', 'LineWidth', 1.3);
end
plot(ax,x,g,'ko','MarkerFaceColor','k'); 
yline(ax,0,'k:');
set(ax,'XTick',x,'XTickLabel',fit.predictorNames);
xlim(ax,[.25 3.75]);
ylim(ax,limits);
ylabel(ax,'Noise gain');
title(ax,titleText);
grid(ax,'on');
box(ax,'off');
txt = '';
for d = 1:numel(g)
  txt = [txt, sprintf('%s: %.2f\n', fit.predictorNames{d}, g(d))]; %#ok<AGROW>
end
txt = [txt, '(95% CI)']; 
text(ax, 0.98, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
  'FontSize', 8);
end

function limits = gainLimits(G)
fits={G.driftRelative,G.combined,G.driftProbe}; values=[];
for k=1:numel(fits), values=[values;fits{k}.gain(:);fits{k}.CI95(:)]; end %#ok<AGROW>
values=values(isfinite(values)); if isempty(values), limits=[-1 1]; return; end
span=max(1,max(values)-min(values)); limits=[min(0,min(values)-.1*span) max(values)+.1*span];
end
