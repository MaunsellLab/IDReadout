function out = analyzeLocalDetectorEvidenceAverage(baseFolder, alpha, winMS)
% analyzeLocalDetectorEvidenceAverage
%
% Pool first-order local-evidence analysis across all valid sessions,
% using the same scan/exclusion logic as kernelAverage.
%
% For each valid NoiseMatrices session:
%   - compute per-trial Echange, EnoChange, isCorrect
% Then:
%   - concatenate trials across sessions
%   - generate pooled plots
%   - retain session identity for later diagnostics if needed
%
% INPUTS
%   baseFolder : root folder (default folderPath())
%   alpha      : probe weight relative to pref (default 0.19)
%   winMS      : [start stop] relative to step onset in ms (default [0 125])
%
% OUTPUT
%   out : struct containing pooled trialwise data and summaries

if nargin < 1 || isempty(baseFolder)
    baseFolder = folderPath();
end
if nargin < 2 || isempty(alpha)
    alpha = 0.19;
end
if nargin < 3 || isempty(winMS)
    winMS = [0 125];
end

dataFolder = baseFolder + "/Data/NoiseMatrices/";
if ~exist(dataFolder, 'dir')
    error('analyzeLocalDetectorEvidenceAverage:MissingFolder', ...
        'Data folder not found: %s', dataFolder);
end

matFiles = dir(fullfile(dataFolder, '*.mat'));

% Storage
allEchange   = [];
allEnoChange = [];
allCorrect   = [];
allSessionID = [];
allFileNames = {};

sessionOut = {};
sessionHeaders = {};

nSessions = 0;

for f = 1:length(matFiles)
    fileName = matFiles(f).name;
    if endsWith(fileName, '_fileInfo.mat')
        continue;
    end

    S = load(fullfile(dataFolder, fileName));
    header = S.header;

    if excludeFile(header)
        continue
    end

    % Same pref-noise restriction as kernelAverage
    if header.prefNoiseCohPC.data ~= 10
        fprintf("Skipping   %s (%d of %d) -- prefNoiseCohPC is %.0f\n", ...
            fileName, f, length(matFiles), header.prefNoiseCohPC.data);
        continue;
    end

    fprintf("Processing %s (%d of %d)\n", fileName, f, length(matFiles));

    try
        tmp = analyzeLocalDetectorEvidence(fileName, ...
            'alpha', alpha, ...
            'winMS', winMS, ...
            'doPlot', false);
    catch ME
        fprintf(2, 'FAILED on %s: %s\n', fileName, ME.message);
        continue;
    end

    nSessions = nSessions + 1;

    sessionOut{nSessions} = tmp;
    sessionHeaders{nSessions} = header;

    nTr = numel(tmp.Echange);

    allEchange   = [allEchange;   tmp.Echange(:)];
    allEnoChange = [allEnoChange; tmp.EnoChange(:)];
    allCorrect   = [allCorrect;   tmp.isCorrect(:)];
    allSessionID = [allSessionID; nSessions * ones(nTr,1)];
    allFileNames{nSessions,1} = fileName;
end

if nSessions == 0
    error('analyzeLocalDetectorEvidenceAverage:NoSessions', ...
        'No valid sessions found.');
end

% ---- Pooled summaries ----
nBinsChange   = 7;
nBinsNoChange = 7;
nGroupsNoChange = 3;
midFrac = 0.5;

[binCtrChange, pCorrByChange, nByChange, edgesChange, binIdxChange] = ...
    quantileBinnedMean(allEchange, allCorrect, nBinsChange);

[groupIdxNoChange, noChangeCuts] = quantileGroups(allEnoChange, nGroupsNoChange);

pCorrSplit = nan(nBinsChange, nGroupsNoChange);
nSplit     = nan(nBinsChange, nGroupsNoChange);

for g = 1:nGroupsNoChange
    for b = 1:nBinsChange
        idx = (groupIdxNoChange == g) & (binIdxChange == b);
        nSplit(b,g) = sum(idx);
        if nSplit(b,g) > 0
            pCorrSplit(b,g) = mean(allCorrect(idx));
        end
    end
end

lo = quantile(allEchange, (1-midFrac)/2);
hi = quantile(allEchange, 1-(1-midFrac)/2);
usePlot3 = allEchange >= lo & allEchange <= hi;

[binCtrNoChange, pErrByNoChange, nByNoChange, edgesNoChange] = ...
    quantileBinnedMean(allEnoChange(usePlot3), ~allCorrect(usePlot3), nBinsNoChange);

% ---- Plots ----
figBase = 60;

figure(figBase); clf;
plot(binCtrChange, pCorrByChange, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
xlabel('Changed-side evidence, Echange');
ylabel('P(correct)');
title(sprintf('Pooled accuracy vs changed-side evidence (%d sessions)', nSessions));
axis tight; box off;

figure(figBase+1); clf; hold on;
mk = {'o-','s-','d-'};
for g = 1:nGroupsNoChange
    style = mk{min(g, numel(mk))};
    plot(binCtrChange, pCorrSplit(:,g), style, 'LineWidth', 1.5, 'MarkerSize', 7);
end
xlabel('Changed-side evidence, Echange');
ylabel('P(correct)');
title(sprintf('Pooled accuracy vs Echange, split by EnoChange (%d sessions)', nSessions));
legend(makeNoChangeLabels(noChangeCuts), 'Location', 'best');
axis tight; box off;

figure(figBase+2); clf;
plot(binCtrNoChange, pErrByNoChange, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
xlabel('Unchanged-side evidence, EnoChange');
ylabel('P(error)');
title(sprintf('Pooled error vs EnoChange, middle %.0f%% of Echange', 100*midFrac));
axis tight; box off;

figure(figBase+3); clf; hold on;

% Change side: P(correct)
plot(binCtrChange, pCorrByChange, 'ko-', ...
    'LineWidth', 1.5, 'MarkerFaceColor', 'k');

% NoChange side: P(error)
plot(binCtrNoChange, pErrByNoChange, 'ro-', ...
    'LineWidth', 1.5, 'MarkerFaceColor', 'r');

xlabel('Evidence (Echange or EnoChange)');
ylabel('Probability');
title(sprintf('Change vs NoChange evidence effects (%d sessions)', nSessions));

legend({'P(correct) vs Echange', 'P(error) vs EnoChange'}, ...
    'Location', 'best');

ylim([0 1]);
box off;
yticks(0:0.05:1);
grid on;
ax = gca;
ax.YGrid = 'on';
ax.XGrid = 'off';

% ---- Output ----
out = struct();
out.alpha = alpha;
out.winMS = winMS;
out.nSessions = nSessions;
out.fileNames = allFileNames;
out.sessionOut = sessionOut;
out.sessionHeaders = sessionHeaders;

out.Echange = allEchange;
out.EnoChange = allEnoChange;
out.isCorrect = allCorrect;
out.sessionID = allSessionID;

out.plot1.binCtrChange = binCtrChange;
out.plot1.pCorrByChange = pCorrByChange;
out.plot1.nByChange = nByChange;
out.plot1.edgesChange = edgesChange;

out.plot2.noChangeCuts = noChangeCuts;
out.plot2.pCorrSplit = pCorrSplit;
out.plot2.nSplit = nSplit;

out.plot3.binCtrNoChange = binCtrNoChange;
out.plot3.pErrByNoChange = pErrByNoChange;
out.plot3.nByNoChange = nByNoChange;
out.plot3.edgesNoChange = edgesNoChange;
out.plot3.usePlot3 = usePlot3;

end

function [binCtr, yMean, nPerBin, edges, binIdx] = quantileBinnedMean(x, y, nBins)
x = x(:);
y = y(:);

ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);

edges = quantile(x, linspace(0,1,nBins+1));
edges(1) = -inf;
edges(end) = inf;

binIdx = discretize(x, edges);
binCtr = nan(nBins,1);
yMean = nan(nBins,1);
nPerBin = zeros(nBins,1);

for b = 1:nBins
    idx = (binIdx == b);
    nPerBin(b) = sum(idx);
    if nPerBin(b) > 0
        yMean(b) = mean(y(idx));
        binCtr(b) = mean(x(idx));
    end
end
end

function [groupIdx, cuts] = quantileGroups(x, nGroups)
x = x(:);
cuts = quantile(x, linspace(0,1,nGroups+1));
cuts(1) = -inf;
cuts(end) = inf;
groupIdx = discretize(x, cuts);
end

function labels = makeNoChangeLabels(cuts)
nG = numel(cuts)-1;
labels = cell(1,nG);
for g = 1:nG
    if g == 1
        labels{g} = 'low E_{noChange}';
    elseif g == nG
        labels{g} = 'high E_{noChange}';
    else
        labels{g} = 'mid E_{noChange}';
    end
end
end