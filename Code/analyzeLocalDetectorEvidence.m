function out = analyzeLocalDetectorEvidence(noiseFile, varargin)
% analyzeLocalDetectorEvidence
%
% First-pass diagnostic for a "local detectors + competitive comparison"
% decision architecture in the 2AFC IDReadout task.
%
% For increment trials only, computes:
%   Echange   = signed integrated evidence on changed side
%   EnoChange = signed integrated evidence on unchanged side
%
% using a scalar evidence model:
%   E = sum(pref(win)) + alpha * sum(probe(win))
%
% where alpha is the probe:pref scale to use for collapsing evidence.
%
% Then generates:
%   1) p(correct) vs Echange
%   2) p(correct) vs Echange, split by EnoChange tertiles
%   3) p(error) vs EnoChange, optionally restricted to middle Echange trials
%
% USAGE
%   out = analyzeLocalDetectorEvidence(noiseFile)
%   out = analyzeLocalDetectorEvidence(noiseFile, 'alpha', 0.19)
%   out = analyzeLocalDetectorEvidence(noiseFile, 'winMS', [0 125])
%
% INPUT
%   noiseFile : path to Data/NoiseMatrices/<base>.mat
%
% NAME-VALUE OPTIONS
%   'alpha'           : probe weight relative to pref (default 0.19)
%   'stepType'        : 2 for INC (default 2)
%   'winMS'           : [start stop] relative to step onset in ms (default [0 125])
%   'nBinsChange'     : number of quantile bins for Echange (default 7)
%   'nBinsNoChange'   : number of quantile bins for EnoChange plot (default 7)
%   'nGroupsNoChange' : number of quantile groups for split curves (default 3)
%   'restrictMidChangeForPlot3' : true/false (default true)
%   'midFrac'         : central fraction of Echange for Plot 3 (default 0.5)
%   'figNum'          : starting figure number (default 40)
%   'doPlot'          : true/false (default true)
%
% OUTPUT
%   out : struct with per-trial evidence and summary stats
%
% ASSUMPTIONS
%   NoiseMatrices file contains:
%     prefNoiseByPatch   [2 x nFrames x nTrials]
%     probeNoiseByPatch  [2 x nFrames x nTrials]
%     changeSidesAll     [1 x nTrials] or [nTrials x 1], 0=RF changed, 1=Opp changed
%     changeIndicesAll   [1 x nTrials] or [nTrials x 1], 1=DEC, 2=INC
%     trialOutcomesAll   [1 x nTrials] or [nTrials x 1], 0=correct, 1=wrong
%     header             struct with timing metadata
%
% NOTES
%   This function does not try to infer the "right" alpha. It just uses the
%   current preferred increment/change estimate as a first-pass scalarization.
%

% ---------------- Parse inputs ----------------
ip = inputParser;
ip.addRequired('noiseFile', @(x) ischar(x) || isstring(x));

ip.addParameter('alpha', 0.19, @isscalar);
ip.addParameter('stepType', 2, @isscalar);               % 2 = INC
ip.addParameter('winMS', [0 125], @(x) isnumeric(x) && numel(x)==2);
ip.addParameter('nBinsChange', 7, @isscalar);
ip.addParameter('nBinsNoChange', 7, @isscalar);
ip.addParameter('nGroupsNoChange', 3, @isscalar);
ip.addParameter('restrictMidChangeForPlot3', true, @islogical);
ip.addParameter('midFrac', 0.5, @(x) isscalar(x) && x>0 && x<=1);
ip.addParameter('figNum', 40, @isscalar);
ip.addParameter('doPlot', true, @islogical);

ip.parse(noiseFile, varargin{:});
P = ip.Results;

% ---------------- Load data ----------------
S = load([char(dataFolderPath()), '/NoiseMatrices/', noiseFile]);

reqFields = {'prefNoiseByPatch','probeNoiseByPatch','changeSidesAll', ...
             'changeIndicesAll','trialOutcomesAll','header'};
for i = 1:numel(reqFields)
    assert(isfield(S, reqFields{i}), 'Missing field: %s', reqFields{i});
end

prefNoiseByPatch  = S.prefNoiseByPatch;
probeNoiseByPatch = S.probeNoiseByPatch;
changeSidesAll    = S.changeSidesAll(:);
changeIndicesAll  = S.changeIndicesAll(:);
trialOutcomesAll  = S.trialOutcomesAll(:);
header            = S.header;

[nPatches, nFrames, nTrials] = size(prefNoiseByPatch);
assert(nPatches == 2, 'Expected 2 patches.');
assert(all(size(probeNoiseByPatch) == [2 nFrames nTrials]), ...
    'probeNoiseByPatch size mismatch.');

assert(numel(changeSidesAll)   == nTrials, 'changeSidesAll length mismatch.');
assert(numel(changeIndicesAll) == nTrials, 'changeIndicesAll length mismatch.');
assert(numel(trialOutcomesAll) == nTrials, 'trialOutcomesAll length mismatch.');

% ---------------- Timing / window ----------------
frameTimesMS = getFrameTimesMS(header, nFrames);
preStepMS = getPreStepMS(header);

winAbsMS = preStepMS + P.winMS;
winMask  = frameTimesMS >= winAbsMS(1) & frameTimesMS < winAbsMS(2);

assert(any(winMask), 'Integration window selects no frames.');

% ---------------- Trial selection ----------------
isDesiredStep = (changeIndicesAll == P.stepType);
isValidOutcome = ismember(trialOutcomesAll, [0 1]);  % 0=correct, 1=wrong

keep = isDesiredStep & isValidOutcome;

trialIdx = find(keep);
nKeep = numel(trialIdx);

Echange   = nan(nKeep,1);
EnoChange = nan(nKeep,1);
isCorrect = nan(nKeep,1);

prefChangeInt    = nan(nKeep,1);
probeChangeInt   = nan(nKeep,1);
prefNoChangeInt  = nan(nKeep,1);
probeNoChangeInt = nan(nKeep,1);

for k = 1:nKeep
    tr = trialIdx(k);

    % Convention assumed:
    %   changeSidesAll(tr) == 0 => RF patch changed => patch 1 changed
    %   changeSidesAll(tr) == 1 => Opp patch changed => patch 2 changed
    changedPatch = 1 + changeSidesAll(tr);
    unchangedPatch = 3 - changedPatch;

    prefChanged   = squeeze(prefNoiseByPatch(changedPatch,   :, tr));
    probeChanged  = squeeze(probeNoiseByPatch(changedPatch,  :, tr));
    prefUnchanged = squeeze(prefNoiseByPatch(unchangedPatch, :, tr));
    probeUnchanged= squeeze(probeNoiseByPatch(unchangedPatch,:, tr));

    % Integrate over window
    pC = sum(prefChanged(winMask));
    qC = sum(probeChanged(winMask));
    pN = sum(prefUnchanged(winMask));
    qN = sum(probeUnchanged(winMask));

    prefChangeInt(k)    = pC;
    probeChangeInt(k)   = qC;
    prefNoChangeInt(k)  = pN;
    probeNoChangeInt(k) = qN;

    % Scalar evidence
    %
    % Sign convention desired:
    %   larger Echange   = more evidence for correct side
    %   larger EnoChange = more misleading evidence on wrong side
    %
    % With your current stimulus conventions, positive pref/probe noise in
    % increment blocks should usually support "change-like" evidence.
    % If later you determine a sign flip is needed for one stream, this is
    % the place to do it.
    Echange(k)   = pC + P.alpha * qC;
    EnoChange(k) = pN + P.alpha * qN;

    isCorrect(k) = (trialOutcomesAll(tr) == 0);
end

% ---------------- Summaries / bins ----------------
[binCtrChange, pCorrByChange, nByChange, edgesChange, binIdxChange] = ...
    quantileBinnedMean(Echange, isCorrect, P.nBinsChange);

[groupIdxNoChange, noChangeCuts] = quantileGroups(EnoChange, P.nGroupsNoChange);

pCorrSplit = nan(P.nBinsChange, P.nGroupsNoChange);
nSplit     = nan(P.nBinsChange, P.nGroupsNoChange);

for g = 1:P.nGroupsNoChange
    for b = 1:P.nBinsChange
        idx = (groupIdxNoChange == g) & (binIdxChange == b);
        nSplit(b,g) = sum(idx);
        if nSplit(b,g) > 0
            pCorrSplit(b,g) = mean(isCorrect(idx));
        end
    end
end

% Plot 3
usePlot3 = true(size(Echange));
if P.restrictMidChangeForPlot3
    lo = quantile(Echange, (1-P.midFrac)/2);
    hi = quantile(Echange, 1-(1-P.midFrac)/2);
    usePlot3 = Echange >= lo & Echange <= hi;
end

[binCtrNoChange, pErrByNoChange, nByNoChange, edgesNoChange, binIdxNoChange] = ...
    quantileBinnedMean(EnoChange(usePlot3), ~isCorrect(usePlot3), P.nBinsNoChange);

% ---------------- Output struct ----------------
out = struct();
out.noiseFile = char(noiseFile);
out.params = P;

out.trialIdx = trialIdx;
out.Echange = Echange;
out.EnoChange = EnoChange;
out.isCorrect = isCorrect;

out.prefChangeInt    = prefChangeInt;
out.probeChangeInt   = probeChangeInt;
out.prefNoChangeInt  = prefNoChangeInt;
out.probeNoChangeInt = probeNoChangeInt;

out.frameTimesMS = frameTimesMS;
out.preStepMS = preStepMS;
out.winAbsMS = winAbsMS;
out.winMask = winMask;

out.plot1.binCtrChange = binCtrChange;
out.plot1.pCorrByChange = pCorrByChange;
out.plot1.nByChange = nByChange;
out.plot1.edgesChange = edgesChange;

out.plot2.groupIdxNoChange = groupIdxNoChange;
out.plot2.noChangeCuts = noChangeCuts;
out.plot2.pCorrSplit = pCorrSplit;
out.plot2.nSplit = nSplit;

out.plot3.usePlot3 = usePlot3;
out.plot3.binCtrNoChange = binCtrNoChange;
out.plot3.pErrByNoChange = pErrByNoChange;
out.plot3.nByNoChange = nByNoChange;
out.plot3.edgesNoChange = edgesNoChange;

% ---------------- Plotting ----------------
if P.doPlot
    % Plot 1
    figure(P.figNum); clf;
    plot(binCtrChange, pCorrByChange, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Changed-side evidence, Echange');
    ylabel('P(correct)');
    title(sprintf('Accuracy vs changed-side evidence (alpha = %.3f)', P.alpha));
    axis tight; box off;

    % Plot 2
    figure(P.figNum+1); clf; hold on;
    mk = {'o-','s-','d-'};
    for g = 1:P.nGroupsNoChange
        style = mk{min(g, numel(mk))};
        plot(binCtrChange, pCorrSplit(:,g), style, 'LineWidth', 1.5, 'MarkerSize', 7);
    end
    xlabel('Changed-side evidence, Echange');
    ylabel('P(correct)');
    title('Accuracy vs changed-side evidence, split by unchanged-side evidence');
    legend(makeNoChangeLabels(noChangeCuts), 'Location', 'best');
    axis tight; box off;

    % Plot 3
    figure(P.figNum+2); clf;
    plot(binCtrNoChange, pErrByNoChange, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Unchanged-side evidence, EnoChange');
    ylabel('P(error)');
    if P.restrictMidChangeForPlot3
        title(sprintf('Error vs unchanged-side evidence (middle %.0f%% of Echange)', 100*P.midFrac));
    else
        title('Error vs unchanged-side evidence');
    end
    axis tight; box off;
end

end

% ==================== Helpers ====================

function frameTimesMS = getFrameTimesMS(header, nFrames)

if ~isfield(header, 'frameRateHz')
    error('Could not determine frameRateHz from header.');
end

Hz = double(header.frameRateHz.data);
dt = 1000 / Hz;
frameTimesMS = (0:nFrames-1) * dt;

end

function preStepMS = getPreStepMS(header)

if ~isfield(header, 'preStepMS')
    error('Could not determine preStepMS from header.');
end

x = header.preStepMS;

if isstruct(x) && isfield(x, 'data')
    x = x.data;
end

x = double(x);
x = x(:)';   % row vector

if isempty(x)
    error('preStepMS is empty.');
end

% Many headers appear to store duplicated values, e.g. [750 750].
% If all entries agree, just take the first.
if numel(x) > 1
    if all(abs(x - x(1)) < 1e-9)
        preStepMS = x(1);
    else
        error('preStepMS contains multiple non-identical values: %s', mat2str(x));
    end
else
    preStepMS = x;
end

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