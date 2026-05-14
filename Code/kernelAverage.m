function kernelAverage(doBootstrap, nBoot, varargin)
% kernelAverage
%
% Compute session-averaged kernels from saved session noise matrices.
% Each session is reprocessed through computeSessionKernels(), then pooled
% across sessions using inverse-variance weighting.
%
%   kernelAverage(folderPath, false, 1000, ...
%       'dataFolder', fullfile(folderPath, 'Data', 'probe45', 'NoiseMatrices'), ...
%       'plotFolder', fullfile(folderPath, 'Plots', 'AverageKernels'));

if nargin < 1 || isempty(doBootstrap)
  doBootstrap = false;
end
if nargin < 2 || isempty(nBoot)
  nBoot = 1000;
end

baseFolder = folderPath();
P = inputParser;
addParameter(P, 'dataFolder', '', @(x) ischar(x) || isstring(x));
addParameter(P, 'plotFolder', '', @(x) ischar(x) || isstring(x));
addParameter(P, 'probeDirDeg', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(P, 'SummarySideType', 'change', @(x) ischar(x) || isstring(x));
parse(P, varargin{:});
R = P.Results;
summarySideType = char(R.SummarySideType);
summarySideTypeNum = sideTypeIndex(summarySideType);

if isempty(R.dataFolder)
  dataFolder = fullfile(baseFolder, 'Data', 'NoiseMatrices');
else
  dataFolder = char(R.dataFolder);
end

if isempty(R.plotFolder)
  plotFolder = fullfile(baseFolder, 'Plots', 'Kernels');
else
  plotFolder = char(R.plotFolder);
end
if ~exist(dataFolder, 'dir')
  error('kernelAverage:MissingFolder', 'Data folder not found: %s', dataFolder);
end
if ~exist(plotFolder, 'dir')
  mkdir(plotFolder);
end

matFiles = dir(fullfile(dataFolder, '*.mat'));
[firstPreStepMS] = integralWindowMS();

% ---- Load valid sessions and recompute per-session kernels ----
sessionDataList  = {};
sessionKernels   = {};
sessionKVars     = {};
sessionKernelsNorm = {};
sessionKVarsNorm   = {};
sessionKStats    = {};
sessionHitStats  = {};
sessionCompStats = {};
sessionProbeHeaders = {};
initialized = false;
nSessions = 0;

for f = 1:length(matFiles)
  fileName = matFiles(f).name;
  if endsWith(fileName, '_fileInfo.mat')
    continue;
  end
  sessionData = load(fullfile(dataFolder, fileName));
  
  assert(isfield(sessionData, 'sessionProbeHeader'), ...
    'kernelAverage:MissingSessionProbeHeader', ...
    'Expected sessionProbeHeader in %s.', fileName);
  
  sessionProbeHeader = sessionData.sessionProbeHeader;
  
  [kernels, kVars, kStats, hitStats, compStats] = computeSessionKernels(sessionData);
  [kernelsNorm, kVarsNorm] = normalizeProbeKernelsToPrefAmplitude(kernels, kVars, sessionProbeHeader);
  if ~initialized
    firstStepMS   = sessionProbeHeader.stepMS.data(1);
    firstVFrames  = size(kernels, 4);
    frameRateHz   = sessionProbeHeader.frameRateHz.data(1);
    msPerVFrame   = 1000.0 / frameRateHz;
    initialized   = true;
  end
  preStepMS = sessionProbeHeader.preStepMS.data(1);
  stepMS    = sessionProbeHeader.stepMS.data(1);
  vFrames   = size(kernels, 4);

  if preStepMS ~= firstPreStepMS || stepMS ~= firstStepMS || vFrames ~= firstVFrames
    error('kernelAverage:IncompatibleFiles', 'Files based on different trial or kernel lengths.');
  end

  sessionDataList{end+1}  = sessionData; %#ok<AGROW>
  sessionKernels{end+1}   = kernels; %#ok<AGROW>
  sessionKVars{end+1}     = kVars; %#ok<AGROW>
  sessionKernelsNorm{end+1} = kernelsNorm; %#ok<AGROW>
  sessionKVarsNorm{end+1}   = kVarsNorm; %#ok<AGROW>  
  sessionKStats{end+1}    = kStats; %#ok<AGROW>
  sessionHitStats{end+1}  = hitStats; %#ok<AGROW>
  sessionCompStats{end+1} = compStats; %#ok<AGROW>
  sessionProbeHeaders{end+1} = sessionProbeHeader; %#ok<AGROW>
  nSessions = nSessions + 1;
end

if nSessions == 0
  error('kernelAverage:NoValidSessions', 'No valid sessions found.');
end

sessionProbeHeader = sessionProbeHeaders{end};  % use last valid metadata record for plotting

% ---- Pool session kernels ----
%
% Raw pooled kernels are used for plotting the measured kernel traces.
% Normalized pooled kernels are used for reported integrals, ratios, and
% scales. Normalization rescales the probe stream to the pref-noise amplitude
% convention before pooling.
[avgKernels, avgKVars] = poolSessionKernels(sessionKernels, sessionKVars, firstVFrames);
[avgKernelsNorm, avgKVarsNorm] = poolSessionKernels(sessionKernelsNorm, sessionKVarsNorm, firstVFrames);
avgHitStats = poolHitStats(sessionHitStats);

avgCompStats = struct;

% Raw comparison statistics
[avgCompStats.rawIntegrals, avgCompStats.rawR, avgCompStats.rawRVar] = ...
    kernelIntegral(avgKernels, avgKVars, msPerVFrame);
[avgCompStats.rawScale, avgCompStats.rawScaleSEM, avgCompStats.rawFitR2, avgCompStats.rawSSE] = ...
    kernelScaleFit(avgKernels, msPerVFrame);

% Normalized comparison statistics
[avgCompStats.normIntegrals, avgCompStats.normR, avgCompStats.normRVar] = ...
    kernelIntegral(avgKernelsNorm, avgKVarsNorm, msPerVFrame);
[avgCompStats.normScale, avgCompStats.normScaleSEM, avgCompStats.normFitR2, avgCompStats.normSSE] = ...
    kernelScaleFit(avgKernelsNorm, msPerVFrame);

% Maintain prior behavior: convert ratio variance to SEM
avgCompStats.rawRVar  = sqrt(avgCompStats.rawRVar);
avgCompStats.normRVar = sqrt(avgCompStats.normRVar);

avgCompStats.normInfo = normalizationInfoFromSessionProbeHeader(sessionProbeHeader);

% Legacy aliases: preserve old downstream behavior for now.
avgCompStats.kIntegrals = avgCompStats.rawIntegrals;
avgCompStats.R          = avgCompStats.rawR;
avgCompStats.RVar       = avgCompStats.rawRVar;
avgCompStats.scale      = avgCompStats.rawScale;
avgCompStats.scaleSEM   = avgCompStats.rawScaleSEM;
avgCompStats.fitR2      = avgCompStats.rawFitR2;
avgCompStats.sse        = avgCompStats.rawSSE;

% ---- Hierarchical bootstrap of scale CI ----
if doBootstrap
  rng(1);
  bootRawScale      = nan([size(avgCompStats.rawScale), nBoot]);
  bootRawIntegrals  = nan([size(avgCompStats.rawIntegrals), nBoot]);

  bootNormScale     = nan([size(avgCompStats.normScale), nBoot]);
  bootNormIntegrals = nan([size(avgCompStats.normIntegrals), nBoot]);
  for b = 1:nBoot
    if mod(b, 25) == 0
      fprintf('     bootstrap %d of %d\n', b, nBoot);
    end
    bootSessionIdx = randi(nSessions, [1 nSessions]);
    bootSessionKernels = cell(1, nSessions);
    bootSessionKVars   = cell(1, nSessions);
    bootSessionKernelsNorm = cell(1, nSessions);
    bootSessionKVarsNorm   = cell(1, nSessions);
    for j = 1:nSessions
      iSession = bootSessionIdx(j);
      thisSession = sessionDataList{iSession};

      nTrials = size(thisSession.prefNoiseByPatch, 3);
      trialIdx = randi(nTrials, [1 nTrials]);
      [bootKernels, bootKVars] = computeSessionKernels(thisSession, trialIdx);
      [bootKernelsNorm, bootKVarsNorm] = ...
              normalizeProbeKernelsToPrefAmplitude(bootKernels, bootKVars, thisSession.sessionProbeHeader);
      bootSessionKernels{j} = bootKernels;
      bootSessionKVars{j}   = bootKVars;

      bootSessionKernelsNorm{j} = bootKernelsNorm;
      bootSessionKVarsNorm{j}   = bootKVarsNorm;
    end
    [bootAvgKernels, bootAvgKVars] = ...
        poolSessionKernels(bootSessionKernels, bootSessionKVars, firstVFrames);
    [bootAvgKernelsNorm, bootAvgKVarsNorm] = ...
        poolSessionKernels(bootSessionKernelsNorm, bootSessionKVarsNorm, firstVFrames);

    bootCompStats = struct;

    [bootCompStats.rawIntegrals, ~, ~] = ...
        kernelIntegral(bootAvgKernels, bootAvgKVars, msPerVFrame);
    [bootCompStats.rawScale, ~, ~, ~] = ...
        kernelScaleFit(bootAvgKernels, msPerVFrame);

    [bootCompStats.normIntegrals, ~, ~] = ...
        kernelIntegral(bootAvgKernelsNorm, bootAvgKVarsNorm, msPerVFrame);
    [bootCompStats.normScale, ~, ~, ~] = ...
        kernelScaleFit(bootAvgKernelsNorm, msPerVFrame);

    bootRawScale(:,:,b)       = bootCompStats.rawScale;
    bootRawIntegrals(:,:,:,b) = bootCompStats.rawIntegrals;

    bootNormScale(:,:,b)       = bootCompStats.normScale;
    bootNormIntegrals(:,:,:,b) = bootCompStats.normIntegrals;
  end
  avgCompStats.bootRawScale = bootRawScale;
  avgCompStats.rawScaleCI.lo = prctile(bootRawScale, 15.865, 3);
  avgCompStats.rawScaleCI.hi = prctile(bootRawScale, 84.135, 3);
  avgCompStats.rawScaleBootSD = std(bootRawScale, 0, 3);

  avgCompStats.bootRawIntegrals = bootRawIntegrals;
  avgCompStats.rawIntegralCI.lo = prctile(bootRawIntegrals, 15.865, 4);
  avgCompStats.rawIntegralCI.hi = prctile(bootRawIntegrals, 84.135, 4);
  avgCompStats.rawIntegralBootSD = std(bootRawIntegrals, 0, 4);

  avgCompStats.bootNormScale = bootNormScale;
  avgCompStats.normScaleCI.lo = prctile(bootNormScale, 15.865, 3);
  avgCompStats.normScaleCI.hi = prctile(bootNormScale, 84.135, 3);
  avgCompStats.normScaleBootSD = std(bootNormScale, 0, 3);

  avgCompStats.bootNormIntegrals = bootNormIntegrals;
  avgCompStats.normIntegralCI.lo = prctile(bootNormIntegrals, 15.865, 4);
  avgCompStats.normIntegralCI.hi = prctile(bootNormIntegrals, 84.135, 4);
  avgCompStats.normIntegralBootSD = std(bootNormIntegrals, 0, 4);

  % Legacy aliases remain raw for now.
  avgCompStats.bootScale = bootRawScale;
  avgCompStats.scaleCI = avgCompStats.rawScaleCI;
  avgCompStats.scaleBootSD = avgCompStats.rawScaleBootSD;

  avgCompStats.bootKIntegrals = bootRawIntegrals;
  avgCompStats.kIntegralCI = avgCompStats.rawIntegralCI;
  avgCompStats.kIntegralBootSD = avgCompStats.rawIntegralBootSD;
end

% ---- Determine probe direction for naming ----
if isempty(R.probeDirDeg)
  if isfield(sessionProbeHeader, 'probeDirDeg')
    probeDirDeg = sessionProbeHeader.probeDirDeg;
  else
    probeDirDeg = [];
  end
else
  probeDirDeg = R.probeDirDeg;
end

if isempty(probeDirDeg)
  plotName = '_AverageKernel.pdf';
  plotTitle = sprintf('%d Session Average', nSessions);
else
  probeTag = sprintf('probe%d', round(probeDirDeg));
  plotName = sprintf('AverageKernel_%s.pdf', probeTag);
  plotTitle = sprintf('\\bf%d° Probe Kernels %d Session Average', probeDirDeg, nSessions);
end

% ---- Plot/export averaged kernels ----
plotKernels(2, plotTitle, sessionProbeHeader, avgKernels(1:5,:,:,:), ...
  avgKVars(1:5,:,:,:), avgCompStats, avgHitStats, probeDirDeg);
pdfFile = fullfile(plotFolder, plotName);
exportgraphics(gcf, pdfFile, 'ContentType', 'vector');

% ---- Save kernel data for summary display of one side type across probe directions ----
averageKernelPlotData = struct();

averageKernelPlotData.summarySideType = summarySideType;
averageKernelPlotData.sideTypeNum = summarySideTypeNum;
averageKernelPlotData.probeDirDeg = probeDirDeg;

averageKernelPlotData.kernels = squeeze(avgKernels(summarySideTypeNum,:,:,:));
averageKernelPlotData.kVars   = squeeze(avgKVars(summarySideTypeNum,:,:,:));

[~, averageKernelPlotData.sideTypeNames] = sideTypeIndex();
averageKernelPlotData.stepTypeNames = {'inc', 'dec'};      % confirm order if needed
averageKernelPlotData.streamTypeNames = {'pref', 'probe'};

averageKernelPlotData.tMS = ((1:firstVFrames) - 1) * msPerVFrame - firstPreStepMS;
averageKernelPlotData.firstPreStepMS = firstPreStepMS;
averageKernelPlotData.firstStepMS = firstStepMS;
averageKernelPlotData.msPerVFrame = msPerVFrame;

averageKernelPlotData.avgCompStats = avgCompStats;
averageKernelPlotData.avgHitStats = avgHitStats;
averageKernelPlotData.nSessions = nSessions;

averageKernelPlotData.meta = struct();
averageKernelPlotData.meta.createdDate = datetime('now');

summaryDataFolder = fullfile(baseFolder, 'Data', probeTag, 'AverageKernels', summarySideType);
validFolder(summaryDataFolder);

save(fullfile(summaryDataFolder, 'AverageKernelPlotData.mat'), 'averageKernelPlotData');

end

%% Pool session kernels using inverse-variance weighting.
function [avgKernels, avgKVars] = poolSessionKernels(sessionKernels, sessionKVars, nFrames)

nSessions = numel(sessionKernels);
nSideTypes = size(sessionKernels{1}, 1);

avgKernels = zeros(nSideTypes, 2, 2, nFrames);
sumWeights = zeros(nSideTypes, 2, 2);

for iSession = 1:nSessions
  kernels = sessionKernels{iSession};
  kVars   = sessionKVars{iSession};

  for sideType = 1:nSideTypes
    for stepType = 1:2
      for streamType = 1:2
        if isfinite(kVars(sideType, stepType, streamType)) && ...
            kVars(sideType, stepType, streamType) > 0
          w = 1.0 / kVars(sideType, stepType, streamType);
          avgKernels(sideType, stepType, streamType, :) = ...
              avgKernels(sideType, stepType, streamType, :) + ...
              kernels(sideType, stepType, streamType, :) * w;
          sumWeights(sideType, stepType, streamType) = ...
              sumWeights(sideType, stepType, streamType) + w;
        end
      end
    end
  end
end

avgKVars = nan(nSideTypes, 2, 2);
for sideType = 1:nSideTypes
  for stepType = 1:2
    for streamType = 1:2
      if sumWeights(sideType, stepType, streamType) > 0
        avgKernels(sideType, stepType, streamType, :) = ...
            avgKernels(sideType, stepType, streamType, :) ./ ...
            sumWeights(sideType, stepType, streamType);
        avgKVars(sideType, stepType, streamType) = ...
            1.0 / sumWeights(sideType, stepType, streamType);
      else
        avgKernels(sideType, stepType, streamType, :) = nan;
        avgKVars(sideType, stepType, streamType) = nan;
      end
    end
  end
end
end

%%
function avgHitStats = poolHitStats(sessionHitStats)
% Sum hit/trial statistics across sessions.

avgHitStats = struct;
avgHitStats.nTrials    = [0 0];
avgHitStats.nHits      = [0 0];
avgHitStats.nRFTrials  = [0 0];
avgHitStats.nRFHits    = [0 0];
avgHitStats.nLeftTrials = [0 0];
avgHitStats.nLeftHits   = [0 0];

for iSession = 1:numel(sessionHitStats)
  hs = sessionHitStats{iSession};
  avgHitStats.nTrials     = avgHitStats.nTrials     + hs.nTrials;
  avgHitStats.nHits       = avgHitStats.nHits       + hs.nHits;
  avgHitStats.nRFTrials   = avgHitStats.nRFTrials   + hs.nRFTrials;
  avgHitStats.nRFHits     = avgHitStats.nRFHits     + hs.nRFHits;

  if isfield(hs, 'nLeftTrials')
    avgHitStats.nLeftTrials = avgHitStats.nLeftTrials + hs.nLeftTrials;
    avgHitStats.nLeftHits   = avgHitStats.nLeftHits   + hs.nLeftHits;
  end
end
end

%%
function [kernelsNorm, kVarsNorm, normInfo] = normalizeProbeKernelsToPrefAmplitude(kernels, kVars, sessionProbeHeader)
% Normalize probe-stream kernels to the pref-noise amplitude convention.
%
% The raw kernels are measured using each stream's actual coherence-noise
% amplitude. Under the standard small-signal assumption, kernel amplitude
% scales with noise variance. Therefore, converting a probe kernel measured
% at probeCohNoisePC to the equivalent prefCohNoisePC convention requires:
%
%   K_probe_norm = K_probe_raw * (prefCohNoisePC / probeCohNoisePC)^2

normInfo = normalizationInfoFromSessionProbeHeader(sessionProbeHeader);
probeNormFactor = normInfo.probeNormFactor;

kernelsNorm = kernels;
kVarsNorm   = kVars;

kernelsNorm(:, :, 2, :) = kernelsNorm(:, :, 2, :) * probeNormFactor;
kVarsNorm(:, :, 2)      = kVarsNorm(:, :, 2) * probeNormFactor^2;
end

%%
function normInfo = normalizationInfoFromSessionProbeHeader(sessionProbeHeader)

prefCohNoisePC  = sessionProbeHeader.prefCohNoisePC.data(1);
probeCohNoisePC = sessionProbeHeader.probeCohNoisePC.data(1);

if ~isfinite(prefCohNoisePC) || ~isfinite(probeCohNoisePC) || probeCohNoisePC <= 0
  error('kernelAverage:BadNoiseAmplitude', ...
    'Invalid pref/probe coherence noise amplitudes: pref=%g, probe=%g.', ...
    prefCohNoisePC, probeCohNoisePC);
end

normInfo = struct();
normInfo.prefCohNoisePC = prefCohNoisePC;
normInfo.probeCohNoisePC = probeCohNoisePC;

nYokedProbeStreams = probeStreamCountFromSessionProbeHeader(sessionProbeHeader);
combinedProbeCohNoisePC = nYokedProbeStreams * probeCohNoisePC;

normInfo.nYokedProbeStreams = nYokedProbeStreams;
normInfo.combinedProbeCohNoisePC = combinedProbeCohNoisePC;
normInfo.probeNormFactor = (prefCohNoisePC / combinedProbeCohNoisePC)^2;
normInfo.method = ...
  'probe kernels multiplied by (prefCohNoisePC/(nYokedProbeStreams*probeCohNoisePC))^2 before computing normalized ratios/scales';
end

%%
function n = probeStreamCountFromSessionProbeHeader(sessionProbeHeader)

assert(isfield(sessionProbeHeader, 'probeDirDeg'), ...
  'kernelAverage:MissingProbeDir', ...
  'Cannot determine probe stream count because sessionProbeHeader.probeDirDeg is missing.');

probeDirDeg = abs(double(sessionProbeHeader.probeDirDeg));

if probeDirDeg > 0 && probeDirDeg < 180
  n = 2;
elseif abs(probeDirDeg - 180) < 1e-9
  n = 1;
else
  error('kernelAverage:UnsupportedProbeDir', ...
    'Unsupported probeDirDeg for probe normalization: %g.', probeDirDeg);
end
end