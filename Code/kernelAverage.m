function kernelAverage(baseFolder, doBootstrap, nBoot)
% kernelAverage
%
% Compute session-averaged kernels from saved session noise matrices.
% Each session is reprocessed through computeSessionKernels(), then pooled
% across sessions using inverse-variance weighting.
%
% Optional hierarchical bootstrap:
%   1) resample sessions with replacement
%   2) resample trials within each selected session with replacement
%   3) recompute session kernels
%   4) recompute pooled kernels
%   5) recompute scale estimates
%
% The plotted kernel traces and SEM bands are the ordinary pooled estimates.
% The plotted scale CIs, when requested, are bootstrap percentile intervals.

if nargin < 1 || isempty(baseFolder)
  baseFolder = folderPath();
end
if nargin < 2 || isempty(doBootstrap)
  doBootstrap = false;
end
if nargin < 3 || isempty(nBoot)
  nBoot = 1000;
end

dataFolder = baseFolder + "/Data/NoiseMatrices/";
if ~exist(dataFolder, 'dir')
  error('kernelAverage:MissingFolder', 'Data folder not found: %s', dataFolder);
end

matFiles = dir(fullfile(dataFolder, '*.mat'));
[firstPreStepMS] = integralWindowMS();

% ---- Load valid sessions and recompute per-session kernels ----
sessionDataList  = {};
sessionKernels   = {};
sessionKVars     = {};
sessionKStats    = {};
sessionHitStats  = {};
sessionCompStats = {};
sessionHeaders   = {};
sideTypeNames = {};
initialized = false;
nSessions = 0;

for f = 1:length(matFiles)
  fileName = matFiles(f).name;
  if endsWith(fileName, '_fileInfo.mat')
    continue;
  end
  sessionData = load(fullfile(dataFolder, fileName));
  header = sessionData.header;
  if excludeFile(header)
    continue
  end

  % Skip sessions with 5% pref noise
  if header.prefNoiseCohPC.data ~= 10
    fprintf("Skipping   %s (%d of %d) -- prefNoiseCohPC is %.0f\n", ...
            fileName, f, length(matFiles), header.prefNoiseCohPC.data);
    continue;
  end

  fprintf("Processing %s (%d of %d)\n", fileName, f, length(matFiles));
  [kernels, kVars, kStats, hitStats, compStats] = computeSessionKernels(sessionData);

  if ~initialized
    firstStepMS   = header.stepMS.data(1);
    firstVFrames  = size(kernels, 4);
    frameRateHz   = header.frameRateHz.data(1);
    msPerVFrame   = 1000.0 / frameRateHz;
    initialized   = true;
  end
  if isempty(sideTypeNames) && isfield(sessionData, 'sideTypeNames')
    sideTypeNames = sessionData.sideTypeNames;
  end  
  preStepMS = header.preStepMS.data(1);
  stepMS    = header.stepMS.data(1);
  vFrames   = size(kernels, 4);

  if preStepMS ~= firstPreStepMS || stepMS ~= firstStepMS || vFrames ~= firstVFrames
    error('kernelAverage:IncompatibleFiles', ...
          'Files based on different trial or kernel lengths.');
  end

  sessionDataList{end+1}  = sessionData; %#ok<AGROW>
  sessionKernels{end+1}   = kernels; %#ok<AGROW>
  sessionKVars{end+1}     = kVars; %#ok<AGROW>
  sessionKStats{end+1}    = kStats; %#ok<AGROW>
  sessionHitStats{end+1}  = hitStats; %#ok<AGROW>
  sessionCompStats{end+1} = compStats; %#ok<AGROW>
  sessionHeaders{end+1}   = header; %#ok<AGROW>

  % Plot this session's kernels
  [~, sessBaseName, ~] = fileparts(header.fileName);
  plotKernels(1, sessBaseName, header, kernels, kVars, compStats, hitStats);

  nSessions = nSessions + 1;
end

if nSessions == 0
  error('kernelAverage:NoValidSessions', 'No valid sessions found.');
end

header = sessionHeaders{end};  % use last valid header for plotting title metadata

% ---- Pool ordinary session kernels ----
[avgKernels, avgKVars] = poolSessionKernels(sessionKernels, sessionKVars, firstVFrames);
avgHitStats = poolHitStats(sessionHitStats);

avgCompStats = struct;
[avgCompStats.kIntegrals, avgCompStats.R, avgCompStats.RVar] = ...
    kernelIntegral(avgKernels, avgKVars, msPerVFrame);
[avgCompStats.scale, avgCompStats.scaleSEM, avgCompStats.fitR2, avgCompStats.sse] = ...
    kernelScaleFit(avgKernels, msPerVFrame);

% Maintain prior behavior: convert ratio variance to SEM
avgCompStats.RVar = sqrt(avgCompStats.RVar);

% ---- Hierarchical bootstrap of scale CI ----
if doBootstrap
  rng(1);

  bootScale      = nan([size(avgCompStats.scale), nBoot]);
  bootKIntegrals = nan([size(avgCompStats.kIntegrals), nBoot]);

  for b = 1:nBoot
    if mod(b, 25) == 0
      fprintf('Bootstrap %d of %d\n', b, nBoot);
    end
    bootSessionIdx = randi(nSessions, [1 nSessions]);

    bootSessionKernels = cell(1, nSessions);
    bootSessionKVars   = cell(1, nSessions);

    for j = 1:nSessions
      iSession = bootSessionIdx(j);
      thisSession = sessionDataList{iSession};

      nTrials = size(thisSession.prefNoiseByPatch, 3);
      trialIdx = randi(nTrials, [1 nTrials]);

      [bootKernels, bootKVars] = computeSessionKernels(thisSession, trialIdx);

      bootSessionKernels{j} = bootKernels;
      bootSessionKVars{j}   = bootKVars;
    end

    [bootAvgKernels, bootAvgKVars] = poolSessionKernels(bootSessionKernels, bootSessionKVars, firstVFrames);

    bootCompStats = struct;
    [bootCompStats.kIntegrals, bootCompStats.R, bootCompStats.RVar] = ...
        kernelIntegral(bootAvgKernels, bootAvgKVars, msPerVFrame);
    [bootCompStats.scale, bootCompStats.scaleSEM, bootCompStats.fitR2, bootCompStats.sse] = ...
        kernelScaleFit(bootAvgKernels, msPerVFrame);

    bootScale(:,:,b) = bootCompStats.scale;
    bootKIntegrals(:,:,:,b) = bootCompStats.kIntegrals;
  end

  avgCompStats.bootScale = bootScale;
  avgCompStats.scaleCI.lo = prctile(bootScale, 15.865, 3);
  avgCompStats.scaleCI.hi = prctile(bootScale, 84.135, 3);
  avgCompStats.scaleBootSD = std(bootScale, 0, 3);

  avgCompStats.bootKIntegrals = bootKIntegrals;
  avgCompStats.kIntegralCI.lo = prctile(bootKIntegrals, 15.865, 4);
  avgCompStats.kIntegralCI.hi = prctile(bootKIntegrals, 84.135, 4);
  avgCompStats.kIntegralBootSD = std(bootKIntegrals, 0, 4);
end

% ---- Plot/export averaged kernels ----
plotKernels(2, sprintf('%d Session Average', nSessions), header, avgKernels, avgKVars, ...
  avgCompStats, avgHitStats);

pdfFile = fullfile(baseFolder, 'Plots', 'Kernels', 'Latest Session Average Kernel.pdf');
exportgraphics(gcf, pdfFile, 'ContentType', 'vector');
fprintf(' Saved session average kernel: %s\n', pdfFile);

end


function [avgKernels, avgKVars] = poolSessionKernels(sessionKernels, sessionKVars, nFrames)
% Pool session kernels using inverse-variance weighting.

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