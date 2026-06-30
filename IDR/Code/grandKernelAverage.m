function grandKernelAverage(doBootstrap, nBoot, varargin)
% grandKernelAverage
%
% Compute session-averaged kernels from saved session noise matrices,
% compiling pref and probe kernels without regard for probe offset.
% Each session is reprocessed through computeSessionKernels(), then pooled
% across all sessions using inverse-variance weighting.
%

% cleanupObj = initProjectPath(); %#ok<NASGU>
if nargin < 1 || isempty(doBootstrap)
  doBootstrap = true;
end
if nargin < 2 || isempty(nBoot)
  nBoot = 500;
end

P = makeParser();
parse(P, varargin{:});
R0 = P.Results;
% check for nested file selection arguments and include them if they exist
if ~isempty(R0.FileSelectionArgs)
  topArgs = removeParameterPair(varargin, 'FileSelectionArgs');
  fileSelectionArgs = R0.FileSelectionArgs;
  P = makeParser();
  parse(P, topArgs{:}, fileSelectionArgs{:});
  R = P.Results;
else
  R = R0;
  fileSelectionArgs = R.FileSelectionArgs;
end

files = {};
fileInfo = table();
probeSessionDirs = dir([char(domainFolder(mfilename('fullpath'))), '/Data/Probe*']);
for p = 1:numel(probeSessionDirs)
    probeSessionDir = [probeSessionDirs(p).folder, '/', probeSessionDirs(p).name, '/ProbeSessions/'];
  [theFiles, theFileInfo] = selectAnalysisFiles(probeSessionDir, fileSelectionArgs{:});
  files = [files; theFiles]; %#ok<AGROW>
  fileInfo = [fileInfo; theFileInfo]; %#ok<AGROW>
end
[files, ia] = unique(files, 'stable');
fileInfo = fileInfo(ia, :);
if R.verbose
  nFiles = numel(fileInfo.fileName);
  if nFiles > 0
    fprintf('  Found %d sessions:\n', nFiles);
    fileNames = fileInfo.("fileName");
    for f = 1:nFiles
      fprintf('     %s\n', fileNames{f});
    end
    fprintf('     Total of %d sessions\n', nFiles);
  end
end
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
[firstPreStepMS] = integralWindowMS();

for f = 1:numel(files)
  filePath = files{f};
  load(filePath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr', 'prefNoiseByPatch', ...
        'probeNoiseByPatch', 'trialOutcomesAll', 'changeSidesAll', 'chosenSidesAll', 'changeIndicesAll');
  sessionData = struct;
  sessionData.sessionHeader = sessionHeader;
  sessionData.sessionProbeHeader = sessionProbeHeader;
  sessionData.sideTypeNames = sideTypeNames;
  sessionData.lr = lr;
  sessionData.prefNoiseByPatch = prefNoiseByPatch;
  sessionData.probeNoiseByPatch = probeNoiseByPatch;
  sessionData.trialOutcomesAll = trialOutcomesAll;
  sessionData.changeSidesAll = changeSidesAll;
  sessionData.changeIndicesAll = changeIndicesAll;

  [kernels, kVars, kStats, hitStats, compStats] = computeSessionKernels(sessionData);
  [kernelsNorm, kVarsNorm] = normalizeProbeKernelsToPrefAmplitude(kernels, kVars, sessionHeader, sessionProbeHeader);
  if ~initialized
    firstPreStepMS   = sessionHeader.preStepMS;
    firstStepMS   = sessionHeader.stepMS;
    firstVFrames  = size(kernels, 4);
    frameRateHz   = sessionHeader.frameRateHz;
    msPerVFrame   = 1000.0 / frameRateHz;
    initialized   = true;
  end
  preStepMS = sessionHeader.preStepMS;
  stepMS    = sessionHeader.stepMS;
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
  fprintf('No valid sessions found.\n');
return;
end

prefNoisePC = nan(nSessions,1);
for iSession = 1:nSessions
  prefNoisePC(iSession) = localDataValue(sessionDataList{iSession}.sessionHeader.prefCohNoisePC);
end
disp(unique(prefNoisePC));
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

avgCompStats.normInfo = normalizationInfoFromHeaders(sessionHeader, sessionProbeHeader);
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
    if b == 1 || mod(b, 25) == 0
      fprintf('      bootstrap %d of %d\n', b, nBoot);
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
      [bootKernelsNorm, bootKVarsNorm, ~] = normalizeProbeKernelsToPrefAmplitude(bootKernels, bootKVars, ...
        thisSession.sessionHeader, thisSession.sessionProbeHeader);
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
  avgCompStats.rawScaleCI.lo = prctile(bootRawScale, 2.5, 3);
  avgCompStats.rawScaleCI.hi = prctile(bootRawScale, 97.5, 3);
  avgCompStats.rawScaleBootSD = std(bootRawScale, 0, 3);

  avgCompStats.bootRawIntegrals = bootRawIntegrals;
  avgCompStats.rawIntegralCI.lo = prctile(bootRawIntegrals, 2.5, 4);
  avgCompStats.rawIntegralCI.hi = prctile(bootRawIntegrals, 97.5, 4);
  avgCompStats.rawIntegralBootSD = std(bootRawIntegrals, 0, 4);

  avgCompStats.bootNormScale = bootNormScale;
  avgCompStats.normScaleCI.lo = prctile(bootNormScale, 2.5, 3);
  avgCompStats.normScaleCI.hi = prctile(bootNormScale, 97.5, 3);
  avgCompStats.normScaleBootSD = std(bootNormScale, 0, 3);

  avgCompStats.bootNormIntegrals = bootNormIntegrals;
  avgCompStats.normIntegralCI.lo = prctile(bootNormIntegrals, 2.5, 4);
  avgCompStats.normIntegralCI.hi = prctile(bootNormIntegrals, 97.5, 4);
  avgCompStats.normIntegralBootSD = std(bootNormIntegrals, 0, 4);

  % Legacy aliases remain raw for now.
  avgCompStats.bootScale = bootRawScale;
  avgCompStats.scaleCI = avgCompStats.rawScaleCI;
  avgCompStats.scaleBootSD = avgCompStats.rawScaleBootSD;

  avgCompStats.bootKIntegrals = bootRawIntegrals;
  avgCompStats.kIntegralCI = avgCompStats.rawIntegralCI;
  avgCompStats.kIntegralBootSD = avgCompStats.rawIntegralBootSD;
end

plotName = sprintf('GrandAverageKernel_%s.pdf', R.Animal);
plotTitle = sprintf('\\bfGrand Kernel Averages Over %d Probe Directions, %d Probe Sessions (%s)', ...
        numel(probeSessionDirs), nSessions, R.Animal);

% ---- Plot/export averaged kernels ----
plotKernels(2, plotTitle, sessionHeader, avgKernels(1:5,:,:,:), ...
        avgKVars(1:5,:,:), avgCompStats, avgHitStats, R.probeDirDeg);
plotFolder = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Plots', 'AcrossProbes', 'Kernels'));
pdfFile = fullfile(plotFolder, plotName);
exportgraphics(gcf, pdfFile, 'ContentType', 'vector');
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
function [kernelsNorm, kVarsNorm, normInfo] = normalizeProbeKernelsToPrefAmplitude(kernels, kVars, sessionHeader, sessionProbeHeader)

normInfo = normalizationInfoFromHeaders(sessionHeader, sessionProbeHeader);
probeNormFactor = normInfo.probeNormFactor;

kernelsNorm = kernels;
kVarsNorm   = kVars;

kernelsNorm(:, :, 2, :) = kernelsNorm(:, :, 2, :) * probeNormFactor;
kVarsNorm(:, :, 2)      = kVarsNorm(:, :, 2) * probeNormFactor^2;
end

%%
function normInfo = normalizationInfoFromHeaders(sessionHeader, sessionProbeHeader)

prefCohNoisePC  = localDataValue(sessionHeader.prefCohNoisePC);
probeCohNoisePC = localDataValue(sessionProbeHeader.probeCohNoisePC);

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

function v = localDataValue(x)
if isstruct(x) && isfield(x, 'data')
  v = x.data;
else
  v = x;
end

if isnumeric(v) && ~isscalar(v)
  v = v(1);
end
end

% -------------------------------------------------------------------------
function p = makeParser()

p = inputParser;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'Bin179With180', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'FileSelectionArgs', {}, @(x) iscell(x));
addParameter(p, 'probeDirDeg', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'SummarySideType', 'change', @(x) ischar(x) || isstring(x));
addParameter(p, 'verbose', true, @islogical);
end

% -------------------------------------------------------------------------
function argsOut = removeParameterPair(argsIn, paramName)
argsOut = {};
k = 1;

while k <= numel(argsIn)
  if (ischar(argsIn{k}) || isstring(argsIn{k})) && strcmpi(char(argsIn{k}), paramName)
    k = k + 2;
  else
    argsOut{end+1} = argsIn{k}; %#ok<AGROW>
    k = k + 1;
  end
end
end
