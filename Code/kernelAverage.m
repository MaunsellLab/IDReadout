function kernelAverage(baseFolder)

% excludedFiles = {'IDReadout_Meetz_20260113', 'IDReadout_Meetz_20260114', 'IDReadout_Meetz_20260114_2', ...
%   'IDReadout_Meetz_20260114_3', 'IDReadout_Meetz_20260114_4', 'IDReadout_Meetz_20260115', ...
%   'IDReadout_Meetz_20260116', 'IDReadout_Meetz_20260209', 'IDReadout_Meetz_20260210'};

excludedFiles = {};

if nargin < 1 || isempty(baseFolder)
  baseFolder = folderPath();
end
dataFolder = baseFolder + "/Kernels/";
if ~exist(dataFolder, 'dir')
  mkdir(dataFolder);
end
summaryFolder = baseFolder + "/Summaries/";
if ~exist(summaryFolder, 'dir')
  mkdir(summaryFolder);
end

matFiles = dir(fullfile(dataFolder, '*.mat'));
[firstPreStepMS] = integralWindowMS();
initialized = false;

% ---- Grand-pooled accumulators ----
avgHitStats = struct('nTrials', zeros(1,2), 'nHits', zeros(1,2), 'nRFTrials', zeros(1,2), 'nRFHits', zeros(1,2));    

for f = 1:length(matFiles)
  fileName = matFiles(f).name;
  [~, baseName, ~] = fileparts(fileName);  % ignore path, take name
  if endsWith(fileName, '_fileInfo.mat')
    continue;
  end
  if exist('excludedFiles', 'var') && ~isempty(excludedFiles) && ismember(baseName, excludedFiles)
    continue;
  end
  S = load(dataFolder + fileName);   % load variables
  header = S.header;
  % Skip the sessions where we ran with 5% noise.  No signal there and it
  % distorts the SEM.
  if header.prefNoiseCohPC.data ~= 10
    fprintf("Skipping   %s (%d of %d) -- prefNoiseCohPC is %.0f\n", ...
      fileName, f, length(matFiles), header.prefNoiseCohPC.data);
    continue;
  end
  fprintf("Processing %s (%d of %d)\n", fileName, f, length(matFiles));

  kernels = S.kernels;
  kVars = S.kVars;
  compStats = S.compStats;
  hitStats = S.hitStats;
  if ~initialized
    firstStepMS = header.stepMS.data(1);
    firstVFrames = size(kernels, 4);

    % ---- Session-averaged (difference-of-means) accumulators ----
    avgKernels = zeros(5, 2, 2, firstVFrames);
    % avgKInts = zeros(5, 2, 2);
    sumWeights = zeros(5, 2, 2);             % summed inverse-variance weights over sessions
    nHits = zeros(5, 2);
    nTrials = zeros(5, 2);
    % sumRWeights = zeros(5, 1, 2);            % we weight R by its inverse variance
    % sumWeightedR = zeros(5, 1, 2);
    nSessions = 0;

    % ---- Grand-pooled accumulators ----
    grandSumCorrect = cell(5, 2, 2);
    grandSumWrong   = cell(5, 2, 2);
    grandNCorrect   = zeros(5, 2, 2);
    grandNWrong     = zeros(5, 2, 2);

    for sideType = 1:5
      for s = 1:2
        for p = 1:2
          grandSumCorrect{sideType, s, p} = zeros(firstVFrames, 1);
          grandSumWrong{sideType, s, p}   = zeros(firstVFrames, 1);
        end
      end
    end
    initialized = true;
  end

  frameRateHz = header.frameRateHz.data(1);
  preStepMS = header.preStepMS.data(1);       % sometimes header values are replicated vectors
  msPerVFrame = 1000.0 / frameRateHz;
  stepMS = header.stepMS.data(1);
  vFrames = size(kernels, 4);
  if preStepMS ~= firstPreStepMS || stepMS ~= firstStepMS || vFrames ~= firstVFrames
    error('Files based on different trial or kernel lengths.');
  end

  % ---- Session-averaged accumulation (inverse-variance weighted kernels) ----
  for sideType = 1:5
    for s = 1:2
      for p = 1:2
        % accumulate weighted kernels only if variance is valid
        if isfinite(kVars(sideType, s, p)) && kVars(sideType, s, p) > 0
          w = 1.0 / kVars(sideType, s, p);
          avgKernels(sideType, s, p, :) = avgKernels(sideType, s, p, :) + kernels(sideType, s, p, :) * w;
          sumWeights(sideType, s, p) = sumWeights(sideType, s, p) + w;
        end
      end
    end
  end

  % ---- Grand-pooled accumulation ----
  kStats = S.kStats;
  for sideType = 1:5
    for s = 1:2
      for p = 1:2
        grandSumCorrect{sideType, s, p} = grandSumCorrect{sideType, s, p} + kStats(sideType, s, p).sumCorrect;
        grandSumWrong{sideType, s, p}   = grandSumWrong{sideType, s, p}   + kStats(sideType, s, p).sumWrong;
        grandNCorrect(sideType, s, p)   = grandNCorrect(sideType, s, p) + kStats(sideType, s, p).nCorrect;
        grandNWrong(sideType, s, p)     = grandNWrong(sideType, s, p)   + kStats(sideType, s, p).nWrong;
      end
    end
  end
  for s = 1:2
    avgHitStats.nTrials(s) = avgHitStats.nTrials(s) + hitStats.nTrials(s);
    avgHitStats.nHits(s) = avgHitStats.nHits(s) + hitStats.nHits(s);
    avgHitStats.nRFTrials(s) = avgHitStats.nRFTrials(s) + hitStats.nRFTrials(s);
    avgHitStats.nRFHits(s) = avgHitStats.nRFHits(s) + hitStats.nRFHits(s);
  end
  % plot this session's kernels
  [~, baseName, ~] = fileparts(header.fileName);
  plotKernels(1, baseName, header, kernels, kVars, compStats, hitStats);
  nSessions = nSessions + 1;
end % end of file loop

% ---- Compute and display session-average from weighted-average kernels ----
avgKVars = nan(5, 2, 2);
for sideType = 1:5
  for s = 1:2
    for p = 1:2
      if sumWeights(sideType, s, p) > 0
        avgKernels(sideType, s, p, :) = avgKernels(sideType, s, p, :) / sumWeights(sideType, s, p);
        % variance of the weighted-average kernel estimate
        avgKVars(sideType, s, p) = 1.0 / sumWeights(sideType, s, p);
      else
        avgKernels(sideType, s, p, :) = nan;
        avgKVars(sideType, s, p) = nan;
      end
    end
  end
end

% Integrals and ratios for the averaged kernels
avgCompStats = struct;
[avgCompStats.kIntegrals, avgCompStats.R, avgCompStats.RVar] = kernelIntegral(avgKernels, avgKVars, msPerVFrame);
[avgCompStats.scale, avgCompStats.scaleSEM, avgCompStats.fitR2, avgCompStats.sse] = ...
                kernelScaleFit(avgKernels, msPerVFrame);

% SEM should be sqrt(variance), not variance itself
for sideType = 1:5
  for s = 1:2
    avgCompStats.RVar = sqrt(avgCompStats.RVar);
  end
end

plotKernels(2, sprintf('%d Session Average (Diff-of-Means)', nSessions), header, avgKernels, avgKVars, ...
  avgCompStats, avgHitStats);
pdfFile = fullfile(baseFolder, 'Plots', 'Kernels', ' Latest Session Average Kernel.pdf');
exportgraphics(gcf, pdfFile, 'ContentType', 'vector');
fprintf(' Saved session average kernel: %s\n', pdfFile);

% ---- Compute and display grand-pooled kernel ----
% Grand-pooled variance uses sigma^2 from the stimulus noise distribution.
% makeKernels/meanPsychKernel now saves stats.sigma2 for each (s,p).
sigma2Hat = nan(5, 2, 2);
% Prefer sigma2 saved by meanPsychKernel (should be constant across sessions).
% Use sigma2 from the most recently loaded session file (S) as an estimate.
if isfield(S, 'kStats') && isstruct(S.kStats)
  for sideType = 1:5
    for s = 1:2
      for p = 1:2
        if isfield(S.kStats(sideType, s, p), 'sigma2')
          sigma2Hat(sideType, s, p) = S.kStats(sideType, s, p).sigma2;
        end
      end
    end
  end
end

grandKernels = nan(5, 2, 2, firstVFrames);
grandKVars   = nan(5, 2, 2);
for sideType = 1:5
  for s = 1:2
    for p = 1:2
      nC = grandNCorrect(sideType, s, p);
      nW = grandNWrong(sideType, s, p);
      if nC > 0 && nW > 0
        meanC = grandSumCorrect{sideType, s, p} ./ nC;
        meanW = grandSumWrong{sideType, s, p}   ./ nW;
        grandKernels(sideType, s, p, :) = meanC - meanW;  
        if ~isnan(sigma2Hat(sideType, s, p))
          grandKVars(sideType, s, p) = sigma2Hat(sideType, s, p) * (1/nC + 1/nW);
        else
          grandKVars(sideType, s, p) = nan;
        end
      end
    end
  end
end

% Compute the integrals and regression scaling for the grand kernels
compStats = struct;
[compStats.kIntegrals, compStats.R, compStats.RVar] = kernelIntegral(grandKernels, grandKVars, msPerVFrame);
[compStats.scale, compStats.scaleSEM, compStats.fitR2, compStats.sse] = kernelScaleFit(grandKernels, msPerVFrame);

plotKernels(3, sprintf("Grand-Pooled Trial Average (%d sessions)", nSessions), header, grandKernels, grandKVars, ...
  compStats, avgHitStats);

pdfFile = fullfile(baseFolder, 'Plots', 'Kernels', ' Latest Trial Average Kernel.pdf');
exportgraphics(gcf, pdfFile, 'ContentType', 'vector');
fprintf(' Saved trials average kernel: %s\n', pdfFile);
end
