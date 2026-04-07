function summary = collectKernelSessionSummary(baseFolder)
% collectKernelSessionSummary
%
% Build a per-session summary table for all kernel classes.
%
% One row is created for each:
%   session x sideType x stepType
%
% Stored quantities are intended to support:
%   1) averaging session scales
%   2) pref-energy weighting
%   3) numerator/denominator aggregation
%   4) pooled-vs-session comparisons
%   5) later diff-vs-change decomposition
%   6) kernel shape coherence analyses
%
% Assumptions:
%   - Saved session files are in Data/NoiseMatrices/*.mat
%   - sessionKernelFromSaved(sessionData) returns:
%       [kernels, kVars, kStats, hitStats, compStats]
%   - kernels has shape: 5 x 2 x 2 x nFrames
%       sideType: 1=diff, 2=change, 3=noChange, 4=RF, 5=Opp
%       stepType: 1=inc, 2=dec
%       streamType: 1=pref, 2=probe
%
% Output:
%   summary is a table with one row per session x sideType x stepType.
%
% Notes:
%   - This function does not pool across sessions.
%   - It computes session-level quantities over the current integration window.
%   - It stores the integration-window kernel vectors in cell arrays so they
%     can be reused later without recomputing.
%
% John can later add extra columns if useful.

if nargin < 1 || isempty(baseFolder)
  baseFolder = folderPath();
end

dataFolder = fullfile(baseFolder, 'Data', 'NoiseMatrices');
if ~exist(dataFolder, 'dir')
  error('collectKernelSessionSummary:MissingFolder', ...
        'Data folder not found: %s', dataFolder);
end

matFiles = dir(fullfile(dataFolder, '*.mat'));
if isempty(matFiles)
  error('collectKernelSessionSummary:NoFiles', ...
        'No .mat files found in %s', dataFolder);
end

% Get integration window definition from existing code
[firstPreStepMS, intStartMS, intDurMS] = integralWindowMS();

% Row accumulator
rows = struct( ...
  'sessionIndex',   {}, ...
  'fileName',       {}, ...
  'sessionName',    {}, ...
  'sideType',       {}, ...
  'stepType',       {}, ...
  'sideLabel',      {}, ...
  'stepLabel',      {}, ...
  'nTrialsStep',    {}, ...
  'nHitsStep',      {}, ...
  'nRFTrialsStep',  {}, ...
  'nRFHitsStep',    {}, ...
  'prefVar',        {}, ...
  'probeVar',       {}, ...
  'poolWeightPref', {}, ...
  'poolWeightProbe',{}, ...
  'msPerVFrame',    {}, ...
  'preStepMS',      {}, ...
  'stepMS',         {}, ...
  'intStartMS',     {}, ...
  'intDurMS',       {}, ...
  'intStartFrame',  {}, ...
  'intEndFrame',    {}, ...
  'prefVec',        {}, ...
  'probeVec',       {}, ...
  'prefIntegral',   {}, ...
  'probeIntegral',  {}, ...
  'prefEnergy',     {}, ...
  'probeEnergy',    {}, ...
  'crossDot',       {}, ...
  'sessionScale',   {}, ...
  'cosSim',         {}, ...
  'fitR2',          {}, ...
  'sse',            {} );

sideLabels = {'diff', 'change', 'noChange', 'RF', 'Opp'};
stepLabels = {'inc', 'dec'};

rowCounter = 0;
validSessionCounter = 0;

firstStepMS = [];
firstNFrames = [];
firstFrameRateHz = [];

for f = 1:length(matFiles)
  fileName = matFiles(f).name;

  if endsWith(fileName, '_fileInfo.mat')
    continue;
  end

  sessionData = load(fullfile(dataFolder, fileName));
  header = sessionData.header;

  if excludeFile(header)
    continue;
  end

  % Maintain current kernelAverage behavior: skip sessions with 5% pref noise
  if header.prefNoiseCohPC.data ~= 10
    fprintf('Skipping %s -- prefNoiseCohPC is %.0f\n', ...
      fileName, header.prefNoiseCohPC.data);
    continue;
  end

  fprintf('Collecting %s (%d of %d)\n', fileName, f, length(matFiles));

  [kernels, kVars, ~, hitStats, compStats] = sessionKernelFromSaved(sessionData);

  preStepMS   = header.preStepMS.data(1);
  stepMS      = header.stepMS.data(1);
  frameRateHz = header.frameRateHz.data(1);
  msPerVFrame = 1000 / frameRateHz;
  nFrames     = size(kernels, 4);

  if isempty(firstStepMS)
    firstStepMS     = stepMS;
    firstNFrames    = nFrames;
    firstFrameRateHz = frameRateHz;
  else
    if preStepMS ~= firstPreStepMS || stepMS ~= firstStepMS || nFrames ~= firstNFrames
      error('collectKernelSessionSummary:IncompatibleFiles', ...
            'Incompatible session lengths or timing in %s.', fileName);
    end
  end

  % Integration window frames
  % Assume kernel time base starts at stimulus onset and preStepMS marks step onset.
  intStartFrame = round((preStepMS + intStartMS) / msPerVFrame) + 1;
  intEndFrame   = round((preStepMS + intStartMS + intDurMS) / msPerVFrame);

  intStartFrame = max(1, intStartFrame);
  intEndFrame   = min(nFrames, intEndFrame);

  if intEndFrame < intStartFrame
    error('collectKernelSessionSummary:BadWindow', ...
          'Integration window is empty for %s.', fileName);
  end

  validSessionCounter = validSessionCounter + 1;

  [~, sessionName, ~] = fileparts(fileName);

  for sideType = 1:5
    for stepType = 1:2

      prefVec  = squeeze(kernels(sideType, stepType, 1, intStartFrame:intEndFrame));
      probeVec = squeeze(kernels(sideType, stepType, 2, intStartFrame:intEndFrame));

      prefVec  = prefVec(:);
      probeVec = probeVec(:);

      prefIntegral = sum(prefVec) * msPerVFrame;
      probeIntegral = sum(probeVec) * msPerVFrame;

      prefEnergy  = sum(prefVec .^ 2);
      probeEnergy = sum(probeVec .^ 2);
      crossDot    = sum(prefVec .* probeVec);

      if prefEnergy > 0
        sessionScale = crossDot / prefEnergy;
      else
        sessionScale = nan;
      end

      if prefEnergy > 0 && probeEnergy > 0
        cosSim = crossDot / sqrt(prefEnergy * probeEnergy);
      else
        cosSim = nan;
      end

      prefVar  = kVars(sideType, stepType, 1);
      probeVar = kVars(sideType, stepType, 2);

      if isfinite(prefVar) && prefVar > 0
        poolWeightPref = 1 / prefVar;
      else
        poolWeightPref = nan;
      end

      if isfinite(probeVar) && probeVar > 0
        poolWeightProbe = 1 / probeVar;
      else
        poolWeightProbe = nan;
      end

      % Pull session-level fit diagnostics if present
      fitR2 = nan;
      sse   = nan;
      if isfield(compStats, 'fitR2')
        fitR2 = compStats.fitR2(sideType, stepType);
      end
      if isfield(compStats, 'sse')
        sse   = compStats.sse(sideType, stepType);
      end

      rowCounter = rowCounter + 1;

      rows(rowCounter).sessionIndex    = validSessionCounter;
      rows(rowCounter).fileName        = string(fileName);
      rows(rowCounter).sessionName     = string(sessionName);

      rows(rowCounter).sideType        = sideType;
      rows(rowCounter).stepType        = stepType;
      rows(rowCounter).sideLabel       = string(sideLabels{sideType});
      rows(rowCounter).stepLabel       = string(stepLabels{stepType});

      rows(rowCounter).nTrialsStep     = hitStats.nTrials(stepType);
      rows(rowCounter).nHitsStep       = hitStats.nHits(stepType);
      rows(rowCounter).nRFTrialsStep   = hitStats.nRFTrials(stepType);
      rows(rowCounter).nRFHitsStep     = hitStats.nRFHits(stepType);

      rows(rowCounter).prefVar         = prefVar;
      rows(rowCounter).probeVar        = probeVar;
      rows(rowCounter).poolWeightPref  = poolWeightPref;
      rows(rowCounter).poolWeightProbe = poolWeightProbe;

      rows(rowCounter).msPerVFrame     = msPerVFrame;
      rows(rowCounter).preStepMS       = preStepMS;
      rows(rowCounter).stepMS          = stepMS;
      rows(rowCounter).intStartMS      = intStartMS;
      rows(rowCounter).intDurMS        = intDurMS;
      rows(rowCounter).intStartFrame   = intStartFrame;
      rows(rowCounter).intEndFrame     = intEndFrame;

      rows(rowCounter).prefVec         = {prefVec};
      rows(rowCounter).probeVec        = {probeVec};

      rows(rowCounter).prefIntegral    = prefIntegral;
      rows(rowCounter).probeIntegral   = probeIntegral;
      rows(rowCounter).prefEnergy      = prefEnergy;
      rows(rowCounter).probeEnergy     = probeEnergy;
      rows(rowCounter).crossDot        = crossDot;
      rows(rowCounter).sessionScale    = sessionScale;
      rows(rowCounter).cosSim          = cosSim;

      rows(rowCounter).fitR2           = fitR2;
      rows(rowCounter).sse             = sse;
    end
  end
end

if isempty(rows)
  error('collectKernelSessionSummary:NoValidSessions', ...
        'No valid sessions found.');
end

summary = struct2table(rows);

% Helpful sort order
summary = sortrows(summary, {'sideType', 'stepType', 'sessionIndex'});

fprintf('\nCollected %d rows from %d valid sessions.\n', ...
        height(summary), validSessionCounter);

% Quick compact summary
fprintf('Integration window: start = %.1f ms, duration = %.1f ms\n', ...
        intStartMS, intDurMS);
fprintf('Frames: %d to %d\n', summary.intStartFrame(1), summary.intEndFrame(1));

for sideType = 1:5
  for stepType = 1:2
    idx = summary.sideType == sideType & summary.stepType == stepType;
    if any(idx)
      fprintf('  %-8s %-3s : nSessions=%2d, meanScale=% .3f, prefEwMean=% .3f\n', ...
        sideLabels{sideType}, stepLabels{stepType}, ...
        sum(idx), ...
        mean(summary.sessionScale(idx), 'omitnan'), ...
        sum(summary.prefEnergy(idx) .* summary.sessionScale(idx), 'omitnan') / ...
        sum(summary.prefEnergy(idx), 'omitnan'));
    end
  end
end

end