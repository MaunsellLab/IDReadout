function rectGainByDirection = fitIDQRectGainByDirection( ...
  trialTable, sessionFits, alignedWeibull, targetPerformance, directionLabels)
% fitIDQRectGainByDirection
%
% Diagnostic-only rectangular noise gain fits by physical drift direction.

nDirs = numel(directionLabels);

dirIndex = (1:nDirs)';
directionLabel = directionLabels(:);

gain = nan(nDirs, 1);
SE = nan(nDirs, 1);
CI95Low = nan(nDirs, 1);
CI95High = nan(nDirs, 1);
zVsFlat = nan(nDirs, 1);
pVsFlat = nan(nDirs, 1);
zVsZero = nan(nDirs, 1);
pVsZero = nan(nDirs, 1);
nTrials = nan(nDirs, 1);
nCorrect = nan(nDirs, 1);
meanXNoiseCorrect = nan(nDirs, 1);
meanXNoiseError = nan(nDirs, 1);
meanDiffCorrectMinusError = nan(nDirs, 1);
nEffectiveCohClipped = nan(nDirs, 1);

fitByDirection = cell(nDirs, 1);

for iDir = 1:nDirs

  idxDir = trialTable.dirIndex == iDir;

  Tdir = trialTable(idxDir, :);

  fit = fitIDQNoiseGain( ...
    Tdir, ...
    'rectNoisePredictor', ...
    sessionFits, ...
    alignedWeibull, ...
    targetPerformance);

  fit.directionIndex = iDir;
  fit.directionLabel = directionLabels(iDir);

  fitByDirection{iDir} = fit;

  gain(iDir) = fit.gain;
  SE(iDir) = fit.SE;
  CI95Low(iDir) = fit.CI95(1);
  CI95High(iDir) = fit.CI95(2);
  zVsFlat(iDir) = fit.zVsFlat;
  pVsFlat(iDir) = fit.pVsFlat;
  zVsZero(iDir) = fit.zVsZero;
  pVsZero(iDir) = fit.pVsZero;
  nTrials(iDir) = fit.nTrials;
  nCorrect(iDir) = fit.nCorrect;
  meanXNoiseCorrect(iDir) = fit.meanXNoiseCorrect;
  meanXNoiseError(iDir) = fit.meanXNoiseError;
  meanDiffCorrectMinusError(iDir) = fit.meanDiffCorrectMinusError;

  if isfield(fit, 'nEffectiveCohClipped')
    nEffectiveCohClipped(iDir) = fit.nEffectiveCohClipped;
  end
end

summaryTable = table( ...
  dirIndex, ...
  directionLabel, ...
  nTrials, ...
  nCorrect, ...
  gain, ...
  SE, ...
  CI95Low, ...
  CI95High, ...
  zVsFlat, ...
  pVsFlat, ...
  zVsZero, ...
  pVsZero, ...
  meanXNoiseCorrect, ...
  meanXNoiseError, ...
  meanDiffCorrectMinusError, ...
  nEffectiveCohClipped);

rectGainByDirection = struct();
rectGainByDirection.predictor = 'rectNoisePredictor';
rectGainByDirection.flatPrediction = 1;
rectGainByDirection.summaryTable = summaryTable;
rectGainByDirection.fitByDirection = fitByDirection;

end