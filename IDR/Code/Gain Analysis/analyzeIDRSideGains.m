function inventory = analyzeIDRSideGains(varargin)
% analyzeIDRSideGains  Analyze change- and no-change-side gains in IDR.
%
% inventory = analyzeIDRSideGains(...)
%
% First implementation stage:
%   Inventory valid INC trials in parent FullSession files, separately for
%   noise and no-noise trials and by nominal step coherence.
%
% The function reads only sessionHeader from each selected FullSession file.
% It does not load the large trials cell array.
%
% Name-value arguments:
%   Animal          Animal name or 'All' (default)
%   ReportExcluded  Forwarded to selectAnalysisFiles (default false)
%   PrintInventory  Print session-by-session counts (default true)
%   SaveSummary      Save the analysis MAT-file (default true)
%   FitPsychometrics  Fit fixed-shape session psychometrics (default true)
%   PlotPsychometrics Export a multipage diagnostic PDF (default true)
%   BetaWeibull      Fixed Weibull exponent (default 2)
%   Lapse            Fixed lapse rate (default 0.01)
%   TargetPerformance Threshold performance (default 0.75)
%   AlignedLapseBounds Bounds for pooled aligned lapse (default [0 0.05])
%   SensitivityMaxRelativeCIWidth Threshold-fit uncertainty cutoff used only
%                     for a secondary sensitivity fit (default 0.75)
%   RefinePsychometrics Refit session thresholds using pooled beta/lapse
%                     estimates (default true)
%   DiagnosticThirdPass Perform a second refinement (three total passes)
%                     to assess convergence (default true)
%   FitSideGains      Assemble probeSession trials and fit side gains
%                     (default true)
%   GainBounds       Lower/upper bounds for gain fits (default [-5 5])
%   AnalyzeCandidateOverlap  Compute and save rectangular IDR candidate-
%                     overlap summaries (default true)
%   FitMatchedGains  Fit four absolute side-by-stream gains separately at
%                     every probe offset (default true)
%   FitSignedGains   Fit offset-specific signed and opponent-constrained
%                     gain models (default true)
%   PreferredNullRatio Fixed preferred:null sensitivity ratio used by the
%                     opponent-constrained models (default 2.42)
%   FitControlledPooling Compare opponent-rectified no-change sum, p-norm,
%                     and hard-max models (default true)
%   PoolingPBounds     Bounds for fitted p (default [1 100])
%   PoolingProfilePValues Fixed-p profile grid (default
%                     [1 1.5 2 3 5 10 20 50 100])
%   ExcludedProbeDirectionsDeg Probe offsets omitted from assembled noise
%                     analyses. The sole 1-deg session is excluded by
%                     default (default 1)
%
% OUTPUT
%   inventory
%       createdAt
%       createdBy
%       animalSelection
%       fullSessionFolder
%       fileInfo
%       trialTable
%       conditionSummary
%       sessionSummary
%
% trialTable contains one row per valid INC trial across selected sessions.
%
% conditionSummary contains one row per:
%   daily session x hasStepNoise x  stepDeltaCohPC
%
% sessionSummary contains overall INC counts and no-noise coherence coverage
% for each daily session.

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) ischar(x) || isstring(x));
addParameter(p, 'ReportExcluded', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'PrintInventory', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'SaveSummary', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'FitPsychometrics', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'PlotPsychometrics', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'BetaWeibull', 2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'Lapse', 0.01, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x < 0.5);
addParameter(p, 'TargetPerformance', 0.75, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0.5 && x < 1);
addParameter(p, 'AlignedLapseBounds', [0 0.05], ...
  @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && ...
  x(1) >= 0 && x(2) > x(1) && x(2) < 0.5);
addParameter(p, 'SensitivityMaxRelativeCIWidth', 0.75, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'RefinePsychometrics', true, ...
  @(x) islogical(x) && isscalar(x));
addParameter(p, 'DiagnosticThirdPass', true, ...
  @(x) islogical(x) && isscalar(x));
addParameter(p, 'FitSideGains', true, ...
  @(x) islogical(x) && isscalar(x));
addParameter(p, 'GainBounds', [-5 5], ...
  @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && x(1) < x(2));
addParameter(p, 'AnalyzeCandidateOverlap', true, ...
  @(x) islogical(x) && isscalar(x));
addParameter(p, 'FitMatchedGains', true, ...
  @(x) islogical(x) && isscalar(x));
addParameter(p, 'FitSignedGains', true, ...
  @(x) islogical(x) && isscalar(x));
addParameter(p, 'PreferredNullRatio', 2.42, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'FitControlledPooling', true, ...
  @(x) islogical(x) && isscalar(x));
addParameter(p, 'PoolingPBounds', [1 100], ...
  @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && ...
  x(1) >= 1 && x(2) > x(1));
addParameter(p, 'PoolingProfilePValues', [1 1.5 2 3 5 10 20 50 100], ...
  @(x) isnumeric(x) && isvector(x) && ~isempty(x) && ...
  all(isfinite(x)) && all(x >= 1));
addParameter(p, 'ExcludedProbeDirectionsDeg', 1, ...
  @(x) isnumeric(x) && isvector(x) && all(isfinite(x)));
parse(p, varargin{:});
opts = p.Results;

domainPath = domainFolder(mfilename('fullpath'));
fullSessionFolder = fullfile(domainPath, 'Data', 'FullSessions');
summaryFolder = validFolder(fullfile(domainPath, 'Data', 'AcrossSessionSummaries'));
plotFolder = validFolder(fullfile(domainPath, 'Plots', 'AcrossSessionSummaries'));

[filePaths, fileInfo] = selectAnalysisFiles( ...
  fullSessionFolder, ...
  'Animal', opts.Animal, ...
  'ReportExcluded', opts.ReportExcluded);

if isempty(filePaths)
  error('analyzeIDRSideGains:NoFiles', ...
    'No FullSession files were selected from %s.', fullSessionFolder);
end

sessionTables = cell(numel(filePaths), 1);

for iFile = 1:numel(filePaths)
  filePath = filePaths{iFile};
  S = load(filePath, 'sessionHeader');

  if ~isfield(S, 'sessionHeader')
    error('analyzeIDRSideGains:MissingSessionHeader', ...
      'File %s does not contain sessionHeader.', filePath);
  end

  H = S.sessionHeader;
  if ~isfield(H, 'behaviorTrialTable')
    error('analyzeIDRSideGains:MissingBehaviorTrialTable', ...
      ['sessionHeader.behaviorTrialTable is missing from %s. ' ...
       'Run updateIDRSessionHeaders first.'], filePath);
  end

  T = H.behaviorTrialTable;
  requireBehaviorTrialVariables(T, filePath);

  % Primary analysis is restricted to valid increment trials.
  idxUse = T.isValid & T.changeIndex == 2;
  T = T(idxUse, :);

  n = height(T);
  T.sessionIndex = repmat(iFile, n, 1);
  T.sessionID = repmat(string(H.fileName), n, 1);
  T.animal = repmat(string(H.animal), n, 1);
  T.filePath = repmat(string(filePath), n, 1);

  % Keep identity variables at the front.
  T = movevars(T, ...
    {'sessionIndex', 'sessionID', 'animal', 'filePath'}, ...
    'Before', 1);

  sessionTables{iFile} = T;
end

trialTable = vertcat(sessionTables{:});

conditionSummary = makeConditionSummary(trialTable);
sessionSummary = makeSessionSummary(trialTable, fileInfo);

inventory = struct();
inventory.createdAt = datetime('now');
inventory.createdBy = mfilename;
inventory.animalSelection = string(opts.Animal);
inventory.fullSessionFolder = fullSessionFolder;
inventory.fileInfo = fileInfo;
inventory.trialTable = trialTable;
inventory.conditionSummary = conditionSummary;
inventory.sessionSummary = sessionSummary;

if opts.FitPsychometrics
  % Pass 1: provisional fixed-shape thresholds.
  pass1SessionFits = fitSessionPsychometrics( ...
    trialTable, sessionSummary, opts.BetaWeibull, opts.Lapse, opts.TargetPerformance);
  fprintf('\nPass 1: provisional session fits (fixed beta %.3f, lapse %.4f)\n', ...
    opts.BetaWeibull, opts.Lapse);
  printPsychometricSummary(pass1SessionFits);

  [pass1Aligned, pass1TrialTable] = fitAlignedPsychometrics( ...
    trialTable, pass1SessionFits, opts.TargetPerformance, ...
    opts.AlignedLapseBounds, opts.SensitivityMaxRelativeCIWidth);
  fprintf('\nPass 1: provisional alignment\n');
  printAlignedPsychometricSummary(pass1Aligned);

  % Defaults when no refinement is requested.
  sessionPsychFits = pass1SessionFits;
  alignedPsychometric = pass1Aligned;
  trialTable = pass1TrialTable;
  sessionFitBeta = opts.BetaWeibull;
  sessionFitLapse = opts.Lapse;
  nPasses = 1;

  pass2SessionFits = [];
  pass2Aligned = [];
  pass3SessionFits = [];
  pass3Aligned = [];

  if opts.RefinePsychometrics
    % Pass 2: first refinement using the pass-1 pooled beta and lapse.
    pass2SessionBeta = pass1Aligned.primary.betaWeibull;
    pass2SessionLapse = pass1Aligned.primary.lapse;
    pass2SessionFits = fitSessionPsychometrics( ...
      trialTable, sessionSummary, pass2SessionBeta, pass2SessionLapse, ...
      opts.TargetPerformance);

    fprintf('\nPass 2: refined session fits (fixed beta %.3f, lapse %.4f)\n', ...
      pass2SessionBeta, pass2SessionLapse);
    printPsychometricSummary(pass2SessionFits);

    [pass2Aligned, pass2TrialTable] = fitAlignedPsychometrics( ...
      trialTable, pass2SessionFits, opts.TargetPerformance, ...
      opts.AlignedLapseBounds, opts.SensitivityMaxRelativeCIWidth);
    fprintf('\nPass 2: refined alignment\n');
    printAlignedPsychometricSummary(pass2Aligned);
    printPsychometricPassChange(pass1Aligned.primary, pass2Aligned.primary, 1, 2);

    sessionPsychFits = pass2SessionFits;
    alignedPsychometric = pass2Aligned;
    trialTable = pass2TrialTable;
    sessionFitBeta = pass2SessionBeta;
    sessionFitLapse = pass2SessionLapse;
    nPasses = 2;

    if opts.DiagnosticThirdPass
      % Pass 3: second refinement to test convergence.
      pass3SessionBeta = pass2Aligned.primary.betaWeibull;
      pass3SessionLapse = pass2Aligned.primary.lapse;
      pass3SessionFits = fitSessionPsychometrics( ...
        trialTable, sessionSummary, pass3SessionBeta, pass3SessionLapse, ...
        opts.TargetPerformance);

      fprintf('\nPass 3: diagnostic session fits (fixed beta %.3f, lapse %.4f)\n', ...
        pass3SessionBeta, pass3SessionLapse);
      printPsychometricSummary(pass3SessionFits);

      [pass3Aligned, pass3TrialTable] = fitAlignedPsychometrics( ...
        trialTable, pass3SessionFits, opts.TargetPerformance, ...
        opts.AlignedLapseBounds, opts.SensitivityMaxRelativeCIWidth);
      fprintf('\nPass 3: diagnostic alignment\n');
      printAlignedPsychometricSummary(pass3Aligned);
      printPsychometricPassChange(pass2Aligned.primary, pass3Aligned.primary, 2, 3);

      sessionPsychFits = pass3SessionFits;
      alignedPsychometric = pass3Aligned;
      trialTable = pass3TrialTable;
      sessionFitBeta = pass3SessionBeta;
      sessionFitLapse = pass3SessionLapse;
      nPasses = 3;
    end
  end

  % Preserve all passes and expose the final pass through the established
  % unqualified fields used by subsequent analyses.
  inventory.sessionPsychFitsPass1 = pass1SessionFits;
  inventory.alignedPsychometricPass1 = pass1Aligned;
  inventory.sessionPsychFitsProvisional = pass1SessionFits;
  inventory.alignedPsychometricProvisional = pass1Aligned;

  if ~isempty(pass2SessionFits)
    inventory.sessionPsychFitsPass2 = pass2SessionFits;
    inventory.alignedPsychometricPass2 = pass2Aligned;
  end
  if ~isempty(pass3SessionFits)
    inventory.sessionPsychFitsPass3 = pass3SessionFits;
    inventory.alignedPsychometricPass3 = pass3Aligned;
  end

  inventory.trialTable = trialTable;
  inventory.sessionPsychFits = sessionPsychFits;
  inventory.alignedPsychometric = alignedPsychometric;
  inventory.psychometricRefinementHistory = makePsychometricRefinementHistory( ...
    opts.BetaWeibull, opts.Lapse, pass1Aligned, ...
    pass2Aligned, pass3Aligned);
  inventory.psychometricSettings = struct( ...
    'provisionalBetaWeibull', opts.BetaWeibull, ...
    'provisionalLapse', opts.Lapse, ...
    'finalSessionFitBetaWeibull', sessionFitBeta, ...
    'finalSessionFitLapse', sessionFitLapse, ...
    'nPasses', nPasses, ...
    'refinementPerformed', opts.RefinePsychometrics, ...
    'diagnosticThirdPassPerformed', nPasses == 3, ...
    'targetPerformance', opts.TargetPerformance, ...
    'trialSelection', 'valid INC no-noise trials only');

  fprintf('\nPsychometric refinement history\n');
  disp(inventory.psychometricRefinementHistory);

  if opts.PlotPsychometrics
    psychPlotPath = fullfile(plotFolder, 'IDR_SideGainSessionPsychometrics.pdf');
    plotSessionPsychometrics(conditionSummary, sessionPsychFits, psychPlotPath, ...
      sessionFitBeta, sessionFitLapse, opts.TargetPerformance);
    inventory.psychometricPlotPath = psychPlotPath;
    fprintf('Saved psychometric diagnostics: %s\n', psychPlotPath);

    alignedPlotPath = fullfile(plotFolder, 'IDR_SideGainAlignedPsychometric.pdf');
    plotAlignedPsychometric(alignedPsychometric, alignedPlotPath);
    inventory.alignedPsychometricPlotPath = alignedPlotPath;
    fprintf('Saved final aligned psychometric: %s\n', alignedPlotPath);

    pass1PlotPath = fullfile( ...
      plotFolder, 'IDR_SideGainAlignedPsychometric_Pass1.pdf');
    plotAlignedPsychometric(pass1Aligned, pass1PlotPath);
    inventory.alignedPsychometricPass1PlotPath = pass1PlotPath;
    fprintf('Saved pass-1 aligned psychometric: %s\n', pass1PlotPath);

    if ~isempty(pass2Aligned)
      pass2PlotPath = fullfile( ...
        plotFolder, 'IDR_SideGainAlignedPsychometric_Pass2.pdf');
      plotAlignedPsychometric(pass2Aligned, pass2PlotPath);
      inventory.alignedPsychometricPass2PlotPath = pass2PlotPath;
      fprintf('Saved pass-2 aligned psychometric: %s\n', pass2PlotPath);
    end
  end
end

if opts.FitSideGains
  if ~opts.FitPsychometrics
    error('analyzeIDRSideGains:PsychometricsRequired', ...
      'FitSideGains requires FitPsychometrics=true.');
  end

  fprintf('\nAssembling IDR probeSession side-gain trials...\n');
  sideGainData = assembleIDRSideGainTrials( ...
    domainPath, inventory.trialTable, inventory.sessionPsychFits, opts.Animal, ...
    opts.ExcludedProbeDirectionsDeg);

  fprintf('  ProbeSession files: %d\n', height(sideGainData.probeFileInfo));
  fprintf('  Parent sessions: %d\n', numel(unique(sideGainData.trialTable.sessionIndex)));
  fprintf('  INC noise trials: %d\n', height(sideGainData.trialTable));

  sideGainPredictors = computeIDRSideGainPredictors(sideGainData);

  if opts.AnalyzeCandidateOverlap
    candidateOverlap = analyzeIDRCandidateOverlap(sideGainData);
    inventory.candidateOverlap = candidateOverlap;

    fprintf('\nIDR rectangular candidate-overlap summary\n');
    displayVars = { ...
      'probeDirDeg','nTrials','nSessions', ...
      'sdChangePreferred','sdChangeProbe','corrChangePrefProbe', ...
      'meanStep','meanStepOverSDChangePreferred', ...
      'meanStepOverSDChangeProbe','q95ChangeMargin', ...
      'q99ChangeMargin','maxChangeMargin','pChangeCrossed'};
    disp(candidateOverlap.offsetSummary(:,displayVars));

    overlapOffsetPath = fullfile( ...
      summaryFolder, 'IDR_CandidateOverlap_ByOffset.csv');
    overlapSessionPath = fullfile( ...
      summaryFolder, 'IDR_CandidateOverlap_BySessionOffset.csv');
    overlapConditionPath = fullfile( ...
      summaryFolder, 'IDR_CandidateOverlap_ByCondition.csv');
    writetable(candidateOverlap.offsetSummary, overlapOffsetPath);
    writetable(candidateOverlap.sessionOffsetSummary, overlapSessionPath);
    writetable(candidateOverlap.conditionSummary, overlapConditionPath);
    inventory.candidateOverlapTablePaths = string([ ...
      overlapOffsetPath; overlapSessionPath; overlapConditionPath]);
  end

  fprintf('\nIDR side-gain predictor diagnostics\n');
  disp(sideGainPredictors.predictorDiagnostics);

  kernelPlotPath = fullfile(plotFolder, 'IDR_SideGainCommonTemporalKernel.pdf');
  plotIDRSideGainCommonKernel(sideGainPredictors, kernelPlotPath);
  inventory.sideGainKernelPlotPath = kernelPlotPath;
  fprintf('Saved common temporal kernel: %s\n', kernelPlotPath);

  finalBeta = inventory.alignedPsychometric.primary.betaWeibull;
  finalLapse = inventory.alignedPsychometric.primary.lapse;
  sideGainFits = fitIDRSideGainModels( ...
    sideGainData.trialTable, sideGainPredictors, ...
    finalBeta, finalLapse, opts.TargetPerformance, opts.GainBounds);

  inventory.sideGainData = sideGainData;
  inventory.sideGainPredictors = sideGainPredictors;
  inventory.sideGainFits = sideGainFits;

  if opts.FitMatchedGains
    matchedGainPlotPath = fullfile( ...
      plotFolder, 'IDR_MatchedAbsoluteGains_ByOffset.pdf');
    matchedGainFits = analyzeIDRMatchedGains( ...
      sideGainData, finalBeta, finalLapse, opts.TargetPerformance, ...
      opts.GainBounds, 'PlotPath', matchedGainPlotPath);
    inventory.matchedGainFits = matchedGainFits;
    inventory.matchedGainPlotPath = matchedGainPlotPath;

    matchedGainTablePath = fullfile( ...
      summaryFolder, 'IDR_MatchedAbsoluteGains_ByOffset.csv');
    writetable(matchedGainFits.summaryTable, matchedGainTablePath);
    inventory.matchedGainTablePath = matchedGainTablePath;

    fprintf('\nIDR offset-specific four-gain fits: rectangular step\n');
    displayVars = { ...
      'probeDirDeg','nTrials','nSessions','fullNLL', ...
      'deltaNLLProbeTerms','lrtPProbeTerms', ...
      'gCP','seCP','gCQ','seCQ','gNP','seNP','gNQ','seNQ'};
    disp(matchedGainFits.summaryTable(:,displayVars));
    fprintf('Saved matched-gain table: %s\n',matchedGainTablePath);
    fprintf('Saved matched-gain plot: %s\n',matchedGainPlotPath);
  end

  if opts.FitSignedGains
    signedGainPlotPath = fullfile( ...
      plotFolder, 'IDR_SignedNoChangeGains_ByOffset.pdf');
    signedGainFits = analyzeIDRSignedGains( ...
      sideGainData, finalBeta, finalLapse, opts.TargetPerformance, ...
      opts.GainBounds, 'PreferredNullRatio', opts.PreferredNullRatio, ...
      'PlotPath', signedGainPlotPath);
    inventory.signedGainFits = signedGainFits;
    inventory.signedGainPlotPath = signedGainPlotPath;

    signedGainTablePath = fullfile( ...
      summaryFolder, 'IDR_SignedGains_ByOffset.csv');
    writetable(signedGainFits.summaryTable, signedGainTablePath);
    inventory.signedGainTablePath = signedGainTablePath;

    fprintf('\nIDR offset-specific signed-gain fits: rectangular step\n');
    displayVars = { ...
      'probeDirDeg','nTrials','deltaNLLNoChangeSigned','pNoChangeSigned', ...
      'gNPpos','seNPpos','gNPneg','seNPneg', ...
      'gNQpos','seNQpos','gNQneg','seNQneg', ...
      'deltaNLLOpponentVsFree','pOpponentVsFree'};
    disp(signedGainFits.summaryTable(:,displayVars));
    fprintf('Saved signed-gain table: %s\n',signedGainTablePath);
    fprintf('Saved signed-gain plot: %s\n',signedGainPlotPath);
  end

  if opts.FitControlledPooling
    if ~opts.FitSignedGains
      error('analyzeIDRSideGains:SignedGainsRequired', ...
        'FitControlledPooling requires FitSignedGains=true.');
    end
    poolingPlotPath = fullfile( ...
      plotFolder, 'IDR_ControlledPooling_ByOffset.pdf');
    controlledPoolingFits = analyzeIDRControlledPooling( ...
      sideGainData, signedGainFits, finalBeta, finalLapse, ...
      opts.TargetPerformance, opts.GainBounds, ...
      'PreferredNullRatio', opts.PreferredNullRatio, ...
      'PBounds', opts.PoolingPBounds, ...
      'ProfilePValues', opts.PoolingProfilePValues, ...
      'PlotPath', poolingPlotPath);
    inventory.controlledPoolingFits = controlledPoolingFits;
    inventory.controlledPoolingPlotPath = poolingPlotPath;

    poolingOffsetTablePath = fullfile( ...
      summaryFolder, 'IDR_ControlledPooling_ByOffset.csv');
    poolingSessionTablePath = fullfile( ...
      summaryFolder, 'IDR_ControlledPooling_BySession.csv');
    poolingOverallTablePath = fullfile( ...
      summaryFolder, 'IDR_ControlledPooling_Overall.csv');
    poolingProfileTablePath = fullfile( ...
      summaryFolder, 'IDR_ControlledPooling_SharedPProfile.csv');
    writetable(controlledPoolingFits.offsetSummary,poolingOffsetTablePath);
    writetable(controlledPoolingFits.sessionSummary,poolingSessionTablePath);
    writetable(controlledPoolingFits.overallSummary,poolingOverallTablePath);
    writetable(controlledPoolingFits.sharedPNorm.profile.table, ...
      poolingProfileTablePath);
    inventory.controlledPoolingTablePaths = string([ ...
      poolingOffsetTablePath;poolingSessionTablePath;poolingOverallTablePath; ...
      poolingProfileTablePath]);

    fprintf('\nIDR controlled opponent-rectified no-change pooling\n');
    disp(controlledPoolingFits.offsetSummary(:,{ ...
      'probeDirDeg','nTrials','sumNLL','pNormNLL','hardMaxNLL', ...
      'sharedPNormNLL','deltaNLLSumMinusMax', ...
      'deltaNLLSumMinusShared','fittedP','maxProbeWinnerFraction'}));
    fprintf('\nIDR controlled pooling overall summary\n');
    disp(controlledPoolingFits.overallSummary);
    fprintf('Shared p-norm fit: p %.5g, NLL %.3f, exit flag %g\n', ...
      controlledPoolingFits.sharedPNorm.p, ...
      controlledPoolingFits.sharedPNorm.negLogLikelihood, ...
      controlledPoolingFits.sharedPNorm.exitflag);
    fprintf('\nIDR shared-p profile likelihood\n');
    disp(controlledPoolingFits.sharedPNorm.profile.table(:,1:6));
    fprintf(['Profile lower bounds on p: one-sided 95%% %.4g; ' ...
      'two-sided 95%% %.4g\n'], ...
      controlledPoolingFits.sharedPNorm.profile.lowerOneSided95, ...
      controlledPoolingFits.sharedPNorm.profile.lowerTwoSided95);
    fprintf('Profile endpoint checks: p=1 minus sum NLL %.6g; p=max minus hard-max NLL %.6g\n', ...
      controlledPoolingFits.sharedPNorm.profile.sumEndpointDifference, ...
      controlledPoolingFits.sharedPNorm.profile.maxEndpointDifference);
    fprintf('Saved pooling tables: %s; %s; %s; %s\n', ...
      poolingOffsetTablePath,poolingSessionTablePath,poolingOverallTablePath, ...
      poolingProfileTablePath);
    fprintf('Saved pooling plot: %s\n',poolingPlotPath);
  end

  printIDRSideGainSummary(sideGainFits);

  gainPlotPath = fullfile(plotFolder, 'IDR_SideGainModelComparison.pdf');
  plotIDRSideGainSummary(sideGainFits, gainPlotPath);
  inventory.sideGainPlotPath = gainPlotPath;
  fprintf('Saved side-gain comparison: %s\n', gainPlotPath);
end

if opts.PrintInventory
  printTrialInventory(sessionSummary, conditionSummary);
end

if opts.SaveSummary
  summaryPath = fullfile(summaryFolder, 'IDR_SideGainAnalysis.mat');
  save(summaryPath, 'inventory', '-v7.3');
  inventory.summaryPath = summaryPath;
  fprintf('\nSaved inventory: %s\n', summaryPath);
end
end

%% ------------------------------------------------------------------------
function D = assembleIDRSideGainTrials( ...
  domainPath, parentTrialTable, sessionPsychFits, animal, excludedProbeDirectionsDeg)
% Assemble disjoint probeSession trials into one INC noise-trial table.
% Both preferred and probe noise are retained. Existing changeNoise and
% noChangeNoise outputs remain aliases for preferred noise so established
% side-gain results are unchanged.

dataPath = fullfile(domainPath, 'Data');
probeDirs = dir(fullfile(dataPath, 'Probe*'));
probeDirs = probeDirs([probeDirs.isdir]);

probeSessionFolders = {};
for iDir = 1:numel(probeDirs)
  candidate = fullfile(probeDirs(iDir).folder, probeDirs(iDir).name, 'ProbeSessions');
  if exist(candidate, 'dir')
    probeSessionFolders{end+1} = candidate; %#ok<AGROW>
  end
end
probeSessionFolders = unique(probeSessionFolders, 'stable');

[filePaths, fileInfo] = selectAnalysisFiles( ...
  probeSessionFolders, ...
  'Animal', animal, ...
  'Bin180Into179', false);

if isempty(filePaths)
  error('analyzeIDRSideGains:NoProbeSessions', ...
    'No probeSession files were selected.');
end

[filePaths, fileInfo, exclusionSummary] = excludeIDRProbeFiles( ...
  filePaths, fileInfo, excludedProbeDirectionsDeg);
if isempty(filePaths)
  error('analyzeIDRSideGains:NoProbeSessionsAfterExclusion', ...
    'No ProbeSession files remain after excluding probe directions: %s.', ...
    mat2str(excludedProbeDirectionsDeg));
end

trialTables = cell(numel(filePaths), 1);
changePrefNoiseCells = cell(numel(filePaths), 1);
noChangePrefNoiseCells = cell(numel(filePaths), 1);
changeProbeEffectiveNoiseCells = cell(numel(filePaths), 1);
noChangeProbeEffectiveNoiseCells = cell(numel(filePaths), 1);
changeProbeCandidateNoiseCells = cell(numel(filePaths), 1);
noChangeProbeCandidateNoiseCells = cell(numel(filePaths), 1);
tMS = [];
stepFrames = [];

for iFile = 1:numel(filePaths)
  filePath = filePaths{iFile};
  S = load(filePath, ...
    'sessionHeader', 'sessionProbeHeader', 'trialIdx', ...
    'prefNoiseByPatch', 'probeNoiseByPatch', 'trialOutcomesAll', ...
    'changeSidesAll', 'changeIndicesAll');

  required = {'sessionHeader','sessionProbeHeader','trialIdx', ...
    'prefNoiseByPatch','probeNoiseByPatch','trialOutcomesAll', ...
    'changeSidesAll','changeIndicesAll'};
  for iReq = 1:numel(required)
    if ~isfield(S, required{iReq})
      error('analyzeIDRSideGains:MissingProbeVariable', ...
        '%s is missing %s.', filePath, required{iReq});
    end
  end

  nTrials = size(S.prefNoiseByPatch, 3);
  if ~isequal(size(S.probeNoiseByPatch), size(S.prefNoiseByPatch))
    error('analyzeIDRSideGains:ProbeNoiseSizeMismatch', ...
      'Preferred and probe noise matrices differ in size in %s.', filePath);
  end
  if numel(S.trialIdx) ~= nTrials || ...
      numel(S.trialOutcomesAll) ~= nTrials || ...
      numel(S.changeSidesAll) ~= nTrials || ...
      numel(S.changeIndicesAll) ~= nTrials
    error('analyzeIDRSideGains:ProbeTrialAlignment', ...
      'ProbeSession arrays are not trial-aligned in %s.', filePath);
  end

  % Use the derived probeSession filename as the authoritative analysis
  % identity. sessionHeader.fileName can retain a typo-laden acquisition name
  % that differs from the converted FullSession filename.
  [~, probeSessionBaseName] = fileparts(filePath);

  % For example:
  %   IDReadout_Meetz_20260410_Probe180
  % becomes:
  %   IDReadout_Meetz_20260410
  parentID = regexprep(string(probeSessionBaseName), ...
    '_Probe[^_]*$', '', 'ignorecase');

  sessionFitIDs = erase(string(sessionPsychFits.sessionID), ".mat");
  parentRow = find(sessionFitIDs == parentID);

  if numel(parentRow) ~= 1
    error('analyzeIDRSideGains:ParentSessionMatch', ...
      ['Expected one parent session fit for analysis name %s, derived from ' ...
      'probeSession file %s; found %d. Header provenance name was %s.'], ...
      parentID, probeSessionBaseName, numel(parentRow), ...
      string(S.sessionHeader.fileName));
  end

  % parentID = string(S.sessionHeader.fileName);
  % parentRow = find(string(sessionPsychFits.sessionID) == parentID + ".mat" | ...
  %                  erase(string(sessionPsychFits.sessionID), ".mat") == parentID);
  % if numel(parentRow) ~= 1
  %   error('analyzeIDRSideGains:ParentSessionMatch', ...
  %     'Expected one parent session fit for %s; found %d.', parentID, numel(parentRow));
  % end
  parentSessionIndex = sessionPsychFits.sessionIndex(parentRow);
  parentThreshold = sessionPsychFits.threshold(parentRow);

  use = S.changeIndicesAll(:) == 2;
  localTrialIdx = double(S.trialIdx(:));
  localOutcome = double(S.trialOutcomesAll(:));
  localChangeSide = double(S.changeSidesAll(:));

  trialIdxUse = localTrialIdx(use);
  outcomeUse = localOutcome(use);
  changeSideUse = localChangeSide(use);
  nUse = sum(use);

  % Join to the compact parent behavior table using sessionIndex + trialIdx.
  parentRows = parentTrialTable.sessionIndex == parentSessionIndex;
  parentSubset = parentTrialTable(parentRows, :);
  [tf, loc] = ismember(trialIdxUse, double(parentSubset.trialIdx));
  if ~all(tf)
    error('analyzeIDRSideGains:ParentTrialJoin', ...
      '%d probe trials in %s were not found in the parent behavior table.', ...
      sum(~tf), filePath);
  end
  P = parentSubset(loc, :);

  if any(~P.hasStepNoise) || any(P.changeIndex ~= 2)
    error('analyzeIDRSideGains:ParentTrialMismatch', ...
      'Joined parent rows do not describe INC noise trials in %s.', filePath);
  end
  if any(double(P.correct) ~= double(outcomeUse == 0))
    error('analyzeIDRSideGains:OutcomeMismatch', ...
      'ProbeSession and parent outcomes disagree in %s.', filePath);
  end

  prefNoise = S.prefNoiseByPatch(:, :, use);
  probeEffectiveNoise = S.probeNoiseByPatch(:, :, use);
  probeDirDeg = double(S.sessionProbeHeader.probeDirDeg);
  nYokedProbeStreams = idrProbeStreamCount(probeDirDeg);
  probeCandidateNoise = probeEffectiveNoise ./ nYokedProbeStreams;
  nFrames = size(prefNoise, 2);
  changePrefNoise = nan(nFrames, nUse);
  noChangePrefNoise = nan(nFrames, nUse);
  changeProbeEffectiveNoise = nan(nFrames, nUse);
  noChangeProbeEffectiveNoise = nan(nFrames, nUse);
  changeProbeCandidateNoise = nan(nFrames, nUse);
  noChangeProbeCandidateNoise = nan(nFrames, nUse);

  changePatch = changeSideUse + 1;
  noChangePatch = 3 - changePatch;
  for iTrial = 1:nUse
    changePrefNoise(:, iTrial) = squeeze( ...
      prefNoise(changePatch(iTrial), :, iTrial)).';
    noChangePrefNoise(:, iTrial) = squeeze( ...
      prefNoise(noChangePatch(iTrial), :, iTrial)).';
    changeProbeEffectiveNoise(:, iTrial) = squeeze( ...
      probeEffectiveNoise(changePatch(iTrial), :, iTrial)).';
    noChangeProbeEffectiveNoise(:, iTrial) = squeeze( ...
      probeEffectiveNoise(noChangePatch(iTrial), :, iTrial)).';
    changeProbeCandidateNoise(:, iTrial) = squeeze( ...
      probeCandidateNoise(changePatch(iTrial), :, iTrial)).';
    noChangeProbeCandidateNoise(:, iTrial) = squeeze( ...
      probeCandidateNoise(noChangePatch(iTrial), :, iTrial)).';
  end

  frameRateHz = double(S.sessionHeader.frameRateHz);
  preStepMS = double(S.sessionHeader.preStepMS);
  stepMS = double(S.sessionHeader.stepMS);
  thisTMS = (0:nFrames-1)' .* (1000 / frameRateHz);
  firstStepFrame = round(preStepMS / (1000 / frameRateHz)) + 1;
  nStepFrames = round(stepMS / (1000 / frameRateHz));
  thisStepFrames = (firstStepFrame:firstStepFrame+nStepFrames-1)';

  if isempty(tMS)
    tMS = thisTMS;
    stepFrames = thisStepFrames;
  elseif numel(tMS) ~= numel(thisTMS) || any(abs(tMS - thisTMS) > 1e-9) || ...
      numel(stepFrames) ~= numel(thisStepFrames) || any(stepFrames ~= thisStepFrames)
    error('analyzeIDRSideGains:ProbeTimingMismatch', ...
      'ProbeSession timing differs in %s.', filePath);
  end

  T = table();
  T.sessionIndex = repmat(parentSessionIndex, nUse, 1);
  T.sessionID = repmat(parentID, nUse, 1);
  T.animal = repmat(string(S.sessionHeader.animal), nUse, 1);
  T.probeSessionIndex = repmat(iFile, nUse, 1);
  T.probeDirDeg = repmat(probeDirDeg, nUse, 1);
  T.nYokedProbeStreams = repmat(nYokedProbeStreams, nUse, 1);
  T.trialIdx = trialIdxUse;
  T.stepCohPC = double(P.stepCohPC);
  T.stepDeltaCohPC = double(P.stepCohPC - P.baseCohPC);
  T.prefCohNoisePC = double(P.prefCohNoisePC);
  T.probeCohNoisePC = double(P.probeCohNoisePC);
  T.correct = double(P.correct);
  T.changeSide = changeSideUse;
  T.sessionThreshold = repmat(parentThreshold, nUse, 1);
  T.sourceFile = repmat(string(filePath), nUse, 1);

  trialTables{iFile} = T;
  changePrefNoiseCells{iFile} = changePrefNoise;
  noChangePrefNoiseCells{iFile} = noChangePrefNoise;
  changeProbeEffectiveNoiseCells{iFile} = changeProbeEffectiveNoise;
  noChangeProbeEffectiveNoiseCells{iFile} = noChangeProbeEffectiveNoise;
  changeProbeCandidateNoiseCells{iFile} = changeProbeCandidateNoise;
  noChangeProbeCandidateNoiseCells{iFile} = noChangeProbeCandidateNoise;
end

trialTable = vertcat(trialTables{:});
changePrefNoiseByFrameTrial = horzcat(changePrefNoiseCells{:});
noChangePrefNoiseByFrameTrial = horzcat(noChangePrefNoiseCells{:});
changeProbeEffectiveNoiseByFrameTrial = horzcat(changeProbeEffectiveNoiseCells{:});
noChangeProbeEffectiveNoiseByFrameTrial = horzcat(noChangeProbeEffectiveNoiseCells{:});
changeProbeCandidateNoiseByFrameTrial = horzcat(changeProbeCandidateNoiseCells{:});
noChangeProbeCandidateNoiseByFrameTrial = horzcat(noChangeProbeCandidateNoiseCells{:});

if height(trialTable) ~= size(changePrefNoiseByFrameTrial, 2)
  error('analyzeIDRSideGains:ConcatenationMismatch', ...
    'Trial table and noise matrices have inconsistent trial counts.');
end

% Trials must be unique across probeSessions within each parent session.
keys = string(trialTable.sessionID) + ":" + string(trialTable.trialIdx);
if numel(unique(keys)) ~= numel(keys)
  error('analyzeIDRSideGains:DuplicateProbeTrials', ...
    'Duplicate parent session/trialIdx rows were found across probeSessions.');
end

trialAccounting = makeIDRSideGainTrialAccounting( ...
  parentTrialTable, trialTable, excludedProbeDirectionsDeg);
nMissing = sum(trialAccounting.nMissingFromProbeSessions);
nExtra = sum(trialAccounting.nExtraInProbeSessions);
fprintf('  Parent/probe trial accounting: %d missing, %d extra\n', ...
  nMissing, nExtra);
if nMissing > 0 || nExtra > 0
  disp(trialAccounting(trialAccounting.difference ~= 0, :));
end

D = struct();
D.trialTable = trialTable;
D.changePrefNoiseByFrameTrial = changePrefNoiseByFrameTrial;
D.noChangePrefNoiseByFrameTrial = noChangePrefNoiseByFrameTrial;
D.changeProbeEffectiveNoiseByFrameTrial = changeProbeEffectiveNoiseByFrameTrial;
D.noChangeProbeEffectiveNoiseByFrameTrial = noChangeProbeEffectiveNoiseByFrameTrial;
D.changeProbeCandidateNoiseByFrameTrial = changeProbeCandidateNoiseByFrameTrial;
D.noChangeProbeCandidateNoiseByFrameTrial = noChangeProbeCandidateNoiseByFrameTrial;
% Backward-compatible aliases used by the established preferred-noise fits.
D.changeNoiseByFrameTrial = changePrefNoiseByFrameTrial;
D.noChangeNoiseByFrameTrial = noChangePrefNoiseByFrameTrial;
D.tMS = tMS;
D.stepFrames = stepFrames;
D.probeFileInfo = fileInfo;
D.animalSelection = string(animal);
D.trialAccounting = trialAccounting;
D.excludedProbeDirectionsDeg = excludedProbeDirectionsDeg(:)';
D.probeFileExclusionSummary = exclusionSummary;
D.probeScalingContract = [ ...
  'probeEffective is the sum over yoked probe streams; ' ...
  'probeCandidate is the physical predictor for one candidate mechanism'];
D.commonKernelDefinition = ...
  'correct-minus-error kernel of changeNoise-minus-noChangeNoise';
end

%% ------------------------------------------------------------------------
function P = computeIDRSideGainPredictors(D)
% Construct rectangular and leave-one-parent-session-out common-weight
% predictors. The same temporal weights are applied to change and no-change
% noise streams.

T = D.trialTable;
changeNoise = D.changeNoiseByFrameTrial;
noChangeNoise = D.noChangeNoiseByFrameTrial;
nFrames = size(changeNoise, 1);
allFrames = (1:nFrames)';
stepFrames = D.stepFrames(:);

types = {'rectStep','rectFull','looKernelStep','looKernelFull'};
P = struct();
P.predictorTypes = types;
P.tMS = D.tMS;
P.stepFrames = stepFrames;
P.commonKernelDefinition = D.commonKernelDefinition;

% Rectangular predictors are means in coherence units.
P.rectStep = makePredictorRecord( ...
  mean(changeNoise(stepFrames, :), 1, 'omitnan')', ...
  mean(noChangeNoise(stepFrames, :), 1, 'omitnan')', ...
  rectangularWeights(nFrames, stepFrames), ...
  'mean noise over step');
P.rectFull = makePredictorRecord( ...
  mean(changeNoise(allFrames, :), 1, 'omitnan')', ...
  mean(noChangeNoise(allFrames, :), 1, 'omitnan')', ...
  rectangularWeights(nFrames, allFrames), ...
  'mean noise over full stimulus');

sessionIDs = unique(T.sessionIndex, 'stable');
nSessions = numel(sessionIDs);
weightsStep = nan(nFrames, nSessions);
weightsFull = nan(nFrames, nSessions);
changePredStep = nan(height(T), 1);
noChangePredStep = nan(height(T), 1);
changePredFull = nan(height(T), 1);
noChangePredFull = nan(height(T), 1);

diagnosticSessionIndex = nan(nSessions, 1);
kernelStepSum = nan(nSessions, 1);
kernelFullSum = nan(nSessions, 1);
effectiveNStep = nan(nSessions, 1);
effectiveNFull = nan(nSessions, 1);

diffNoise = changeNoise - noChangeNoise;

for iSession = 1:nSessions
  sessionIndex = sessionIDs(iSession);
  leaveOut = T.sessionIndex == sessionIndex;
  train = ~leaveOut;

  correctTrain = train & T.correct == 1;
  errorTrain = train & T.correct == 0;

  if ~any(correctTrain) || ~any(errorTrain)
    error('analyzeIDRSideGains:LOOOutcomeClass', ...
      'LOO training data lack an outcome class for session %d.', sessionIndex);
  end

  kernel = mean(diffNoise(:, correctTrain), 2, 'omitnan') - ...
           mean(diffNoise(:, errorTrain), 2, 'omitnan');

  wStep = normalizedKernelWeights(kernel, stepFrames);
  wFull = normalizedKernelWeights(kernel, allFrames);
  weightsStep(:, iSession) = wStep;
  weightsFull(:, iSession) = wFull;

  changePredStep(leaveOut) = sum( ...
    changeNoise(:, leaveOut) .* wStep, 1, 'omitnan')';
  noChangePredStep(leaveOut) = sum( ...
    noChangeNoise(:, leaveOut) .* wStep, 1, 'omitnan')';
  changePredFull(leaveOut) = sum( ...
    changeNoise(:, leaveOut) .* wFull, 1, 'omitnan')';
  noChangePredFull(leaveOut) = sum( ...
    noChangeNoise(:, leaveOut) .* wFull, 1, 'omitnan')';

  diagnosticSessionIndex(iSession) = sessionIndex;
  kernelStepSum(iSession) = sum(kernel(stepFrames), 'omitnan');
  kernelFullSum(iSession) = sum(kernel(allFrames), 'omitnan');
  effectiveNStep(iSession) = 1 / sum(wStep.^2, 'omitnan');
  effectiveNFull(iSession) = 1 / sum(wFull.^2, 'omitnan');
end

P.looKernelStep = makePredictorRecord( ...
  changePredStep, noChangePredStep, weightsStep, ...
  'LOO common diff-kernel weights over step');
P.looKernelFull = makePredictorRecord( ...
  changePredFull, noChangePredFull, weightsFull, ...
  'LOO common diff-kernel weights over full stimulus');

P.weightDiagnostics = table( ...
  diagnosticSessionIndex, kernelStepSum, kernelFullSum, ...
  effectiveNStep, effectiveNFull, ...
  'VariableNames', {'sessionIndex','kernelStepSum','kernelFullSum', ...
  'effectiveNStep','effectiveNFull'});

% Pooled common kernel for display only; trial predictors remain LOO.
correctAll = T.correct == 1;
errorAll = T.correct == 0;
P.pooledCommonKernel = ...
  mean(diffNoise(:, correctAll), 2, 'omitnan') - ...
  mean(diffNoise(:, errorAll), 2, 'omitnan');
P.pooledCommonWeightsStep = normalizedKernelWeights( ...
  P.pooledCommonKernel, stepFrames);
P.pooledCommonWeightsFull = normalizedKernelWeights( ...
  P.pooledCommonKernel, allFrames);
P.predictorDiagnostics = makeIDRSideGainPredictorDiagnostics(P);
end

%% ------------------------------------------------------------------------
function A = makeIDRSideGainTrialAccounting( ...
  parentTrialTable, probeTrialTable, excludedProbeDirectionsDeg)

sessionIndex = unique(parentTrialTable.sessionIndex, 'stable');
nSessions = numel(sessionIndex);
sessionID = strings(nSessions,1);
nParentNoiseINC = zeros(nSessions,1);
nProbeNoiseINC = zeros(nSessions,1);
nMissingFromProbeSessions = zeros(nSessions,1);
nExtraInProbeSessions = zeros(nSessions,1);
nExcludedProbeTrials = zeros(nSessions,1);
for i = 1:nSessions
  idx = sessionIndex(i);

  allParentNoiseINC = parentTrialTable( ...
    parentTrialTable.sessionIndex == idx & ...
    parentTrialTable.isValid & ...
    parentTrialTable.changeIndex == 2 & ...
    parentTrialTable.hasStepNoise, :);

  nExcludedProbeTrials(i) = sum(ismember( ...
    allParentNoiseINC.probeDirDeg, excludedProbeDirectionsDeg));

  parent = allParentNoiseINC(~ismember( ...
    allParentNoiseINC.probeDirDeg, excludedProbeDirectionsDeg), :);
  probe = probeTrialTable(probeTrialTable.sessionIndex == idx, :);

  if ~isempty(parent)
    sessionID(i) = string(parent.sessionID(1));
  elseif ~isempty(probe)
    sessionID(i) = string(probe.sessionID(1));
  end

  parentKeys = unique(double(parent.trialIdx));
  probeKeys = unique(double(probe.trialIdx));
  nParentNoiseINC(i) = numel(parentKeys);
  nProbeNoiseINC(i) = numel(probeKeys);
  nMissingFromProbeSessions(i) = numel(setdiff(parentKeys, probeKeys));
  nExtraInProbeSessions(i) = numel(setdiff(probeKeys, parentKeys));
end

difference = nProbeNoiseINC - nParentNoiseINC;
A = table(sessionIndex, sessionID, nParentNoiseINC, nProbeNoiseINC, ...
  difference, nMissingFromProbeSessions, nExtraInProbeSessions, ...
  nExcludedProbeTrials);
end

%% ------------------------------------------------------------------------
function [filePaths, fileInfo, summary] = excludeIDRProbeFiles( ...
  filePaths, fileInfo, excludedProbeDirectionsDeg)

n = numel(filePaths);
probeDirDeg = nan(n,1);
excluded = false(n,1);
for i = 1:n
  S = load(filePaths{i}, 'sessionProbeHeader');
  if ~isfield(S, 'sessionProbeHeader') || ...
      ~isfield(S.sessionProbeHeader, 'probeDirDeg')
    error('analyzeIDRSideGains:MissingProbeDirection', ...
      'Missing sessionProbeHeader.probeDirDeg in %s.', filePaths{i});
  end
  probeDirDeg(i) = double(S.sessionProbeHeader.probeDirDeg);
  excluded(i) = ismember(probeDirDeg(i), excludedProbeDirectionsDeg);
end

summary = table(string(filePaths(:)), probeDirDeg, excluded, ...
  'VariableNames', {'filePath','probeDirDeg','excluded'});
if any(excluded)
  fprintf('  Explicitly excluded %d ProbeSession file(s) at offset(s) %s deg.\n', ...
    sum(excluded), mat2str(unique(probeDirDeg(excluded))'));
end
filePaths = filePaths(~excluded);
fileInfo = fileInfo(~excluded,:);
end

%% ------------------------------------------------------------------------
function n = idrProbeStreamCount(probeDirDeg)

probeDirDeg = abs(double(probeDirDeg));
if probeDirDeg > 0 && probeDirDeg < 180
  n = 2;
elseif abs(probeDirDeg - 180) < 1e-9
  n = 1;
else
  error('analyzeIDRSideGains:UnsupportedProbeDirection', ...
    'Cannot determine the number of probe candidates at offset %g deg.', ...
    probeDirDeg);
end
end

%% ------------------------------------------------------------------------
function D = makeIDRSideGainPredictorDiagnostics(P)

types = string(P.predictorTypes(:));
n = numel(types);
meanChange = nan(n,1);
sdChange = nan(n,1);
meanNoChange = nan(n,1);
sdNoChange = nan(n,1);
changeNoChangeCorrelation = nan(n,1);

for i = 1:n
  R = P.(char(types(i)));
  x = R.change;
  y = R.noChange;
  valid = isfinite(x) & isfinite(y);
  meanChange(i) = mean(x(valid));
  sdChange(i) = std(x(valid));
  meanNoChange(i) = mean(y(valid));
  sdNoChange(i) = std(y(valid));
  C = corrcoef(x(valid), y(valid));
  if numel(C) == 4
    changeNoChangeCorrelation(i) = C(1,2);
  end
end

D = table(types, meanChange, sdChange, meanNoChange, sdNoChange, ...
  changeNoChangeCorrelation, ...
  'VariableNames', {'predictorType','meanChange','sdChange', ...
  'meanNoChange','sdNoChange','changeNoChangeCorrelation'});
end

%% ------------------------------------------------------------------------
function plotIDRSideGainCommonKernel(P, pdfPath)

tMS = P.tMS(:);
kernel = P.pooledCommonKernel(:);
wStep = P.pooledCommonWeightsStep(:);
wFull = P.pooledCommonWeightsFull(:);
stepFrames = P.stepFrames(:);

fig = figure('Color','w','Visible','off', ...
  'Units','inches','Position',[1 1 10 7.5]);
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
title(tl, ['IDR common temporal weighting: ' ...
  'change-minus-no-change preferred noise']);

ax = nexttile(tl);
plot(ax,tMS,kernel,'LineWidth',1.4);
yline(ax,0,'k:');
xline(ax,tMS(stepFrames(1)),'k--');
xline(ax,tMS(stepFrames(end)),'k--');
ylabel(ax,'Kernel');
title(ax,'Pooled correct-minus-error common kernel');
grid(ax,'on');
box(ax,'off');

ax = nexttile(tl);
plot(ax,tMS,wStep,'LineWidth',1.4);
yline(ax,0,'k:');
xline(ax,tMS(stepFrames(1)),'k--');
xline(ax,tMS(stepFrames(end)),'k--');
ylabel(ax,'Weight');
title(ax,sprintf('Step-normalized weights; effective N = %.1f', ...
  1/sum(wStep.^2,'omitnan')));
grid(ax,'on');
box(ax,'off');

ax = nexttile(tl);
plot(ax,tMS,wFull,'LineWidth',1.4);
yline(ax,0,'k:');
xline(ax,tMS(stepFrames(1)),'k--');
xline(ax,tMS(stepFrames(end)),'k--');
xlabel(ax,'Time from stimulus onset (ms)');
ylabel(ax,'Weight');
title(ax,sprintf('Full-interval normalized weights; effective N = %.1f', ...
  1/sum(wFull.^2,'omitnan')));
grid(ax,'on');
box(ax,'off');

exportgraphics(fig,pdfPath,'ContentType','vector');
close(fig);
end

%% ------------------------------------------------------------------------
function R = makePredictorRecord(changePredictor, noChangePredictor, weights, description)
R = struct();
R.change = changePredictor(:);
R.noChange = noChangePredictor(:);
R.weights = weights;
R.description = description;
end

%% ------------------------------------------------------------------------
function w = rectangularWeights(nFrames, useFrames)
w = zeros(nFrames, 1);
w(useFrames) = 1 / numel(useFrames);
end

%% ------------------------------------------------------------------------
function w = normalizedKernelWeights(kernel, useFrames)
w = zeros(size(kernel));
k = kernel(useFrames);
kSum = sum(k, 'omitnan');
if ~isfinite(kSum) || abs(kSum) < 1e-12
  error('analyzeIDRSideGains:BadKernelSum', ...
    'Common LOO kernel sum is zero or nonfinite.');
end
w(useFrames) = k ./ kSum;
end

%% ------------------------------------------------------------------------
function F = fitIDRSideGainModels(T, P, betaWeibull, lapse, targetPerformance, gainBounds)

types = P.predictorTypes;
F = struct();
F.betaWeibull = betaWeibull;
F.lapse = lapse;
F.targetPerformance = targetPerformance;
F.gainBounds = gainBounds;
F.predictorTypes = types;

for iType = 1:numel(types)
  type = types{iType};
  R = P.(type);
  X = [R.change, R.noChange];

  common = struct();
  common.stepCoh = double(T.stepDeltaCohPC);
  common.correct = double(T.correct);
  common.sessionThreshold = double(T.sessionThreshold);
  common.betaWeibull = betaWeibull;
  common.lapse = lapse;
  common.targetPerformance = targetPerformance;

  unrestricted = fitOneIDRSideGainModel( ...
    type + "_twoGain", ["Change","No-change"], X, [1; 0], common, gainBounds);

  equalX = sum(X, 2);
  equalGain = fitOneIDRSideGainModel( ...
    type + "_equalGain", "Common", equalX, 0.5, common, gainBounds);

  changeOnly = fitOneIDRSideGainModel( ...
    type + "_changeOnly", "Change", X(:,1), 1, common, gainBounds);

  contrast = [1 -1];
  deltaGain = contrast * unrestricted.gain;
  deltaVar = contrast * unrestricted.covariance * contrast';
  if isfinite(deltaVar) && deltaVar > 0
    deltaSE = sqrt(deltaVar);
    deltaCI95 = deltaGain + [-1 1] * 1.96 * deltaSE;
    deltaZ = deltaGain / deltaSE;
  else
    deltaSE = nan;
    deltaCI95 = [nan nan];
    deltaZ = nan;
  end

  result = struct();
  result.predictorType = type;
  result.predictorDescription = R.description;
  result.unrestricted = unrestricted;
  result.equalGain = equalGain;
  result.changeOnly = changeOnly;
  result.deltaGain = deltaGain;
  result.deltaGainSE = deltaSE;
  result.deltaGainCI95 = deltaCI95;
  result.deltaGainZ = deltaZ;
  result.deltaNLLSeparateVsEqual = ...
    equalGain.negLogLikelihood - unrestricted.negLogLikelihood;
  result.deltaNLLChangeOnlyVsTwoGain = ...
    changeOnly.negLogLikelihood - unrestricted.negLogLikelihood;
  F.(type) = result;
end

nTypes = numel(types);
predictorType = strings(nTypes,1);
nll = nan(nTypes,1);
gChange = nan(nTypes,1);
gNoChange = nan(nTypes,1);
seChange = nan(nTypes,1);
seNoChange = nan(nTypes,1);
deltaGain = nan(nTypes,1);
deltaGainSE = nan(nTypes,1);
deltaNLLSeparateVsEqual = nan(nTypes,1);
deltaNLLChangeOnlyVsTwoGain = nan(nTypes,1);

for iType = 1:nTypes
  type = types{iType};
  R = F.(type);
  predictorType(iType) = type;
  nll(iType) = R.unrestricted.negLogLikelihood;
  gChange(iType) = R.unrestricted.gain(1);
  gNoChange(iType) = R.unrestricted.gain(2);
  seChange(iType) = R.unrestricted.SE(1);
  seNoChange(iType) = R.unrestricted.SE(2);
  deltaGain(iType) = R.deltaGain;
  deltaGainSE(iType) = R.deltaGainSE;
  deltaNLLSeparateVsEqual(iType) = R.deltaNLLSeparateVsEqual;
  deltaNLLChangeOnlyVsTwoGain(iType) = R.deltaNLLChangeOnlyVsTwoGain;
end

deltaNLLFromBest = nll - min(nll);
F.comparisonTable = table( ...
  predictorType, nll, deltaNLLFromBest, ...
  gChange, seChange, gNoChange, seNoChange, ...
  deltaGain, deltaGainSE, ...
  deltaNLLSeparateVsEqual, deltaNLLChangeOnlyVsTwoGain);
end

%% ------------------------------------------------------------------------
function fit = fitOneIDRSideGainModel( ...
  modelName, predictorNames, X, gain0, common, gainBounds)

X = double(X);
valid = all(isfinite(X),2) & isfinite(common.stepCoh) & ...
  isfinite(common.correct) & isfinite(common.sessionThreshold);

X = X(valid,:);
stepCoh = common.stepCoh(valid);
correct = common.correct(valid);
sessionThreshold = common.sessionThreshold(valid);

nParameters = size(X,2);
lb = repmat(gainBounds(1), nParameters, 1);
ub = repmat(gainBounds(2), nParameters, 1);
gain0 = gain0(:);

objective = @(g) idrSideGainNLL( ...
  g, X, stepCoh, correct, sessionThreshold, ...
  common.betaWeibull, common.lapse, common.targetPerformance);

opts = optimoptions('fmincon', ...
  'Display','off', 'Algorithm','interior-point', ...
  'OptimalityTolerance',1e-8, 'StepTolerance',1e-10, ...
  'MaxFunctionEvaluations',5000, 'MaxIterations',2000);

[gainHat, nll, exitflag] = fmincon( ...
  objective, gain0, [], [], [], [], lb, ub, [], opts);

gainHat = gainHat(:);
H = finiteDifferenceHessian(objective, gainHat);
covariance = nan(nParameters);
SE = nan(nParameters,1);
CI95 = nan(nParameters,2);

if all(isfinite(H),'all')
  H = (H + H') / 2;
  if all(eig(H) > 0)
    covariance = inv(H);
    variance = diag(covariance);
    ok = isfinite(variance) & variance > 0;
    SE(ok) = sqrt(variance(ok));
    CI95(ok,:) = gainHat(ok) + [-1 1] .* (1.96 * SE(ok));
  end
end

fit = struct();
fit.model = char(modelName);
fit.predictorNames = string(predictorNames(:))';
fit.gain = gainHat;
fit.SE = SE;
fit.CI95 = CI95;
fit.hessian = H;
fit.covariance = covariance;
fit.negLogLikelihood = nll;
fit.exitflag = exitflag;
fit.nParameters = nParameters;
fit.nTrials = numel(correct);
fit.nEffectiveCohClipped = sum(stepCoh + X*gainHat < 0);
end

%% ------------------------------------------------------------------------
function nll = idrSideGainNLL( ...
  gain, X, stepCoh, correct, sessionThreshold, ...
  betaWeibull, lapse, targetPerformance)

effectiveCoh = stepCoh + X * gain(:);
effectiveCoh = max(effectiveCoh, 0);

p = fixedShapeWeibullP( ...
  effectiveCoh, sessionThreshold, betaWeibull, lapse, targetPerformance);
p = min(max(p, eps), 1-eps);
nll = -sum(correct .* log(p) + (1-correct).*log(1-p));
end

%% ------------------------------------------------------------------------
function printIDRSideGainSummary(F)

fprintf('\nIDR change/no-change gain fits\n');
disp(F.comparisonTable);

types = F.predictorTypes;
for iType = 1:numel(types)
  type = types{iType};
  R = F.(type);
  U = R.unrestricted;
  fprintf(['  %-14s NLL %.3f: gChange %.4f +/- %.4f, ' ...
    'gNoChange %.4f +/- %.4f, delta %.4f +/- %.4f, ' ...
    'DeltaNLL(equal-separate) %.3f\n'], ...
    type, U.negLogLikelihood, ...
    U.gain(1), U.SE(1), U.gain(2), U.SE(2), ...
    R.deltaGain, R.deltaGainSE, ...
    R.deltaNLLSeparateVsEqual);
end
end

%% ------------------------------------------------------------------------
function plotIDRSideGainSummary(F, pdfPath)

types = F.predictorTypes;
nTypes = numel(types);

fig = figure('Color','w','Visible','off', ...
  'Units','inches','Position',[1 1 11 8.5]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
title(tl,'IDR change/no-change side gains','FontWeight','bold');

allCI = [];
for iType = 1:nTypes
  allCI = [allCI; F.(types{iType}).unrestricted.CI95(:)]; %#ok<AGROW>
end
allCI = allCI(isfinite(allCI));
if isempty(allCI)
  gainYLim = [-1 1];
else
  yMin = min([0; allCI]);
  yMax = max([0; allCI]);
  pad = .1 * max(yMax-yMin, 1e-6);
  gainYLim = [yMin-pad yMax+pad];
end

for iType = 1:nTypes
  ax = nexttile(tl);
  R = F.(types{iType});
  U = R.unrestricted;

  bar(ax,1:2,U.gain,'BarWidth',0.55);
  hold(ax,'on');
  for j = 1:2
    plot(ax,[j j],U.CI95(j,:),'k-','LineWidth',1.3);
  end
  plot(ax,1:2,U.gain,'ko','MarkerFaceColor','k', ...
    'MarkerSize',5,'LineStyle','none');
  yline(ax,0,'k:');
  set(ax,'XTick',1:2,'XTickLabel',{'Change','No-change'});
  ylim(ax,gainYLim);
  ylabel(ax,'Noise gain');
  title(ax,types{iType},'Interpreter','none');
  grid(ax,'on');
  box(ax,'off');

  txt = sprintf(['NLL %.2f\n' ...
    '\\DeltaNLL equal %.2f\n' ...
    '\\Delta g %.3f \\pm %.3f'], ...
    U.negLogLikelihood, R.deltaNLLSeparateVsEqual, ...
    R.deltaGain, R.deltaGainSE);
  text(ax,.98,.98,txt,'Units','normalized', ...
    'HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontSize',8);
end

exportgraphics(fig,pdfPath,'ContentType','vector');
close(fig);
end

%% ------------------------------------------------------------------------
function requireBehaviorTrialVariables(T, filePath)

requiredVars = { ...
  'trialIdx', ...
  'trialEnd', ...
  'trialCertify', ...
  'isValid', ...
  'correct', ...
  'changeIndex', ...
  'hasStepNoise', ...
  'baseCohPC', ...
  'stepCohPC', ...
  'stepDeltaCohPC', ...
  'prefCohNoisePC', ...
  'probeCohNoisePC', ...
  'probeDirDeg', ...
  'changeSide', ...
  'chosenSide' ...
  };

missingVars = setdiff(requiredVars, T.Properties.VariableNames);
if ~isempty(missingVars)
  error('analyzeIDRSideGains:MissingBehaviorVariables', ...
    'behaviorTrialTable in %s is missing: %s', ...
    filePath, strjoin(missingVars, ', '));
end
end

%% ------------------------------------------------------------------------
function conditionSummary = makeConditionSummary(T)

if isempty(T)
  conditionSummary = table();
  return;
end

sessionIndex = double(T.sessionIndex);
sessionID = string(T.sessionID);
animal = string(T.animal);
hasStepNoise = logical(T.hasStepNoise);
stepDeltaCohPC = double(T.stepDeltaCohPC);
correct = double(T.correct);

valid = isfinite(sessionIndex) & strlength(sessionID) > 0 & ...
  strlength(animal) > 0 & isfinite(stepDeltaCohPC) & isfinite(correct);

sessionIndex = sessionIndex(valid);
sessionID = sessionID(valid);
animal = animal(valid);
hasStepNoise = hasStepNoise(valid);
stepDeltaCohPC = stepDeltaCohPC(valid);
correct = correct(valid);

[G, groupSessionIndex, groupSessionID, groupAnimal, ...
  groupHasStepNoise, groupstepDeltaCohPC] = findgroups( ...
  sessionIndex, sessionID, animal, hasStepNoise, stepDeltaCohPC);

nTrials = splitapply(@numel, correct, G);
nCorrect = splitapply(@sum, correct, G);
pCorrect = nCorrect ./ nTrials;

conditionSummary = table( ...
  groupSessionIndex, groupSessionID, groupAnimal, ...
  groupHasStepNoise, groupstepDeltaCohPC, ...
  nTrials, nCorrect, pCorrect, ...
  'VariableNames', { ...
    'sessionIndex', 'sessionID', 'animal', ...
    'hasStepNoise', 'stepDeltaCohPC', ...
    'nTrials', 'nCorrect', 'pCorrect'});

conditionSummary = sortrows(conditionSummary, ...
  {'sessionIndex', 'hasStepNoise', 'stepDeltaCohPC'}, ...
  {'ascend', 'ascend', 'ascend'});
end

%% ------------------------------------------------------------------------
function sessionSummary = makeSessionSummary(T, fileInfo)

nSessions = height(fileInfo);

sessionIndex = (1:nSessions)';
sessionID = strings(nSessions, 1);
animal = strings(nSessions, 1);
filePath = strings(nSessions, 1);

nINCTrials = zeros(nSessions, 1);
nNoiseINCTrials = zeros(nSessions, 1);
nNoNoiseINCTrials = zeros(nSessions, 1);
nNoiseCohLevels = zeros(nSessions, 1);
nNoNoiseCohLevels = zeros(nSessions, 1);
minNoiseCohPC = nan(nSessions, 1);
maxNoiseCohPC = nan(nSessions, 1);
minNoNoiseCohPC = nan(nSessions, 1);
maxNoNoiseCohPC = nan(nSessions, 1);
pCorrectINC = nan(nSessions, 1);
pCorrectNoiseINC = nan(nSessions, 1);
pCorrectNoNoiseINC = nan(nSessions, 1);
brackets75NoNoise = false(nSessions, 1);

for iSession = 1:nSessions
  idxSession = T.sessionIndex == iSession;
  Ts = T(idxSession, :);

  sessionID(iSession) = string(fileInfo.fileName{iSession});
  animal(iSession) = string(fileInfo.animal{iSession});
  filePath(iSession) = string(fileInfo.filePath{iSession});

  nINCTrials(iSession) = height(Ts);
  pCorrectINC(iSession) = mean(Ts.correct, 'omitnan');

  idxNoise = Ts.hasStepNoise;
  idxNoNoise = ~Ts.hasStepNoise;

  nNoiseINCTrials(iSession) = sum(idxNoise);
  nNoNoiseINCTrials(iSession) = sum(idxNoNoise);
  pCorrectNoiseINC(iSession) = mean(Ts.correct(idxNoise), 'omitnan');
  pCorrectNoNoiseINC(iSession) = mean(Ts.correct(idxNoNoise), 'omitnan');

  noiseCoh = unique(Ts.stepDeltaCohPC(idxNoise & isfinite(Ts.stepDeltaCohPC)));
  noNoiseCoh = unique(Ts.stepDeltaCohPC(idxNoNoise & isfinite(Ts.stepDeltaCohPC)));

  nNoiseCohLevels(iSession) = numel(noiseCoh);
  nNoNoiseCohLevels(iSession) = numel(noNoiseCoh);

  if ~isempty(noiseCoh)
    minNoiseCohPC(iSession) = min(noiseCoh);
    maxNoiseCohPC(iSession) = max(noiseCoh);
  end
  if ~isempty(noNoiseCoh)
    minNoNoiseCohPC(iSession) = min(noNoiseCoh);
    maxNoNoiseCohPC(iSession) = max(noNoiseCoh);

    % Descriptive diagnostic only: do the observed no-noise proportions
    % include values on both sides of 75%?
    [G, ~] = findgroups(Ts.stepDeltaCohPC(idxNoNoise));
    pByCoh = splitapply(@mean, Ts.correct(idxNoNoise), G);
    brackets75NoNoise(iSession) = ...
      any(pByCoh <= 0.75) && any(pByCoh >= 0.75);
  end
end

sessionSummary = table( ...
  sessionIndex, sessionID, animal, filePath, ...
  nINCTrials, nNoiseINCTrials, nNoNoiseINCTrials, ...
  nNoiseCohLevels, nNoNoiseCohLevels, ...
  minNoiseCohPC, maxNoiseCohPC, ...
  minNoNoiseCohPC, maxNoNoiseCohPC, ...
  pCorrectINC, pCorrectNoiseINC, pCorrectNoNoiseINC, ...
  brackets75NoNoise);
end

%% ------------------------------------------------------------------------
function printTrialInventory(sessionSummary, conditionSummary)

fprintf('\nIDR side-gain trial inventory: valid INC trials\n');
fprintf('Selected daily sessions: %d\n', height(sessionSummary));
fprintf('Total INC trials: %d\n', sum(sessionSummary.nINCTrials));
fprintf('  Noise trials: %d\n', sum(sessionSummary.nNoiseINCTrials));
fprintf('  No-noise trials: %d\n', sum(sessionSummary.nNoNoiseINCTrials));

fprintf(['Sessions with no-noise trials: %d/%d\n'], ...
  sum(sessionSummary.nNoNoiseINCTrials > 0), height(sessionSummary));
fprintf(['Sessions whose observed no-noise points bracket 75%%: %d/%d\n\n'], ...
  sum(sessionSummary.brackets75NoNoise), height(sessionSummary));

for iSession = 1:height(sessionSummary)
  S = sessionSummary(iSession, :);
  fprintf('%s  %s\n', S.animal, S.sessionID);
  fprintf('  INC total: %d; noise: %d; no-noise: %d\n', ...
    S.nINCTrials, S.nNoiseINCTrials, S.nNoNoiseINCTrials);

  C = conditionSummary(conditionSummary.sessionIndex == S.sessionIndex, :);

  Cnoise = C(C.hasStepNoise, :);
  if isempty(Cnoise)
    fprintf('  Noise:    none\n');
  else
    fprintf('  Noise:\n');
    for iRow = 1:height(Cnoise)
      fprintf('    %7.2f%% coh: %4d trials, %4d correct, %5.1f%%\n', ...
        Cnoise.stepDeltaCohPC(iRow), Cnoise.nTrials(iRow), ...
        Cnoise.nCorrect(iRow), 100 * Cnoise.pCorrect(iRow));
    end
  end

  CnoNoise = C(~C.hasStepNoise, :);
  if isempty(CnoNoise)
    fprintf('  No noise: none\n');
  else
    fprintf('  No noise:\n');
    for iRow = 1:height(CnoNoise)
      fprintf('    %7.2f%% coh: %4d trials, %4d correct, %5.1f%%\n', ...
        CnoNoise.stepDeltaCohPC(iRow), CnoNoise.nTrials(iRow), ...
        CnoNoise.nCorrect(iRow), 100 * CnoNoise.pCorrect(iRow));
    end
  end

  fprintf('  Observed no-noise points bracket 75%%: %s\n\n', ...
    string(S.brackets75NoNoise));
end
end

%% ------------------------------------------------------------------------
function F = fitSessionPsychometrics(T, sessionSummary, beta, lapse, target)
% Fit one threshold per daily session using valid INC no-noise trials only.

n = height(sessionSummary);
sessionIndex = sessionSummary.sessionIndex;
sessionID = sessionSummary.sessionID;
animal = sessionSummary.animal;
nTrials = zeros(n,1);
nCohLevels = zeros(n,1);
threshold = nan(n,1);
CI95Low = nan(n,1);
CI95High = nan(n,1);
negLogLikelihood = nan(n,1);
exitflag = nan(n,1);
thresholdInsideSampledRange = false(n,1);
profileLowerBoundHit = false(n,1);
profileUpperBoundHit = false(n,1);
maxAbsBinResidual = nan(n,1);

for k = 1:n
  idx = T.sessionIndex == sessionIndex(k) & ~T.hasStepNoise;
  coh = double(T.stepDeltaCohPC(idx));
  correct = double(T.correct(idx));
  use = isfinite(coh) & isfinite(correct) & coh >= 0;
  coh = coh(use);
  correct = correct(use);
  nTrials(k) = numel(correct);
  nCohLevels(k) = numel(unique(coh));
  if nTrials(k) == 0 || nCohLevels(k) < 2
    continue
  end

  fit = fitOneFixedShapePsychometric(coh, correct, beta, lapse, target);
  threshold(k) = fit.threshold;
  CI95Low(k) = fit.CI95(1);
  CI95High(k) = fit.CI95(2);
  negLogLikelihood(k) = fit.negLogLikelihood;
  exitflag(k) = fit.exitflag;
  thresholdInsideSampledRange(k) = threshold(k) >= min(coh) && threshold(k) <= max(coh);
  profileLowerBoundHit(k) = fit.profileLowerBoundHit;
  profileUpperBoundHit(k) = fit.profileUpperBoundHit;

  [G, level] = findgroups(coh);
  pObs = splitapply(@mean, correct, G);
  pFit = fixedShapeWeibullP(level, threshold(k), beta, lapse, target);
  maxAbsBinResidual(k) = max(abs(pObs-pFit));
end

CI95Width = CI95High-CI95Low;
relativeCI95Width = CI95Width./threshold;
F = table(sessionIndex, sessionID, animal, nTrials, nCohLevels, threshold, ...
  CI95Low, CI95High, CI95Width, relativeCI95Width, negLogLikelihood, exitflag, ...
  thresholdInsideSampledRange, profileLowerBoundHit, profileUpperBoundHit, ...
  maxAbsBinResidual);
end

%% ------------------------------------------------------------------------
function fit = fitOneFixedShapePsychometric(coh, correct, beta, lapse, target)

positive = coh(coh>0);
if isempty(positive)
  error('analyzeIDRSideGains:NoPositiveCoherence', ...
    'Psychometric fit requires positive coherence.');
end
lb = max(1e-3,min(positive)/20);
ub = max(5*max(coh),100*lb);
objective = @(logThreshold) psychometricNLL(exp(logThreshold),coh,correct,beta,lapse,target);
opts = optimset('Display','off','TolX',1e-9);
[logHat,nll,exitflag] = fminbnd(objective,log(lb),log(ub),opts);
hat = exp(logHat);

delta95 = 1.920729410347062;
logGrid = linspace(log(lb),log(ub),1600)';
grid = exp(logGrid);
profile = arrayfun(objective,logGrid);
inside = profile <= nll+delta95;

if any(inside)
  i1 = find(inside,1,'first');
  i2 = find(inside,1,'last');
  lo = profileCrossing(grid,profile,i1,-1,nll+delta95);
  hi = profileCrossing(grid,profile,i2,+1,nll+delta95);
else
  i1 = 1; i2 = numel(grid); lo = nan; hi = nan;
end

fit = struct('threshold',hat,'CI95',[lo hi],'negLogLikelihood',nll, ...
  'exitflag',exitflag,'thresholdBounds',[lb ub], ...
  'profileLowerBoundHit',i1==1,'profileUpperBoundHit',i2==numel(grid));
end

%% ------------------------------------------------------------------------
function xCross = profileCrossing(x,y,iInside,direction,targetY)
iOutside = iInside+direction;
if iOutside<1 || iOutside>numel(x)
  xCross=x(iInside); return
end
x1=x(iInside); x2=x(iOutside); y1=y(iInside); y2=y(iOutside);
if ~isfinite(y1) || ~isfinite(y2) || y1==y2
  xCross=x1; return
end
f=(targetY-y1)/(y2-y1);
f=min(max(f,0),1);
xCross=x1+f*(x2-x1);
end

%% ------------------------------------------------------------------------
function nll = psychometricNLL(threshold,coh,correct,beta,lapse,target)
p=fixedShapeWeibullP(coh,threshold,beta,lapse,target);
p=min(max(p,eps),1-eps);
nll=-sum(correct.*log(p)+(1-correct).*log(1-p));
end

%% ------------------------------------------------------------------------
function p = fixedShapeWeibullP(coh,threshold,beta,lapse,target)
% Weibull rising from 0.5 to 1-lapse, parameterized by target threshold.
if target<=0.5 || target>=1-lapse
  error('analyzeIDRSideGains:BadTargetPerformance', ...
    'Target performance must lie between 0.5 and 1-lapse.');
end
ratio=(1-lapse-target)/(0.5-lapse);
alpha=threshold/((-log(ratio))^(1/beta));
coh=max(double(coh),0);
p=1-lapse-(0.5-lapse).*exp(-(coh./alpha).^beta);
end

%% ------------------------------------------------------------------------
function printPsychometricSummary(F)
valid=isfinite(F.threshold);
fprintf('\nFixed-shape session psychometric fits: no-noise INC trials\n');
fprintf('  Finite thresholds: %d/%d\n',sum(valid),height(F));
fprintf('  Threshold inside sampled range: %d/%d\n', ...
  sum(F.thresholdInsideSampledRange & valid),sum(valid));
if any(valid)
  fprintf(['  Threshold: mean %.2f, SD %.2f, median %.2f, range %.2f-%.2f\n'], ...
    mean(F.threshold(valid)),std(F.threshold(valid)),median(F.threshold(valid)), ...
    min(F.threshold(valid)),max(F.threshold(valid)));
  fprintf('  Median 95%% profile interval width: %.2f coherence\n', ...
    median(F.CI95Width(valid),'omitnan'));
  fprintf('  Median relative interval width: %.2f\n', ...
    median(F.relativeCI95Width(valid),'omitnan'));
  fprintf('  Profile interval touched a search bound: %d/%d\n', ...
    sum((F.profileLowerBoundHit|F.profileUpperBoundHit)&valid),sum(valid));
end
end

%% ------------------------------------------------------------------------
function plotSessionPsychometrics(C,F,pdfPath,beta,lapse,target)

if isfile(pdfPath), delete(pdfPath); end
nRows=3; nCols=4; nPerPage=nRows*nCols;
nPages=ceil(height(F)/nPerPage);

for page=1:nPages
  fig=figure('Color','w','Visible','off','Units','inches','Position',[1 1 11 8.5]);
  tl=tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
  first=(page-1)*nPerPage+1;
  last=min(page*nPerPage,height(F));
  title(tl,sprintf('IDR no-noise INC psychometrics: sessions %d-%d of %d', ...
    first,last,height(F)),'FontWeight','bold');

  for k=first:last
    ax=nexttile(tl);
    plotOneSessionPsychometric(ax,C,F(k,:),beta,lapse,target);
  end
  if page==1
    exportgraphics(fig,pdfPath,'ContentType','vector');
  else
    exportgraphics(fig,pdfPath,'ContentType','vector','Append',true);
  end
  close(fig);
end
end

%% ------------------------------------------------------------------------
function plotOneSessionPsychometric(ax,C,F,beta,lapse,target)

hold(ax,'on');
D=C(C.sessionIndex==F.sessionIndex & ~C.hasStepNoise,:);
if isempty(D)
  axis(ax,'off'); title(ax,F.sessionID,'Interpreter','none'); return
end

se=sqrt(D.pCorrect.*(1-D.pCorrect)./D.nTrials);
errorbar(ax,D.stepDeltaCohPC,D.pCorrect,se,'ko','MarkerFaceColor','k', ...
  'MarkerSize',4,'LineStyle','none','CapSize',0);

if isfinite(F.threshold)
  values=[D.stepDeltaCohPC;F.CI95High;F.threshold];
  xMax=max(values(isfinite(values)))*1.08;
  xGrid=linspace(0,xMax,300);
  plot(ax,xGrid,fixedShapeWeibullP(xGrid,F.threshold,beta,lapse,target), ...
    '-','LineWidth',1.2);
  xline(ax,F.threshold,':');
  yline(ax,target,':');
end

for j=1:height(D)
  text(ax,D.stepDeltaCohPC(j),D.pCorrect(j),sprintf(' %d',D.nTrials(j)), ...
    'FontSize',6,'VerticalAlignment','bottom');
end
title(ax,sprintf('%s | %s',F.animal,erase(F.sessionID,'.mat')), ...
  'Interpreter','none','FontSize',6.5);
xlabel(ax,'Step coherence (%)'); ylabel(ax,'P(correct)');
ylim(ax,[0.45 1.02]); grid(ax,'on'); box(ax,'off');
txt=sprintf(['n=%d, levels=%d\nc_{75}=%.2f [%.2f, %.2f]\n' ...
  'inside range=%s\nmax residual=%.3f'],F.nTrials,F.nCohLevels, ...
  F.threshold,F.CI95Low,F.CI95High,string(F.thresholdInsideSampledRange), ...
  F.maxAbsBinResidual);
text(ax,0.02,0.98,txt,'Units','normalized','VerticalAlignment','top','FontSize',6);
end

%% ------------------------------------------------------------------------
function [alignedPsychometric, T] = fitAlignedPsychometrics( ...
  T, sessionPsychFits, targetPerformance, lapseBounds, sensitivityMaxRelativeCIWidth)
% Align each session by its fitted c75 threshold, then fit a common Weibull
% to all valid INC no-noise trials. A secondary fit excludes sessions whose
% threshold profile interval is extremely broad; the all-session fit remains
% primary.

thresholdBySession = nan(height(sessionPsychFits), 1);
thresholdBySession(sessionPsychFits.sessionIndex) = sessionPsychFits.threshold;
T.alignedCoh = T.stepDeltaCohPC ./ thresholdBySession(T.sessionIndex);

idxBase = ~T.hasStepNoise & isfinite(T.alignedCoh) & isfinite(T.correct);
fitAll = fitAlignedWeibull(T.alignedCoh(idxBase), T.correct(idxBase), ...
  targetPerformance, lapseBounds);
fitAll.label = 'All sessions';
fitAll.nSessions = numel(unique(T.sessionIndex(idxBase)));
fitAll.binned = binPsychometric(T.alignedCoh(idxBase), T.correct(idxBase), 12);

useSession = isfinite(sessionPsychFits.threshold) & ...
  isfinite(sessionPsychFits.relativeCI95Width) & ...
  sessionPsychFits.relativeCI95Width <= sensitivityMaxRelativeCIWidth;
useSessionIndex = sessionPsychFits.sessionIndex(useSession);
idxSensitivity = idxBase & ismember(T.sessionIndex, useSessionIndex);

fitSensitivity = fitAlignedWeibull( ...
  T.alignedCoh(idxSensitivity), T.correct(idxSensitivity), ...
  targetPerformance, lapseBounds);
fitSensitivity.label = sprintf('Relative CI width <= %.2f', ...
  sensitivityMaxRelativeCIWidth);
fitSensitivity.nSessions = numel(unique(T.sessionIndex(idxSensitivity)));
fitSensitivity.binned = binPsychometric( ...
  T.alignedCoh(idxSensitivity), T.correct(idxSensitivity), 12);

alignedPsychometric = struct();
alignedPsychometric.targetPerformance = targetPerformance;
alignedPsychometric.lapseBounds = lapseBounds(:)';
alignedPsychometric.sensitivityMaxRelativeCIWidth = sensitivityMaxRelativeCIWidth;
alignedPsychometric.primary = fitAll;
alignedPsychometric.sensitivity = fitSensitivity;
alignedPsychometric.nExcludedSensitivitySessions = sum(~useSession);
alignedPsychometric.excludedSensitivitySessionIDs = ...
  sessionPsychFits.sessionID(~useSession);
end

%% ------------------------------------------------------------------------
function fit = fitAlignedWeibull(coh, correct, targetPerformance, lapseBounds)

coh = double(coh(:));
correct = double(correct(:));
valid = isfinite(coh) & coh >= 0 & isfinite(correct);
coh = coh(valid);
correct = correct(valid);

if numel(correct) < 10
  error('analyzeIDRSideGains:TooFewAlignedTrials', ...
    'Too few trials for aligned Weibull fit.');
end

% Parameters are log(alpha), log(beta), and lapse. The initial alpha is
% chosen so that c75=1 at beta=2 and lapse=0.01.
beta0 = 2;
lapse0 = min(max(0.01, lapseBounds(1) + 1e-6), lapseBounds(2) - 1e-6);
alpha0 = alphaForThreshold(1, targetPerformance, beta0, lapse0);
x0 = [log(alpha0); log(beta0); lapse0];

lb = [log(0.05); log(0.25); lapseBounds(1)];
ub = [log(5); log(10); lapseBounds(2)];
objective = @(x) alignedWeibullNLL(x, coh, correct);

opts = optimoptions('fmincon', 'Display', 'off', ...
  'Algorithm', 'interior-point', 'OptimalityTolerance', 1e-9, ...
  'StepTolerance', 1e-10, 'MaxFunctionEvaluations', 5000, ...
  'MaxIterations', 2000);
[xHat, nll, exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, [], opts);

alpha = exp(xHat(1));
betaWeibull = exp(xHat(2));
lapse = xHat(3);
threshold = thresholdForAlpha(alpha, targetPerformance, betaWeibull, lapse);

% Observed-information covariance in fitted coordinates:
%   x = [log(alpha), log(beta), lapse].
hessian = finiteDifferenceHessian(objective, xHat);
covarianceFitCoordinates = nan(3);
covariance = nan(3);
correlation = nan(3);
SE = nan(3, 1);
CI95 = nan(3, 2);

if all(isfinite(hessian), 'all')
  hessian = (hessian + hessian') / 2;
  if all(eig(hessian) > 0)
    covarianceFitCoordinates = inv(hessian);

    % Transform covariance to [threshold, betaWeibull, lapse].
    transform = @(x) [ ...
      thresholdForAlpha(exp(x(1)), targetPerformance, exp(x(2)), x(3)); ...
      exp(x(2)); ...
      x(3)];
    J = finiteDifferenceJacobian(transform, xHat);
    covariance = J * covarianceFitCoordinates * J';
    covariance = (covariance + covariance') / 2;

    variance = diag(covariance);
    validVariance = isfinite(variance) & variance > 0;
    SE(validVariance) = sqrt(variance(validVariance));
    estimates = [threshold; betaWeibull; lapse];
    CI95(validVariance, :) = estimates(validVariance) + ...
      [-1 1] .* (1.96 * SE(validVariance));

    d = sqrt(diag(covariance));
    correlation = covariance ./ (d * d');
  end
end

fit = struct();
fit.alpha = alpha;
fit.betaWeibull = betaWeibull;
fit.lapse = lapse;
fit.threshold = threshold;
fit.thresholdPerformance = targetPerformance;
fit.negLogLikelihood = nll;
fit.exitflag = exitflag;
fit.nTrials = numel(correct);
fit.meanCorrect = mean(correct);
fit.fitCoordinateNames = ["logAlpha", "logBetaWeibull", "lapse"];
fit.parameterNames = ["threshold", "betaWeibull", "lapse"];
fit.hessian = hessian;
fit.covarianceFitCoordinates = covarianceFitCoordinates;
fit.covariance = covariance;
fit.correlation = correlation;
fit.SE = SE;
fit.CI95 = CI95;
end

%% ------------------------------------------------------------------------
function nll = alignedWeibullNLL(x, coh, correct)
alpha = exp(x(1));
betaWeibull = exp(x(2));
lapse = x(3);
p = weibullPAlpha(coh, alpha, betaWeibull, lapse);
p = min(max(p, eps), 1 - eps);
nll = -sum(correct .* log(p) + (1 - correct) .* log(1 - p));
end

%% ------------------------------------------------------------------------
function p = weibullPAlpha(coh, alpha, betaWeibull, lapse)
coh = max(double(coh), 0);
p = 1 - lapse - (0.5 - lapse) .* ...
  exp(-power(coh ./ alpha, betaWeibull));
end

%% ------------------------------------------------------------------------
function alpha = alphaForThreshold(threshold, targetPerformance, betaWeibull, lapse)
ratio = (1 - lapse - targetPerformance) / (0.5 - lapse);
alpha = threshold ./ ((-log(ratio)) .^ (1 ./ betaWeibull));
end

%% ------------------------------------------------------------------------
function threshold = thresholdForAlpha(alpha, targetPerformance, betaWeibull, lapse)
ratio = (1 - lapse - targetPerformance) / (0.5 - lapse);
threshold = alpha .* ((-log(ratio)) .^ (1 ./ betaWeibull));
end

%% ------------------------------------------------------------------------
function B = binPsychometric(coh, correct, nBins)

edges = linspace(min(coh), max(coh), nBins + 1);
binCenter = 0.5 * (edges(1:end-1) + edges(2:end));
nTrials = zeros(nBins, 1);
nCorrect = zeros(nBins, 1);
pCorrect = nan(nBins, 1);

for iBin = 1:nBins
  idx = coh >= edges(iBin) & coh < edges(iBin + 1);
  if iBin == nBins
    idx = coh >= edges(iBin) & coh <= edges(iBin + 1);
  end
  nTrials(iBin) = sum(idx);
  nCorrect(iBin) = sum(correct(idx));
  if nTrials(iBin) > 0
    pCorrect(iBin) = nCorrect(iBin) / nTrials(iBin);
  end
end

B = table(binCenter(:), nTrials, nCorrect, pCorrect, ...
  'VariableNames', {'alignedCoh', 'nTrials', 'nCorrect', 'pCorrect'});
B = B(B.nTrials > 0, :);
end

%% ------------------------------------------------------------------------
function H = finiteDifferenceHessian(objective, x)

x = x(:);
n = numel(x);
H = nan(n);
h = 3e-4 .* max(1, abs(x));
f0 = objective(x);

for i = 1:n
  ei = zeros(n, 1);
  ei(i) = h(i);
  fp = objective(x + ei);
  fm = objective(x - ei);
  H(i, i) = (fp - 2 * f0 + fm) / h(i)^2;

  for j = i+1:n
    ej = zeros(n, 1);
    ej(j) = h(j);
    fpp = objective(x + ei + ej);
    fpm = objective(x + ei - ej);
    fmp = objective(x - ei + ej);
    fmm = objective(x - ei - ej);
    H(i, j) = (fpp - fpm - fmp + fmm) / (4 * h(i) * h(j));
    H(j, i) = H(i, j);
  end
end
end

%% ------------------------------------------------------------------------
function J = finiteDifferenceJacobian(fun, x)

x = x(:);
y0 = fun(x);
nOut = numel(y0);
nIn = numel(x);
J = nan(nOut, nIn);
h = 1e-5 .* max(1, abs(x));

for i = 1:nIn
  ei = zeros(nIn, 1);
  ei(i) = h(i);
  J(:, i) = (fun(x + ei) - fun(x - ei)) / (2 * h(i));
end
end

%% ------------------------------------------------------------------------
function printPsychometricPassChange(previousFit, currentFit, previousPass, currentPass)

fprintf(['  Change pass %d -> %d: threshold %+.5f, beta %+.5f, ' ...
  'lapse %+.6f, NLL %+.3f\n'], ...
  previousPass, currentPass, ...
  currentFit.threshold - previousFit.threshold, ...
  currentFit.betaWeibull - previousFit.betaWeibull, ...
  currentFit.lapse - previousFit.lapse, ...
  currentFit.negLogLikelihood - previousFit.negLogLikelihood);
end

%% ------------------------------------------------------------------------
function H = makePsychometricRefinementHistory( ...
  provisionalBeta, provisionalLapse, pass1Aligned, pass2Aligned, pass3Aligned)

pass = [];
sessionFitBeta = [];
sessionFitLapse = [];
alignedThreshold = [];
alignedBeta = [];
alignedLapse = [];
thresholdSE = [];
betaSE = [];
lapseSE = [];
negLogLikelihood = [];
nSessions = [];
nTrials = [];

appendPass(1, provisionalBeta, provisionalLapse, pass1Aligned.primary);
if ~isempty(pass2Aligned)
  appendPass(2, pass1Aligned.primary.betaWeibull, ...
    pass1Aligned.primary.lapse, pass2Aligned.primary);
end
if ~isempty(pass3Aligned)
  appendPass(3, pass2Aligned.primary.betaWeibull, ...
    pass2Aligned.primary.lapse, pass3Aligned.primary);
end

deltaThreshold = [nan; diff(alignedThreshold)];
deltaBeta = [nan; diff(alignedBeta)];
deltaLapse = [nan; diff(alignedLapse)];
deltaNLL = [nan; diff(negLogLikelihood)];

H = table(pass, sessionFitBeta, sessionFitLapse, ...
  alignedThreshold, alignedBeta, alignedLapse, ...
  thresholdSE, betaSE, lapseSE, negLogLikelihood, ...
  deltaThreshold, deltaBeta, deltaLapse, deltaNLL, ...
  nSessions, nTrials);

  function appendPass(passNumber, fixedBeta, fixedLapse, F)
    pass(end+1,1) = passNumber; %#ok<AGROW>
    sessionFitBeta(end+1,1) = fixedBeta; %#ok<AGROW>
    sessionFitLapse(end+1,1) = fixedLapse; %#ok<AGROW>
    alignedThreshold(end+1,1) = F.threshold; %#ok<AGROW>
    alignedBeta(end+1,1) = F.betaWeibull; %#ok<AGROW>
    alignedLapse(end+1,1) = F.lapse; %#ok<AGROW>
    if isfield(F, 'SE') && numel(F.SE) == 3
      thresholdSE(end+1,1) = F.SE(1); %#ok<AGROW>
      betaSE(end+1,1) = F.SE(2); %#ok<AGROW>
      lapseSE(end+1,1) = F.SE(3); %#ok<AGROW>
    else
      thresholdSE(end+1,1) = nan; %#ok<AGROW>
      betaSE(end+1,1) = nan; %#ok<AGROW>
      lapseSE(end+1,1) = nan; %#ok<AGROW>
    end
    negLogLikelihood(end+1,1) = F.negLogLikelihood; %#ok<AGROW>
    nSessions(end+1,1) = F.nSessions; %#ok<AGROW>
    nTrials(end+1,1) = F.nTrials; %#ok<AGROW>
  end
end

%% ------------------------------------------------------------------------
function printAlignedPsychometricSummary(P)

fprintf('\nAligned pooled no-noise INC psychometric fits\n');
printOneAlignedFit(P.primary);
fprintf('  Sensitivity exclusion: %d sessions\n', ...
  P.nExcludedSensitivitySessions);
printOneAlignedFit(P.sensitivity);
end

%% ------------------------------------------------------------------------
function printOneAlignedFit(F)
fprintf(['  %s: n=%d trials, %d sessions, threshold %.3f, ' ...
  'beta %.3f, lapse %.4f, NLL %.3f, exitflag %d\n'], ...
  F.label, F.nTrials, F.nSessions, F.threshold, ...
  F.betaWeibull, F.lapse, F.negLogLikelihood, F.exitflag);

if isfield(F, 'CI95') && all(isfinite(F.CI95), 'all')
  fprintf(['    95%% Wald intervals: threshold [%.3f, %.3f], ' ...
    'beta [%.3f, %.3f], lapse [%.4f, %.4f]\n'], ...
    F.CI95(1,1), F.CI95(1,2), ...
    F.CI95(2,1), F.CI95(2,2), ...
    F.CI95(3,1), F.CI95(3,2));
end
end

%% ------------------------------------------------------------------------
function plotAlignedPsychometric(P, pdfPath)

fig = figure('Color', 'w', 'Visible', 'off', ...
  'Units', 'inches', 'Position', [1 1 10 4.5]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'IDR aligned no-noise INC psychometrics', 'FontWeight', 'bold');

plotOneAlignedPsychometric(nexttile(tl), P.primary);
plotOneAlignedPsychometric(nexttile(tl), P.sensitivity);

exportgraphics(fig, pdfPath, 'ContentType', 'vector');
close(fig);
end

%% ------------------------------------------------------------------------
function plotOneAlignedPsychometric(ax, F)
hold(ax, 'on');
B = F.binned;
plot(ax, B.alignedCoh, B.pCorrect, 'ko', ...
  'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineStyle', 'none');
for iBin = 1:height(B)
  text(ax, B.alignedCoh(iBin), B.pCorrect(iBin), ...
    sprintf(' %d', B.nTrials(iBin)), 'FontSize', 7, ...
    'VerticalAlignment', 'bottom');
end

xMax = max([max(B.alignedCoh) F.threshold]) * 1.08;
xGrid = linspace(0, xMax, 400);
pGrid = weibullPAlpha(xGrid, F.alpha, F.betaWeibull, F.lapse);
plot(ax, xGrid, pGrid, '-', 'LineWidth', 1.4);
xline(ax, 1, '--', 'Alignment target', 'LabelVerticalAlignment', 'bottom');
xline(ax, F.threshold, ':', sprintf('fit %.3f', F.threshold));
yline(ax, F.thresholdPerformance, ':');

xlabel(ax, 'Coherence / session c_{75}');
ylabel(ax, 'P(correct)');
title(ax, F.label, 'Interpreter', 'none');
ylim(ax, [0.48 1.01]);
xlim(ax, [0 xMax]);
grid(ax, 'on');
box(ax, 'off');

if isfield(F, 'SE') && all(isfinite(F.SE))
  txt = sprintf(['%d trials, %d sessions\n' ...
    'threshold=%.3f \\pm %.3f\n' ...
    '\\beta=%.3f \\pm %.3f\n' ...
    '\\lambda=%.4f \\pm %.4f'], ...
    F.nTrials, F.nSessions, F.threshold, F.SE(1), ...
    F.betaWeibull, F.SE(2), F.lapse, F.SE(3));
else
  txt = sprintf(['%d trials, %d sessions\n' ...
    'threshold=%.3f\n\\beta=%.3f, \\lambda=%.4f'], ...
    F.nTrials, F.nSessions, F.threshold, F.betaWeibull, F.lapse);
end
text(ax, 0.98, 0.02, txt, 'Units', 'normalized', ...
  'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
  'FontSize', 8);
end
