function results = prefKernelSessionVariance(baseFolder)
% prefKernelSessionVariance
%
% Measure probeSession-to-probeSession variability in the absolute
% preferred-direction kernel integral, pooling probeSessions across all
% probe offsets.
%
% The measurement is:
%   side type   = diff (change patch minus no-change patch)
%   step type   = increment
%   stream type = preferred direction
%
% OUTPUT
%   results : one row per probeSession, with variables:
%       fileName
%       probeDirDeg
%       nTrials
%       signedIntegralPC
%       integralPC
%       integralSEMPC
%
% The function prints the across-probeSession mean and SD of the absolute
% integral and produces:
%   1. Histogram of integral amplitudes
%   2. Integral amplitude versus number of trials
%
% Kernel integrals and SEMs are in percentage-points coherence.
%
% Example:
%   results = prefKernelSessionVariance();

% cleanupObj = initProjectPath(); %#ok<NASGU>
if nargin < 1 || isempty(baseFolder)
  baseFolder = domainFolder(mfilename('fullpath'));
end

% Kernel array conventions used by computeSessionKernels:
%   sideType:   2 = change
%   stepType:   1 = decrement, 2 = increment
%   streamType: 1 = preferred, 2 = probe
sideType = 2;
stepType  = 2;
prefStream   = 1;

% Locate every saved probeSession, regardless of offset.
searchPattern = fullfile(baseFolder, 'Data', 'Probe*', 'ProbeSessions', '*.mat');
files = dir(searchPattern);
if isempty(files)
  error('prefKernelSessionVariance:NoFiles', 'No probeSession files found using:\n%s', searchPattern);
end

nFiles = numel(files);
fileName          = strings(nFiles, 1);
probeDirDeg       = nan(nFiles, 1);
nTrials           = nan(nFiles, 1);
signedIntegralPC  = nan(nFiles, 1);
integralPC        = nan(nFiles, 1);
integralSEMPC     = nan(nFiles, 1);
keep              = false(nFiles, 1);

[preStepMS, intStartMS, intDurMS] = integralWindowMS();

fprintf('Examining %d probeSession files...\n', nFiles);

for f = 1:nFiles
  filePath = fullfile(files(f).folder, files(f).name);

  S = load(filePath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr', 'prefNoiseByPatch', ...
    'probeNoiseByPatch', 'trialOutcomesAll', 'changeSidesAll', 'chosenSidesAll', 'changeIndicesAll');
  requiredFields = {'sessionProbeHeader', 'prefNoiseByPatch', 'probeNoiseByPatch','trialOutcomesAll', ...
    'changeSidesAll', 'chosenSidesAll', 'changeIndicesAll'};
  missing = requiredFields(~isfield(S, requiredFields));
  if ~isempty(missing)
    warning('prefKernelSessionVariance:IncompleteFile', 'Skipping %s; missing: %s', ...
      files(f).name, strjoin(missing, ', '));
    continue
  end

  % Current files should contain sessionHeader. For compatibility, use
  % sessionProbeHeader for timing information if necessary.
  if isfield(S, 'sessionHeader') && ~isempty(S.sessionHeader)
    sessionHeader = S.sessionHeader;
  else
    sessionHeader = S.sessionProbeHeader;
  end

  sessionData = struct;
  sessionData.sessionHeader      = sessionHeader;
  sessionData.sessionProbeHeader = S.sessionProbeHeader;

  if isfield(S, 'sideTypeNames')
    sessionData.sideTypeNames = S.sideTypeNames;
  end

  if isfield(S, 'lr')
    sessionData.lr = S.lr;
  end

  sessionData.prefNoiseByPatch  = S.prefNoiseByPatch;
  sessionData.probeNoiseByPatch = S.probeNoiseByPatch;
  sessionData.trialOutcomesAll  = S.trialOutcomesAll;
  sessionData.changeSidesAll    = S.changeSidesAll;
  sessionData.changeIndicesAll  = S.changeIndicesAll;

  [kernels, kVars, ~, hitStats] = computeSessionKernels(sessionData);

  frameRateHz = headerScalar(sessionHeader, 'frameRateHz');
  msPerVFrame = 1000 / frameRateHz;

  % Match the indexing convention in kernelIntegral().
  firstIndex = round((preStepMS + intStartMS) / msPerVFrame);
  lastIndex  = round((preStepMS + intStartMS + intDurMS) / ...
                     msPerVFrame);
  intIndices = firstIndex:lastIndex;

  if firstIndex < 1 || lastIndex > size(kernels, 4)
    warning('prefKernelSessionVariance:BadIntegralWindow', ...
      'Skipping %s; integration window lies outside the kernel.', ...
      files(f).name);
    continue
  end

  nIntegralBins = numel(intIndices);

  kernelTrace = squeeze(kernels(sideType, stepType, prefStream, intIndices));
  thisSignedIntegral = mean(kernelTrace, 'omitnan');

  % kVars is the variance of each kernel bin. The variance of the mean
  % over N bins is approximated by kVars/N, as in kernelIntegral().
  thisKernelVar = kVars(sideType, stepType, prefStream);
  thisIntegralSEM = sqrt(thisKernelVar / nIntegralBins);
  if ~isfinite(thisSignedIntegral) || ~isfinite(thisIntegralSEM)
    warning('prefKernelSessionVariance:NonfiniteEstimate', ...
      'Skipping %s; integral or SEM is nonfinite.', files(f).name);
    continue
  end
  fileName(f)         = string(files(f).name);
  probeDirDeg(f)      = double(S.sessionProbeHeader.probeDirDeg);
  signedIntegralPC(f) = thisSignedIntegral;
  integralPC(f)       = abs(thisSignedIntegral);
  integralSEMPC(f)    = thisIntegralSEM;

  % Increment trials contributing to the increment kernel.
  if isfield(hitStats, 'nTrials') && ...
      numel(hitStats.nTrials) >= stepType
    nTrials(f) = hitStats.nTrials(stepType);
  else
    nTrials(f) = sum(S.changeIndicesAll == stepType);
  end

  keep(f) = true;
end

% Remove skipped or unusable probeSessions.
fileName         = fileName(keep);
probeDirDeg      = probeDirDeg(keep);
nTrials          = nTrials(keep);
signedIntegralPC = signedIntegralPC(keep);
integralPC       = integralPC(keep);
integralSEMPC    = integralSEMPC(keep);

if isempty(integralPC)
  error('prefKernelSessionVariance:NoValidSessions', ...
    'No valid probeSessions were found.');
end

results = table(fileName, probeDirDeg, nTrials, ...
                signedIntegralPC, integralPC, integralSEMPC);

% Sort primarily by parent/session filename and then by probe direction.
results = sortrows(results, {'fileName', 'probeDirDeg'});

nProbeSessions = height(results);
meanIntegral   = mean(results.integralPC);
sdIntegral     = std(results.integralPC);

fprintf('\nPreferred-direction increment kernel integrals\n');
fprintf('  ProbeSessions: %d\n', nProbeSessions);
fprintf('  Mean:          %.4f %% coherence\n', meanIntegral);
fprintf('  SD:            %.4f %% coherence\n', sdIntegral);
fprintf('  CV:            %.3f\n', sdIntegral / meanIntegral);
fprintf('  Median SEM:    %.4f %% coherence\n\n', ...
        median(results.integralSEMPC));

%% Histogram

figure('Name', 'Preferred-kernel session variability', ...
       'Color', 'w');

tiledlayout(1, 2, 'TileSpacing', 'compact', ...
                  'Padding', 'compact');

nexttile;

histogram(results.integralPC, 'BinMethod', 'fd');

xMax = max(results.integralPC + results.integralSEMPC);
if ~isfinite(xMax) || xMax <= 0
  xMax = 1;
else
  xMax = 1.05 * xMax;
end

xlim([0 xMax]);
xline(meanIntegral, '--', ...
  sprintf('Mean = %.3f', meanIntegral), ...
  'LabelVerticalAlignment', 'middle');

xlabel('Absolute preferred-kernel integral (% coherence)');
ylabel('Number of probeSessions');
title(sprintf('Mean %.3f, SD %.3f, n = %d', ...
      meanIntegral, sdIntegral, nProbeSessions));
box off;

%% Integral versus trial count

nexttile;

errorbar(results.nTrials, results.integralPC, ...
         results.integralSEMPC, ...
         'o', ...
         'LineStyle', 'none', ...
         'CapSize', 0);

yline(0, ':');
ylim([0 xMax]);

xlabel('Number of increment trials');
ylabel('Absolute preferred-kernel integral (% coherence)');
title('ProbeSession integral versus trial count');
box off;
[~, sideTypeNames] = sideTypeIndex();
sgtitle(sprintf('Preferred-direction increment %s kernel amplitude', sideTypeNames{sideType}));

end


function value = headerScalar(header, fieldName)
% Return a scalar from either a direct numeric header field or a legacy
% header field containing a .data value.

if ~isfield(header, fieldName)
  error('prefKernelSessionVariance:MissingHeaderField', ...
    'Header lacks field %s.', fieldName);
end

value = header.(fieldName);

if isstruct(value)
  if ~isfield(value, 'data')
    error('prefKernelSessionVariance:BadHeaderField', ...
      'Header field %s is a struct without a data field.', fieldName);
  end
  value = value.data;
end

value = double(value(1));

if ~isfinite(value)
  error('prefKernelSessionVariance:BadHeaderValue', ...
    'Header field %s is not finite.', fieldName);
end
end