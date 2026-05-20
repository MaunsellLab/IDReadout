function [kernels, kVars, kStats, hitStats, compStats, trialOutcomesByCond, changeSidesByCond, changeIndicesByCond] = ...
  computeSessionKernels(sessionData, trialIdx)
% computeSessionKernels
% Recomputes the full nSideTypes x 2 x 2 kernel set from saved patchwise noise matrices.
% If trialIdx is provided, trials are resampled/reordered accordingly.
% This function is used both for ordinary per-session analysis and for
% within-session bootstrap resampling.
%
% INPUT
%   sessionData.prefNoiseByPatch : 2 x m x nTrials
%   sessionData.probeNoiseByPatch: 2 x m x nTrials
%   sessionData.trialOutcomesAll : 1 x nTrials
%   sessionData.changeSidesAll   : 1 x nTrials
%   sessionData.changeIndicesAll : 1 x nTrials   (1=decrement, 2=increment)
%   sessionData.sessionProbeHeader : metadata struct for this derived probe session
%   sessionData.lr               : optional session L/R mapping struct
%   sessionData.sideTypeNames    : optional cell array of side-type names
%
%   trialIdx                     : optional bootstrap resample indices
%                                  if omitted, use all trials
%
% OUTPUT
%   kernels(sideType, stepType, streamType, frame)
%       sideType   : 1=diff, 2=change, 3=noChange, 4=L, 5=R, 6=RF, 7=Opp
%       stepType   : 1=decrement, 2=increment
%       streamType : 1=pref, 2=probe

if nargin < 2 || isempty(trialIdx)
  trialIdx = 1:size(sessionData.prefNoiseByPatch, 3);
end
assert(isfield(sessionData, 'sessionProbeHeader') && ~isempty(sessionData.sessionProbeHeader), ...
  'computeSessionKernels:MissingSessionProbeHeader', 'sessionData.sessionProbeHeader is required.');

sessionProbeHeader = sessionData.sessionProbeHeader;
prefNoiseAll     = sessionData.prefNoiseByPatch(:, :, trialIdx);
probeNoiseAll    = sessionData.probeNoiseByPatch(:, :, trialIdx);
trialOutcomesAll = sessionData.trialOutcomesAll(trialIdx);
changeSidesAll   = sessionData.changeSidesAll(trialIdx);
changeIndicesAll = sessionData.changeIndicesAll(trialIdx);

% Pull amplitudes from the authoritative per-probe session metadata.
% probeCohNoisePC is the single-stream probe noise amplitude stored in the
% experimental file. probeNoiseByPatch has already been converted to the
% effective yoked-probe perturbation, so kernel estimation must use the
% corresponding combined probe amplitude.
prefCohNoisePC  = sessionProbeHeader.prefCohNoisePC.data(1);
probeCohNoisePC = sessionProbeHeader.probeCohNoisePC.data(1);

if ~isfinite(prefCohNoisePC) || ~isfinite(probeCohNoisePC) || probeCohNoisePC <= 0
  error('computeSessionKernels:BadNoiseAmplitude', ...
    'Invalid pref/probe coherence noise amplitudes: pref=%g, probe=%g.', ...
    prefCohNoisePC, probeCohNoisePC);
end

nYokedProbeStreams = probeStreamCountFromSessionProbeHeader(sessionProbeHeader);
combinedProbeCohNoisePC = nYokedProbeStreams * probeCohNoisePC;

% side types
[~, sideTypeNames] = sideTypeIndex();
nSideTypes = numel(sideTypeNames);

% sessionwise L/R mapping
if isfield(sessionData, 'lr') && ~isempty(sessionData.lr)
  lr = sessionData.lr;
else
  % backward compatibility for older saved sessionData structs
  if isfield(sessionData, 'trials') && ~isempty(sessionData.trials)
    lr = sessionLRMap(sessionData.trials);
  else
    error('computeSessionKernels:MissingLRMap', ...
      ['sessionData.lr is missing, and no trials field is available to ' ...
       'reconstruct the L/R mapping.']);
  end
end

m = size(prefNoiseAll, 2);
prefNoiseAll  = forcePatchNoise3D(prefNoiseAll, m);
probeNoiseAll = forcePatchNoise3D(probeNoiseAll, m);

kernels = nan(nSideTypes, 2, 2, m);
kVars   = nan(nSideTypes, 2, 2);
kStats  = repmat(struct('nRFCorrect', 0, 'nRFWrong', 0, 'nCorrect', 0, 'nWrong', 0, ...
  'sumCorrect', [], 'sumWrong', [], 'sigma2', nan, 'emptyOutcomeClass', nan), nSideTypes, 2, 2);

trialOutcomesByCond = cell(nSideTypes, 2);
changeSidesByCond   = cell(nSideTypes, 2);
changeIndicesByCond = cell(nSideTypes, 2);

for sideType = 1:nSideTypes
  if sideType == 1
    ampScale = sqrt(2);   % diff kernels
  else
    ampScale = 1;
  end

  for s = 1:2   % 1=DEC, 2=INC
    useTrials = (changeIndicesAll == s);
    prefNoiseStep  = prefNoiseAll(:, :, useTrials);
    probeNoiseStep = probeNoiseAll(:, :, useTrials);
    prefNoiseStep  = forcePatchNoise3D(prefNoiseStep, m);
    probeNoiseStep = forcePatchNoise3D(probeNoiseStep, m);
    trialOutcomesByCond{sideType, s} = trialOutcomesAll(useTrials);
    changeSidesByCond{sideType, s}   = changeSidesAll(useTrials);
    changeIndicesByCond{sideType, s} = changeIndicesAll(useTrials);

    [prefMat, probeMat] = selectSideTypeMatrices( ...
      prefNoiseStep, probeNoiseStep, changeSidesByCond{sideType, s}, sideType, lr);

    [kernels(sideType, s, 1, :), kVars(sideType, s, 1), kStats(sideType, s, 1)] = meanPsychKernel( ...
      prefMat, trialOutcomesByCond{sideType, s}, changeSidesByCond{sideType, s}, ampScale * prefCohNoisePC);

    [kernels(sideType, s, 2, :), kVars(sideType, s, 2), kStats(sideType, s, 2)] = meanPsychKernel( ...
      probeMat, trialOutcomesByCond{sideType, s}, changeSidesByCond{sideType, s}, ampScale * combinedProbeCohNoisePC);
  end
end

% The number of trials and hits overall are the same for all sideTypes.
% We also save counts for RF and for Left.
hitStats = struct;
for s = 1:2
  hitStats.nTrials(s) = numel(trialOutcomesByCond{1, s});
  hitStats.nHits(s)   = sum(trialOutcomesByCond{1, s} == 0);

  % RF counts
  idxRF = (changeSidesByCond{1, s} == 0);
  hitStats.nRFTrials(s) = sum(idxRF);
  hitStats.nRFHits(s)   = sum(trialOutcomesByCond{1, s} == 0 & idxRF);

  % Left counts
  if lr.leftIsRF
    idxLeft = (changeSidesByCond{1, s} == 0);   % RF is Left
  else
    idxLeft = (changeSidesByCond{1, s} == 1);   % Opp is Left
  end
  hitStats.nLeftTrials(s) = sum(idxLeft);
  hitStats.nLeftHits(s)   = sum(trialOutcomesByCond{1, s} == 0 & idxLeft);
end
msPerVFrame = 1000.0 / sessionProbeHeader.frameRateHz.data(1);

% ---- Comparison statistics ----
%
% Raw statistics use the kernels exactly as measured, with pref and probe
% kernels expressed using their actual session noise amplitudes.
%
% Normalized statistics rescale the probe kernels to the pref-noise
% amplitude convention. This is the appropriate comparison for probe/pref
% ratios and scales when pref and probe noise amplitudes differ.
%
% If the measured probe kernel amplitude scales with noise variance, then
% converting a probe kernel measured at probeCohNoisePC to the equivalent
% prefCohNoisePC convention requires:
%
%   K_probe_norm = K_probe_raw * (prefCohNoisePC / probeCohNoisePC)^2
%
% The returned raw kernels are not modified.

compStats = struct;

% Raw compStats
[compStats.rawIntegrals, compStats.rawR, compStats.rawRVar] = ...
  kernelIntegral(kernels, kVars, msPerVFrame);
[compStats.rawScale, compStats.rawScaleSEM, compStats.rawFitR2, compStats.rawSSE] = ...
  kernelScaleFit(kernels, msPerVFrame);

% Normalized compStats
probeNormFactor = (prefCohNoisePC / combinedProbeCohNoisePC)^2;
kernelsNorm = kernels;
kVarsNorm   = kVars;

kernelsNorm(:, :, 2, :) = kernelsNorm(:, :, 2, :) * probeNormFactor;
kVarsNorm(:, :, 2)      = kVarsNorm(:, :, 2) * probeNormFactor^2;

[compStats.normIntegrals, compStats.normR, compStats.normRVar] = kernelIntegral(kernelsNorm, kVarsNorm, msPerVFrame);
[compStats.normScale, compStats.normScaleSEM, compStats.normFitR2, compStats.normSSE] = ...
  kernelScaleFit(kernelsNorm, msPerVFrame);

compStats.normInfo = struct();
compStats.normInfo.prefCohNoisePC = prefCohNoisePC;
compStats.normInfo.probeCohNoisePC = probeCohNoisePC;
compStats.normInfo.probeNormFactor = probeNormFactor;
compStats.normInfo.nYokedProbeStreams = nYokedProbeStreams;
compStats.normInfo.combinedProbeCohNoisePC = combinedProbeCohNoisePC;
compStats.normInfo.method = ...
  ['probeNoiseByPatch stores effective yoked probe noise; probe kernels are ' ...
   'normalized by (prefCohNoisePC/(nYokedProbeStreams*probeCohNoisePC))^2 ' ...
   'before computing normalized ratios/scales'];
% Legacy aliases: preserve old downstream behavior for now.
compStats.kIntegrals = compStats.rawIntegrals;
compStats.R          = compStats.rawR;
compStats.RVar       = compStats.rawRVar;
compStats.scale      = compStats.rawScale;
compStats.scaleSEM   = compStats.rawScaleSEM;
compStats.fitR2      = compStats.rawFitR2;
compStats.sse        = compStats.rawSSE;

end

%%
function n = probeStreamCountFromSessionProbeHeader(sessionProbeHeader)
% Number of yoked probe streams represented by the probe-noise variable.
% Current convention:
%   0 < probeDirDeg < 180 : paired yoked streams at +/- probeDirDeg
%   probeDirDeg == 180   : legacy single opposite-direction stream

assert(isfield(sessionProbeHeader, 'probeDirDeg'), ...
  'computeSessionKernels:MissingProbeDir', ...
  'Cannot determine probe stream count because sessionProbeHeader.probeDirDeg is missing.');

probeDirDeg = abs(double(sessionProbeHeader.probeDirDeg));

if probeDirDeg > 0 && probeDirDeg < 180
  n = 2;
elseif abs(probeDirDeg - 180) < 1e-9
  n = 1;
else
  error('computeSessionKernels:UnsupportedProbeDir', ...
    'Unsupported probeDirDeg for probe normalization: %g.', probeDirDeg);
end
end

function X = forcePatchNoise3D(X, nFrames)
% forcePatchNoise3D  Preserve 2 x m x nTrials shape when nTrials == 1.
%
% MATLAB drops trailing singleton dimensions in many indexing operations.
% Patchwise noise matrices should always be treated as:
%
%   2 x nFrames x nTrials

if isempty(X)
  X = nan(2, nFrames, 0);
  return;
end

if ndims(X) == 2
  assert(size(X,1) == 2 && size(X,2) == nFrames, ...
    'computeSessionKernels:BadPatchNoiseShape', ...
    'Expected patch noise matrix to be 2 x %d or 2 x %d x nTrials.', ...
    nFrames, nFrames);

  X = reshape(X, 2, nFrames, 1);
end
end