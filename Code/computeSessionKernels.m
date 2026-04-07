function [kernels, kVars, kStats, hitStats, compStats, trialOutcomesByCond, changeSidesByCond, changeIndicesByCond] = ...
  computeSessionKernels(sessionData, trialIdx)
% sessionKernelFromSaved
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
%   sessionData.header           : header struct
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

prefNoiseAll     = sessionData.prefNoiseByPatch(:, :, trialIdx);
probeNoiseAll    = sessionData.probeNoiseByPatch(:, :, trialIdx);
trialOutcomesAll = sessionData.trialOutcomesAll(trialIdx);
changeSidesAll   = sessionData.changeSidesAll(trialIdx);
changeIndicesAll = sessionData.changeIndicesAll(trialIdx);

% pull amplitudes from header
prefCohNoisePC  = sessionData.header.blockStatus.data.prefCohNoisePC;
probeCohNoisePC = sessionData.header.blockStatus.data.probeCohNoisePC;

% side types
if isfield(sessionData, 'sideTypeNames') && ~isempty(sessionData.sideTypeNames)
  sideTypeNames = sessionData.sideTypeNames;
else
  sideTypeNames = {'diff', 'change', 'noChange', 'L', 'R', 'RF', 'Opp'};
end
nSideTypes = numel(sideTypeNames);

% sessionwise L/R mapping
if isfield(sessionData, 'lr') && ~isempty(sessionData.lr)
  lr = sessionData.lr;
else
  % backward compatibility for older saved sessionData structs
  if isfield(sessionData, 'trials') && ~isempty(sessionData.trials)
    lr = sessionLRMap(sessionData.trials);
  else
    error('sessionKernelFromSaved:MissingLRMap', ...
      ['sessionData.lr is missing, and no trials field is available to ' ...
       'reconstruct the L/R mapping.']);
  end
end

m = size(prefNoiseAll, 2);

kernels = nan(nSideTypes, 2, 2, m);
kVars   = nan(nSideTypes, 2, 2);
kStats  = repmat(struct('nRFCorrect', 0, 'nRFWrong', 0, 'nCorrect', 0, 'nWrong', 0, ...
  'sumCorrect', [], 'sumWrong', [], 'sigma2', nan), nSideTypes, 2, 2);

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

    trialOutcomesByCond{sideType, s} = trialOutcomesAll(useTrials);
    changeSidesByCond{sideType, s}   = changeSidesAll(useTrials);
    changeIndicesByCond{sideType, s} = changeIndicesAll(useTrials);

    [prefMat, probeMat] = selectSideTypeMatrices( ...
      prefNoiseStep, probeNoiseStep, changeSidesByCond{sideType, s}, sideType, lr);

    [kernels(sideType, s, 1, :), kVars(sideType, s, 1), kStats(sideType, s, 1)] = meanPsychKernel( ...
      prefMat, trialOutcomesByCond{sideType, s}, changeSidesByCond{sideType, s}, ampScale * prefCohNoisePC);

    [kernels(sideType, s, 2, :), kVars(sideType, s, 2), kStats(sideType, s, 2)] = meanPsychKernel( ...
      probeMat, trialOutcomesByCond{sideType, s}, changeSidesByCond{sideType, s}, ampScale * probeCohNoisePC);
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

% hitStats = struct;
% for s = 1:2
%   idxStep = (changeIndicesAll == s);
%   hitStats.nTrials(s) = sum(idxStep);
%   hitStats.nHits(s)   = sum(trialOutcomesAll(idxStep) == 0);
% 
%   idxRF = idxStep & (changeSidesAll == 0);
%   hitStats.nRFTrials(s) = sum(idxRF);
%   hitStats.nRFHits(s)   = sum(trialOutcomesAll(idxRF) == 0);
% end

msPerVFrame = 1000.0 / sessionData.header.frameRateHz.data(1);

compStats = struct;
[compStats.kIntegrals, compStats.R, compStats.RVar] = kernelIntegral(kernels, kVars, msPerVFrame);
[compStats.scale, compStats.scaleSEM, compStats.fitR2, compStats.sse] = kernelScaleFit(kernels, msPerVFrame);

end