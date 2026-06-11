function [kernels, kVars, kStats] = computeSessionKernelsFromNoiseMatrices( ...
  prefNoiseByPatch, probeNoiseByPatch, ...
  trialOutcomesAll, changeSidesAll, changeIndicesAll, ...
  prefCohNoisePC, probeCohNoisePC)
% computeSessionKernelsFromNoiseMatrices
%
% Recompute session kernels from saved noise matrices, mirroring makeKernels.
%
% Inputs:
%   prefNoiseByPatch   2 x T x nTrials
%   probeNoiseByPatch  2 x T x nTrials
%   trialOutcomesAll   1 x nTrials or nTrials x 1
%   changeSidesAll     1 x nTrials or nTrials x 1
%   changeIndicesAll   1 x nTrials or nTrials x 1   (1=DEC, 2=INC)
%   prefCohNoisePC     scalar
%   probeCohNoisePC    scalar
%
% Outputs:
%   kernels            5 x 2 x 2 x T
%   kVars              5 x 2 x 2
%   kStats             5 x 2 x 2 cell array

trialOutcomesAll = trialOutcomesAll(:)';
changeSidesAll   = changeSidesAll(:)';
changeIndicesAll = changeIndicesAll(:)';

m = size(prefNoiseByPatch, 2);

kernels = nan(5, 2, 2, m);
kVars   = nan(5, 2, 2);
kStats  = cell(5, 2, 2);

for sideType = 1:5
  if sideType == 1
    ampScale = sqrt(2);
  else
    ampScale = 1;
  end

  for s = 1:2   % 1=DEC, 2=INC
    useTrials = (changeIndicesAll == s);

    prefNoiseStep  = prefNoiseByPatch(:,:,useTrials);
    probeNoiseStep = probeNoiseByPatch(:,:,useTrials);
    trialOutcomesStep = trialOutcomesAll(useTrials);
    changeSidesStep   = changeSidesAll(useTrials);

    [prefNoise, probeNoise] = selectSideTypeMatrices( ...
      prefNoiseStep, probeNoiseStep, changeSidesStep, sideType);

    [kernels(sideType, s, 1, :), kVars(sideType, s, 1), kStats{sideType, s, 1}] = ...
      meanPsychKernel(prefNoise, trialOutcomesStep, ...
      changeSidesStep, ampScale * prefCohNoisePC);

    [kernels(sideType, s, 2, :), kVars(sideType, s, 2), kStats{sideType, s, 2}] = ...
      meanPsychKernel(probeNoise, trialOutcomesStep, ...
      changeSidesStep, ampScale * probeCohNoisePC);
  end
end

end