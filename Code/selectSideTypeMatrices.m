function [prefNoise, probeNoise] = selectSideTypeMatrices( ...
    prefNoiseByPatch, probeNoiseByPatch, changeSides, sideType, lr)
% selectSideTypeMatrices
%
% Select trial x frame matrices for one sideType.
%
% Inputs:
%   prefNoiseByPatch   2 x m x nTrials   (patch 1 = RF, patch 2 = Opp)
%   probeNoiseByPatch  2 x m x nTrials
%   changeSides        1 x nTrials or nTrials x 1   (0 = RF change, 1 = Opp change)
%   sideType           1..7
%   lr                 struct from sessionLRMap()
%
% sideTypes:
%   1 diff
%   2 change
%   3 noChange
%   4 L
%   5 R
%   6 RF
%   7 Opp
%
% Outputs:
%   prefNoise          m x nTrials
%   probeNoise         m x nTrials

  SIDE_DIFF     = 1;
  SIDE_CHANGE   = 2;
  SIDE_NOCHANGE = 3;
  SIDE_LEFT     = 4;
  SIDE_RIGHT    = 5;
  SIDE_RF       = 6;
  SIDE_OPP      = 7;

  [nPatches, m, nTrials] = size(prefNoiseByPatch);
  if nPatches ~= 2
    error('selectSideTypeMatrices:BadInput', ...
      'Expected prefNoiseByPatch to have size 2 x m x nTrials.');
  end
  if ~isequal(size(probeNoiseByPatch), [2, m, nTrials])
    error('selectSideTypeMatrices:BadInput', ...
      'prefNoiseByPatch and probeNoiseByPatch must have the same size.');
  end

  changeSides = reshape(changeSides, 1, []);
  if numel(changeSides) ~= nTrials
    error('selectSideTypeMatrices:BadInput', ...
      'Number of changeSides entries must match nTrials.');
  end

  % Trialwise patch identities
  changePatch   = 1 + changeSides;   % RF->1, Opp->2
  noChangePatch = 3 - changePatch;   % RF<->Opp

  % Sessionwise L/R mapping onto RF/Opp
  if lr.leftIsRF
    leftPatch  = 1;   % RF
    rightPatch = 2;   % Opp
  else
    leftPatch  = 2;   % Opp
    rightPatch = 1;   % RF
  end

  prefNoise  = nan(m, nTrials);
  probeNoise = nan(m, nTrials);

  switch sideType
    case SIDE_DIFF
      for t = 1:nTrials
        prefNoise(:, t)  = squeeze(prefNoiseByPatch(changePatch(t), :, t) - ...
                                   prefNoiseByPatch(noChangePatch(t), :, t)).';
        probeNoise(:, t) = squeeze(probeNoiseByPatch(changePatch(t), :, t) - ...
                                   probeNoiseByPatch(noChangePatch(t), :, t)).';
      end

    case SIDE_CHANGE
      for t = 1:nTrials
        prefNoise(:, t)  = squeeze(prefNoiseByPatch(changePatch(t), :, t)).';
        probeNoise(:, t) = squeeze(probeNoiseByPatch(changePatch(t), :, t)).';
      end

    case SIDE_NOCHANGE
      for t = 1:nTrials
        prefNoise(:, t)  = squeeze(prefNoiseByPatch(noChangePatch(t), :, t)).';
        probeNoise(:, t) = squeeze(probeNoiseByPatch(noChangePatch(t), :, t)).';
      end

    case SIDE_LEFT
      for t = 1:nTrials
        prefNoise(:, t)  = squeeze(prefNoiseByPatch(leftPatch, :, t)).';
        probeNoise(:, t) = squeeze(probeNoiseByPatch(leftPatch, :, t)).';
      end

    case SIDE_RIGHT
      for t = 1:nTrials
        prefNoise(:, t)  = squeeze(prefNoiseByPatch(rightPatch, :, t)).';
        probeNoise(:, t) = squeeze(probeNoiseByPatch(rightPatch, :, t)).';
      end

    case SIDE_RF
      prefNoise  = squeeze(prefNoiseByPatch(1, :, :));
      probeNoise = squeeze(probeNoiseByPatch(1, :, :));

    case SIDE_OPP
      prefNoise  = squeeze(prefNoiseByPatch(2, :, :));
      probeNoise = squeeze(probeNoiseByPatch(2, :, :));

    otherwise
      error('selectSideTypeMatrices:BadSideType', ...
        'Unknown sideType value: %d', sideType);
  end

  % Ensure m x nTrials orientation if squeeze collapsed dimensions oddly
  if size(prefNoise, 1) ~= m && size(prefNoise, 2) == m
    prefNoise = prefNoise.';
  end
  if size(probeNoise, 1) ~= m && size(probeNoise, 2) == m
    probeNoise = probeNoise.';
  end
end