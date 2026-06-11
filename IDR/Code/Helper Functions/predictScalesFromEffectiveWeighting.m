function pred = predictScalesFromEffectiveWeighting(offsetsDeg, mtModel, modelName, params, fixedBaseline)
% predictScalesFromEffectiveWeighting
%
% Predict normalized psychophysical scales from an effective MT-to-choice
% weighting function and a fixed MT forward model.
%
% Model:
%   R(delta)      = sum_phi wEff(phi) * rMT(phi; delta)
%   S_pred(delta) = R(delta) / R(0)
%
% Inputs:
%   offsetsDeg    : vector of probe offsets to predict
%   mtModel       : struct defining the MT forward model
%   modelName     : effective-weighting model name
%   params        : parameters for modelName
%   fixedBaseline : optional fixed DOG baseline
%
% Output:
%   pred          : predicted normalized scale at each offset
%
% Notes:
%   - The only fitted scientific object is wEff(phi).
%   - The MT model enters only as forward-model machinery used to map
%     wEff(phi) onto predicted scale.

if nargin < 5
    fixedBaseline = [];
end

offsetsDeg = offsetsDeg(:)';
phiDeg     = mtModel.phiDeg(:)';

wEff = makeEffectiveWeighting(phiDeg, modelName, params, fixedBaseline);

R0 = sum(wEff .* mtModel.shiftedResponseFun(phiDeg, 0));
if ~isfinite(R0) || abs(R0) < 1e-9
    pred = nan(size(offsetsDeg));
    return;
end

pred = nan(size(offsetsDeg));
for i = 1:numel(offsetsDeg)
  delta = offsetsDeg(i);
  R = sum(wEff .* mtModel.shiftedResponseFun(phiDeg, delta));
  pred(i) = R ./ R0;
end

end