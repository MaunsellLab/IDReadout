function pred = predictScalesFromReadout(offsetsDeg, mt, modelName, params)
% predictScalesFromReadout
%
% Predict normalized behavioral scales:
%   R(delta) = sum_phi w(phi) * m(phi; delta)
%   S_pred(delta) = R(delta) / R(0)
%
% Inputs:
%   offsetsDeg : vector of probe offsets to predict
%   mt         : struct from makeMTKernel
%   modelName  : readout model name
%   params     : parameters for modelName
%
% Output:
%   pred       : predicted scale at each offset

offsetsDeg = offsetsDeg(:)';
phiDeg     = mt.phiDeg(:)';

w = makeReadout(phiDeg, modelName, params);

R0 = sum(w .* mt.shiftedGFun(phiDeg, 0));

pred = nan(size(offsetsDeg));
for i = 1:numel(offsetsDeg)
    delta = offsetsDeg(i);
    R = sum(w .* mt.shiftedGFun(phiDeg, delta));
    pred(i) = R ./ R0;
end

end