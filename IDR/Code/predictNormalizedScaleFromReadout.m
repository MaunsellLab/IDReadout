function predScale = predictNormalizedScaleFromReadout(params, offsetsDeg, mtModel, varargin)
% predictNormalizedScaleFromReadout  Predict probe/preferred behavioral scale.
%
% S(delta) = <a, Delta m_eff(delta)> / <a, Delta m_eff(0)>

p = inputParser;
addParameter(p, 'TemplateMode', 'signed', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
templateMode = lower(char(string(p.Results.TemplateMode)));

[aPhi, ~] = evaluateReadoutDOG(mtModel.phiDeg, params);
refTemplate = mtReadoutTemplate(0, mtModel, 'TemplateMode', templateMode);
refVal = sum(aPhi .* refTemplate);

offsetsDeg = offsetsDeg(:)';
predScale = nan(size(offsetsDeg));
for i = 1:numel(offsetsDeg)
    probeTemplate = mtReadoutTemplate(offsetsDeg(i), mtModel, ...
        'TemplateMode', templateMode);
    probeVal = sum(aPhi .* probeTemplate);
    if isfinite(refVal) && abs(refVal) > 0
        predScale(i) = probeVal / refVal;
    end
end
end
