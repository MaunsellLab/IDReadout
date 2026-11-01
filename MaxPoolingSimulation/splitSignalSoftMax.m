function value = splitSignalSoftMax(response,temperature)
% SPLITSIGNALSOFTMAX  Stable hard-max or log-sum-exp landscape summary.

if temperature<=0
    value=max(response,[],2);
    return;
end
scaled=response/single(temperature);
peak=max(scaled,[],2);
value=single(temperature)*(peak+log(sum(exp(bsxfun(@minus,scaled,peak)),2)));
end
