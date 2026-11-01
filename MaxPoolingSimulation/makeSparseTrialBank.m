function bank = makeSparseTrialBank(n, nTrials, cfg, seed)
% MAKESPARSETRIALBANK  Draw fixed random inputs for one n condition.
%
% Rows are trials. Columns are perfectly aligned stream/mechanism directions.
% Changed and no-change patches have independent stimulus and mechanism noise.

oldState = rng;
cleanup = onCleanup(@() rng(oldState));
rng(seed, 'twister');

directionsDeg = (0:n-1) * (360/n);
delta = directionsDeg(:) - directionsDeg(:)'; % stream j by mechanism m
[W, meanT] = centeredMTSensitivity(delta, cfg.TuningSigmaDeg, cfg.CenterTuning);

xChanged = cfg.ExternalNoiseSD * randn(nTrials, n);
xNoChange = cfg.ExternalNoiseSD * randn(nTrials, n);
epsilonChanged = cfg.InternalNoiseSD * randn(nTrials, n);
epsilonNoChange = cfg.InternalNoiseSD * randn(nTrials, n);

bank = struct();
bank.n = n;
bank.nTrials = nTrials;
bank.directionsDeg = directionsDeg;
bank.offsetDeg = min(directionsDeg, 360-directionsDeg);
bank.W = W;
bank.meanT = meanT;
bank.xChanged = xChanged;
bank.xNoChange = xNoChange;
bank.predictors = xChanged - xNoChange;
if isscalar(cfg.StreamPedestal)
    pedestal = repmat(cfg.StreamPedestal,1,n);
else
    pedestal = reshape(cfg.StreamPedestal,1,[]);
    if numel(pedestal) ~= n
        error('StreamPedestal must be scalar or have one value per stream.');
    end
end
pedestalResponse = pedestal*W;
bank.pedestal = pedestal;
bank.pedestalResponse = pedestalResponse;
bank.responseChangedBase = bsxfun(@plus,xChanged*W + epsilonChanged,pedestalResponse);
bank.responseNoChange = bsxfun(@plus,xNoChange*W + epsilonNoChange,pedestalResponse);
bank.seed = seed;
end
