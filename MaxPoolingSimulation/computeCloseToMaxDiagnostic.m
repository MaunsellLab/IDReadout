function diagnostic = computeCloseToMaxDiagnostic(sim, cfg)
% COMPUTECLOSETOMAXDIAGNOSTIC  Paired mechanism perturbation versus max gap.

R = sim.responseChanged;
[maxValue, winner] = max(R, [], 2);
RwithoutWinner = R;
linearWinner = sub2ind(size(R), (1:size(R,1))', winner);
RwithoutWinner(linearWinner) = -Inf;
secondValue = max(RwithoutWinner, [], 2);

n = size(R,2);
gap = bsxfun(@minus, R, maxValue);
for m = 1:n
    isWinner = winner == m;
    gap(isWinner,m) = R(isWinner,m) - secondValue(isWinner);
end

delta = cfg.DiagnosticPerturbation;
influence = zeros(size(R));
for m = 1:n
    influence(:,m) = (max(maxValue, R(:,m)+delta)-maxValue)/delta;
end

gapScale = max([cfg.InternalNoiseSD, cfg.ExternalNoiseSD, eps]);
edges = cfg.DiagnosticGapEdgesSD * gapScale;
centers = (edges(1:end-1)+edges(2:end))/2;
meanInfluence = nan(size(centers));
counts = zeros(size(centers));
flatGap = gap(:);
flatInfluence = influence(:);
for b = 1:numel(centers)
    take = flatGap >= edges(b) & flatGap < edges(b+1);
    counts(b) = sum(take);
    if counts(b) > 0
        meanInfluence(b) = mean(flatInfluence(take));
    end
end

diagnostic = struct();
diagnostic.n = sim.n;
diagnostic.stepCoherence = sim.stepCoherence;
diagnostic.gapBinCenters = centers;
diagnostic.meanCounterfactualInfluence = meanInfluence;
diagnostic.binCounts = counts;
diagnostic.perturbation = delta;
diagnostic.gapDefinition = 'R_m - max_{k~=m} R_k';
end

