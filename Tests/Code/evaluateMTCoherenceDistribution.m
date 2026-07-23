function populationResponse = evaluateMTCoherenceDistribution( ...
  componentDirectionsDeg,componentCoherence,mtModel)
% evaluateMTCoherenceDistribution  MT response to coherence distributions.
%
% componentCoherence is observations-by-components. A vector is treated as
% one observation. Output is observations-by-MT-neurons.

directions = double(componentDirectionsDeg(:));
coherence = double(componentCoherence);
if isvector(coherence)
  coherence = coherence(:)';
end
if size(coherence,2)~=numel(directions)
  error('evaluateMTCoherenceDistribution:SizeMismatch', ...
    'Coherence columns must match the number of component directions.');
end
if any(~isfinite(coherence),'all')
  error('evaluateMTCoherenceDistribution:NonfiniteCoherence', ...
    'Coherence values must be finite.');
end
templates = mtPopulationTemplate(directions,mtModel);
populationResponse = coherence*templates;
end
