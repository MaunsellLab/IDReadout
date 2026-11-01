function G = coherenceTuning(mechanismDeg,streamDeg,sigmaDeg,cfg)
%COHERENCETUNING Directional response produced by redistributing random dots.
% In the corrected fixed-density RDK model, coherent dots replace uniformly
% directed random dots. Their net response is therefore G(theta-d)-mean(G).
% Configurations predating this contract retain the historical additive mode.
mode="additive";
if nargin>=4 && isfield(cfg,'coherenceResponseMode')
  mode=string(cfg.coherenceResponseMode);
end
assert(isscalar(mode),'coherenceResponseMode must be scalar.');
G=maxpoolfit.circularTuning(mechanismDeg,streamDeg,sigmaDeg);
switch mode
  case "additive"
    % Historical model: coherent streams were added without removing dots
    % from the isotropic random pool.
  case "random-dot-redistribution"
    G=G-maxpoolfit.randomDirectionMean(sigmaDeg);
  otherwise
    error('maxpoolfit:InvalidCoherenceResponseMode', ...
      'Unknown coherenceResponseMode "%s".',mode);
end
end
