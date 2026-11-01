function value = randomDirectionMean(sigmaDeg)
%RANDOMDIRECTIONMEAN Mean circular-Gaussian response to random directions.
% The random-dot pool is uniform over the 360-degree direction circle.
sigmaDeg=double(sigmaDeg);
assert(isscalar(sigmaDeg)&&isfinite(sigmaDeg)&&sigmaDeg>0, ...
  'sigmaDeg must be a positive finite scalar.');
value=(sigmaDeg*sqrt(2*pi)/360)* ...
  erf(180/(sqrt(2)*sigmaDeg));
end
