function effectiveRank = circularNoiseEffectiveRank(nMechanisms,cfg)
% CIRCULARNOISEEFFECTIVERANK  Participation ratio of the implemented field.
%
% The calculation includes the removal and rescaling of the circular mean in
% generateCircularGPNoise.  It therefore describes the covariance actually
% supplied to the max-pooling mechanisms, not merely the nominal von-Mises
% kernel specified by InternalCorrelationLengthDeg.

if cfg.InternalCorrelationLengthDeg<=0
    if cfg.RemoveInternalCircularMean
        effectiveRank=nMechanisms-1;
    else
        effectiveRank=nMechanisms;
    end
    return;
end

spacing=360/nMechanisms;
d=(0:nMechanisms-1)*spacing;
kappa=(180/(pi*cfg.InternalCorrelationLengthDeg))^2;
rho=exp(kappa*(cosd(d)-1));

if cfg.RemoveInternalCircularMean
    commonFraction=mean(rho);
    rho=(rho-commonFraction)/(1-commonFraction);
end

eigenvalues=max(real(fft(rho)),0);
effectiveRank=sum(eigenvalues)^2/sum(eigenvalues.^2);
end
