function noise = generateCircularGPNoise(nTrials,nMechanisms,cfg)
% GENERATECIRCULARGPNOISE  FFT-generated smooth zero-mean circular noise field.

if cfg.InternalNoiseSD==0
    noise=zeros(nTrials,nMechanisms,'single');
    return;
end
if cfg.InternalCorrelationLengthDeg<=0
    noise=single(cfg.InternalNoiseSD*randn(nTrials,nMechanisms));
    if cfg.RemoveInternalCircularMean
        noise=noise-mean(noise,2);
        noise=noise/sqrt(1-1/nMechanisms);
    end
    return;
end

spacing=360/nMechanisms;
d=(0:nMechanisms-1)*spacing;
d=min(d,360-d);
kappa=(180/(pi*cfg.InternalCorrelationLengthDeg))^2;
firstRow=exp(kappa*(cosd(d)-1));
lambda=max(real(fft(firstRow)),0);
sqrtLambda=single(sqrt(lambda));

noise=zeros(nTrials,nMechanisms,'single');
batch=cfg.NoiseGenerationBatchSize;
for first=1:batch:nTrials
    last=min(first+batch-1,nTrials);
    z=randn(last-first+1,nMechanisms,'single');
    field=real(ifft(fft(z,[],2).*sqrtLambda,[],2));
    if cfg.RemoveInternalCircularMean
        field=field-mean(field,2);
    end
    noise(first:last,:)=field;
end

if cfg.RemoveInternalCircularMean
    varianceAfterRemoval=max(1-lambda(1)/nMechanisms,eps);
else
    varianceAfterRemoval=1;
end
noise=noise*single(cfg.InternalNoiseSD/sqrt(varianceAfterRemoval));
end

