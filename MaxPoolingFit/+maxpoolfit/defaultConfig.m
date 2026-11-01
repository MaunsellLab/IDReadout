function cfg = defaultConfig()
%DEFAULTCONFIG Frozen first-fit configuration for the max-pooling model.

cfg.modelVersion = "DenseHardMax-v5-shared-prestep-adaptation-random-dot-redistribution";
cfg.cacheSchemaVersion = "MaxPoolingTrialCache-v2";
cfg.decisionRule = "absolute-step-max";
cfg.coherenceResponseMode = "random-dot-redistribution";

cfg.sigmaMTDeg = 37.5;
cfg.sigmaAdaptDeg = 37.5;
cfg.noiseCorrLengthDeg = 30;
cfg.independentNoiseFraction = 0;
cfg.noiseVarianceExponent = 0;
cfg.directionSpacingDeg = 2;
cfg.directionGridDeg = 0:cfg.directionSpacingDeg:(360-cfg.directionSpacingDeg);

cfg.adaptation = 0.56;
cfg.adaptationSelection = "shared 256-draw profile-likelihood optimum";
cfg.adaptationExponent = 0;
cfg.adaptationGainStrength = 0;
cfg.adaptationReferencePC = 40;
cfg.noiseBoundsPC = [0.25 40];
cfg.nNoiseDraws = 64;
cfg.finalNoiseDraws = 256;
cfg.noiseSeed = 1729;
cfg.probabilityFloor = 1e-8;
cfg.batchSize = 128;

cfg.useAllCrossPairs = true;
cfg.useAntitheticDraws = true;
cfg.primarySelector = @(T) T.valid & T.hasStepNoise & T.isIncrement;
cfg.verbose = true;
end
