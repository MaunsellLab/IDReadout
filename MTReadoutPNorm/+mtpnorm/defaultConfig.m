function cfg = defaultConfig()
%DEFAULTCONFIG Initial full-response, pooled pre/post-difference model.
cfg.modelVersion = "MTReadoutPNorm-v2-local";
cfg.animal = "Meetz";
cfg.sigmaMTDeg = 37.5;
cfg.directionSpacingDeg = 1;
cfg.noiseCorrLengthDeg = 30;
cfg.spontaneousResponse = 0;
cfg.preWindowMS = [500 750];
cfg.stepWindowMS = [750 1000];
cfg.fitDraws = 128;
cfg.finalDraws = 512;
cfg.seed = 1729;
cfg.noiseBounds = [0 40];
cfg.noiseGrid = [0 .5 1 2 4 8 16 32 40];
cfg.probabilityFloor = 1e-8;
cfg.trialBatchSize = 128;
cfg.checkpointFile = "";
cfg.keepRawStepNoise = false;
% Units: a 100%-coherent preferred stimulus contributes 100 above the
% optional spontaneous offset. Noise SD is in these SAME response units,
% per unit per 250-ms estimate; it is not directly a coherence threshold.
end
