%EXAMPLESYNTHETICFIT End-to-end executable example.
cfg=maxpoolfit.defaultConfig();
cfg.nNoiseDraws=32;
cfg.finalNoiseDraws=64;
cfg.batchSize=64;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',12);
truth=struct('adaptation',cfg.adaptation,'noiseSDPC',[13.5 8 10]);
cache=maxpoolfit.simulateChoices(cache,cfg.adaptation,truth.noiseSDPC,cfg);
fit=maxpoolfit.fitModel(cache,cfg);
prediction=maxpoolfit.predict(cache,fit,cfg); %#ok<NASGU>
disp(fit)
