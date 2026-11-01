function cfg = splitSignalDefaults()
% SPLITSIGNALDEFAULTS  Defaults for the split-signal design simulation.

cfg=struct();
cfg.SeparationsDeg=[90 120];
cfg.AllocationFractions=[0.50 0.625 0.75];
cfg.TotalStepCoherences=[12 16 20 24 28 36 44 52 60];
cfg.StepNoiseHalfAmplitudes=[5 7.071];
cfg.PrestepPedestalFractions=[1 0.5 0];
cfg.IncludeBalancedPrestep=[false true];

cfg.BankSpacingDeg=1;
cfg.TuningSigmaDeg=37.5;
cfg.CenterTuning=true;
cfg.InternalNoiseSD=8;
cfg.InternalCorrelationLengthDeg=30;
cfg.RemoveInternalCircularMean=true;
cfg.NoiseGenerationBatchSize=1000;
cfg.NNoiseFrames=5;
cfg.AdaptationFractions=[0 4/60];
cfg.TemplateCancellationFraction=1;
cfg.SoftmaxTemperatures=[0 2 4 6 8 12 16];
cfg.PrimarySoftmaxTemperature=8;

cfg.TargetPerformance=0.75;
cfg.NThresholdTrials=12000;
cfg.NDesignTrials=8000;
cfg.NPowerReplicates=250;
cfg.TrialsPerSession=2000;
cfg.MaxSessions=20;
cfg.PowerCriterionDeltaLL=5;
cfg.TargetPower=0.90;
cfg.Seed=731944;

cfg.PrimaryStepNoiseHalfAmplitude=7.071;
cfg.PrimaryPrestepPedestalFraction=1;
cfg.PrimaryAdaptationFraction=4/60;
cfg.PrimaryBalancePrestep=false;
cfg.PrimaryPowerSeparationDeg=90;

cfg.CreateFigures=true;
cfg.Visible='on';
cfg.SaveOutputs=true;
cfg.OutputDir=fullfile(pwd,'SplitSignalDesign_Output');
end
