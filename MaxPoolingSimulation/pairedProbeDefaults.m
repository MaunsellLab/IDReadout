function cfg = pairedProbeDefaults()
% PAIREDPROBEDEFAULTS  Defaults for runPairedProbeOffsetSweep.

cfg=struct();
cfg.OffsetsDeg=[0 5 10 15 20 25 30 45 60 75 90 105 120 135 150 160 170 175 179 180];
cfg.BankSpacingDeg=1;
cfg.TuningSigmaDeg=37.5;
cfg.CenterTuning=true;

cfg.ProbePedestal=7.071;
cfg.DriftPedestal=0;
cfg.BaseDriftCoherence=0;
cfg.DriftNoiseSD=10/sqrt(5);
cfg.ProbeNoiseSD=7.071/sqrt(5);

cfg.FixedStepCoherence=20;
cfg.CalibrateInternalNoise=true;
cfg.NoiseCalibrationOffsetDeg=90;
cfg.NoiseCalibrationYoked=false;
cfg.InternalNoiseSD=1;
cfg.InternalCorrelationLengthDeg=30;
cfg.EffectiveNoiseRank=NaN; % Computed by the primary runner after parsing.
cfg.RemoveInternalCircularMean=true;
cfg.NoiseGenerationBatchSize=2000;

cfg.TargetPerformance=0.75;
cfg.NTrials=20000;
cfg.NCalibrationTrials=12000;
cfg.CalibrationTolerance=2e-3;
cfg.CalibrationMaxIterations=35;
cfg.CalibrationInitialUpper=1;

cfg.CounterfactualPerturbation=0.25;
cfg.WinnerAssociationRadiusDeg=10;
cfg.Seed=91327;
cfg.Visible='on';
cfg.CreateFigures=true;
cfg.SaveOutputs=true;
cfg.OutputDir=fullfile(pwd,'PairedProbeOffsetSweep_Output');
end
