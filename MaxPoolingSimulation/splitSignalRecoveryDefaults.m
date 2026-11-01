function cfg = splitSignalRecoveryDefaults()
% SPLITSIGNALRECOVERYDEFAULTS  Defaults for fair split-signal model recovery.

base=splitSignalDefaults();
cfg=base;
cfg.SeparationDeg=90;
cfg.AllocationFractions=[0.50 0.625 0.75];
cfg.SignalLevelMultipliers=[0.80 1.00 1.20];
cfg.StepNoiseHalfAmplitude=7.071;
cfg.PrestepPedestalFraction=1;
cfg.AdaptationFraction=4/60;
cfg.BalancePrestep=false;

cfg.PoolingModels={'max','vector','energy'};
cfg.TrueInternalNoiseSD=8;
cfg.NCalibrationTrials=6000;
cfg.NStimulusTrials=9000;
cfg.NLinkCalibrationTrials=12000;
cfg.NRecoveryReplicates=100;
cfg.NGainReplicates=50;
cfg.TrialsPerSession=2000;
cfg.MaxRecoverySessions=10;
cfg.TargetRecoveryProbability=0.90;
cfg.ThresholdTolerance=0.15;
cfg.ThresholdMaxIterations=18;
cfg.ThresholdInitialBounds=[4 80];
cfg.MaximumTotalStepCoherence=85;
cfg.Seed=842791;
cfg.Verbose=true;
cfg.TrueLapse=0;
cfg.SelectionCriterion='bic';

cfg.CreateFigures=true;
cfg.Visible='on';
cfg.SaveOutputs=true;
cfg.OutputDir=fullfile(pwd,'SplitSignalModelRecovery_Output');
end
