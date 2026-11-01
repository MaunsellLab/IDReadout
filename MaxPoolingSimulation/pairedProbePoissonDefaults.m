function cfg = pairedProbePoissonDefaults()
% PAIREDPROBEPOISSONDEFAULTS  Defaults for the pooled-Poisson variance stage.

cfg=struct();
paired=pairedProbeDefaults();
paired.OffsetsDeg=[10 25 45 60 75 90 120 135 170 180];
paired.InternalCorrelationLengthDeg=30;
paired.NTrials=12000;
paired.NCalibrationTrials=8000;
paired.Visible='off';
paired.CreateFigures=false;
paired.SaveOutputs=false;
cfg.PairedConfig=paired;

cfg.BaseDriftCoherences=[7.071 40];
cfg.ReferenceBaseDriftCoherence=7.071;
cfg.RandomMotionRateSpPerSec=40;
cfg.AlignedRateSlopeSpPerSecPerPct=0.5;
cfg.PreStepDurationSec=0.750;
cfg.StepDurationSec=0.250;
cfg.MinimumRateSpPerSec=0.1;
cfg.FanoFactor=1;

cfg.PoolCalibrationInitialUpper=1;
cfg.PoolCalibrationMaximum=1e6;
cfg.IncludeHomoscedasticReference=true;
cfg.Visible='on';
cfg.CreateFigures=true;
cfg.SaveOutputs=true;
cfg.OutputDir=fullfile(pwd,'PairedProbePoissonSweep_Output');
end
