function cfg = temporalMaxDefaults()
% TEMPORALMAXDEFAULTS  Defaults for max-before-temporal-difference simulation.

paired=pairedProbeDefaults();
paired.OffsetsDeg=[10 25 45 60 75 90 120 135 170 180];
paired.InternalCorrelationLengthDeg=30;
paired.NTrials=12000;
paired.NCalibrationTrials=8000;
paired.Visible='off';
paired.CreateFigures=false;
paired.SaveOutputs=false;

cfg=struct();
cfg.PairedConfig=paired;
cfg.BaseDriftCoherences=[7.071 40];
cfg.ReferenceBaseDriftCoherence=7.071;
cfg.PreStepDurationSec=0.750;
cfg.StepDurationSec=0.250;
cfg.InternalStepNoiseSD=1;
cfg.CalibrateInternalNoise=true;
cfg.WinnerDirectionBinWidthDeg=5;
cfg.WinnerDirectionPlotOffsetDeg=120;
cfg.CreateFigures=true;
cfg.Visible='on';
cfg.SaveOutputs=true;
cfg.OutputDir=fullfile(pwd,'PairedProbeTemporalMax_Output');
end
