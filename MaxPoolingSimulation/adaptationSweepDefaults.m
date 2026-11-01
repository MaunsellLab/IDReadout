function cfg = adaptationSweepDefaults()
% ADAPTATIONSWEEPDEFAULTS  Defaults for prestep firing-dependent adaptation.

paired=pairedProbeDefaults();
paired.OffsetsDeg=[10 20 45 120 170];
paired.FixedStepCoherence=15.3;
paired.NoiseCalibrationOffsetDeg=120;
paired.NTrials=12000;
paired.NCalibrationTrials=6000;
paired.CalibrationTolerance=3e-3;
paired.CalibrationMaxIterations=28;
paired.Visible='off';
paired.CreateFigures=false;
paired.SaveOutputs=false;

cfg=struct();
cfg.PairedConfig=paired;
cfg.BaseDriftCoherences=[7.071 40];
cfg.AftereffectAt60Values=0:8;
cfg.AdaptationReferenceCoherence=60;
cfg.AdaptationWidthDeg=paired.TuningSigmaDeg;
cfg.ReferenceBaseDriftCoherence=7.071;
cfg.ReferenceAftereffectAt60=0;
cfg.EmpiricalIDSOffsetsDeg=[10 20 45 120 170];
cfg.EmpiricalIDSThresholds=[13.0 11.9 16.2 15.3 11.2];
cfg.EmpiricalIDSSessionCounts=[NaN NaN NaN 17 3];
cfg.EmpiricalIDRThresholdRange=[21 22];
cfg.ExpectedBaseThresholdIncrease=5;
cfg.WinnerDirectionBinWidthDeg=5;
cfg.WinnerDirectionPlotOffsetDeg=120;
cfg.WinnerDirectionPlotBase=40;
cfg.WinnerDirectionPlotAftereffects=[0 4 7];
cfg.CreateFigures=true;
cfg.Visible='on';
cfg.SaveOutputs=true;
cfg.OutputDir=fullfile(pwd,'PairedProbeAdaptationSweep_Output');
end
