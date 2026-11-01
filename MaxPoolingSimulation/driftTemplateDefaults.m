function cfg = driftTemplateDefaults()
% DRIFTTEMPLATEDEFAULTS  Defaults for perceptual drift-template cancellation.

paired=pairedProbeDefaults();
paired.OffsetsDeg=[120 170 180];
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
cfg.CancellationFractions=[0 0.5 0.75 0.825 0.875 0.925 0.975 1];
cfg.ReferenceBaseDriftCoherence=7.071;
cfg.ReferenceCancellationFraction=1;
cfg.WinnerDirectionBinWidthDeg=5;
cfg.WinnerDirectionPlotOffsetDeg=120;
cfg.WinnerDirectionPlotBase=40;
cfg.WinnerDirectionPlotFractions=[0 0.825 1];
cfg.CreateFigures=true;
cfg.Visible='on';
cfg.SaveOutputs=true;
cfg.OutputDir=fullfile(pwd,'PairedProbeDriftTemplate_Output');
end
