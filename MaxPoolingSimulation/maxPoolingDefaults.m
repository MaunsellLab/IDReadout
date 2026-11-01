function cfg = maxPoolingDefaults()
% MAXPOOLINGDEFAULTS  Defaults for runMaxPoolingNSweep.
%
% Noise units are normalized: ExternalNoiseSD = 1 is one expected SD of a
% kernel-weighted scalar coherence-noise integral.  InternalNoiseSD is in the
% corresponding response units after the centered MT sensitivity is applied.

cfg = struct();
cfg.NValues = 2:20;
cfg.NTrials = 50000;
cfg.NCalibrationTrials = 30000;
cfg.ReferenceN = 6;
cfg.TargetPerformance = 0.75;
cfg.TuningSigmaDeg = 37.5;
cfg.CenterTuning = true;
cfg.ExternalNoiseSD = 1;
cfg.InternalNoiseSD = 1;
cfg.StreamPedestal = 0;
cfg.Seed = 73191;

cfg.CalibrationTolerance = 2e-3;
cfg.CalibrationMaxIterations = 35;
cfg.CalibrationInitialUpper = 1;

cfg.DiagnosticN = 8;
cfg.DiagnosticMode = 'Matched';
cfg.DiagnosticPerturbation = 0.10;
cfg.DiagnosticGapEdgesSD = -6:0.20:2;

cfg.PlotNValues = [2 3 4 6 8 10 12 16 20];
cfg.Visible = 'on';
cfg.SaveOutputs = true;
cfg.OutputDir = fullfile(pwd, 'MaxPoolingNSweep_Output');
end
