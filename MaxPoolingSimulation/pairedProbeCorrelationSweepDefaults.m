function cfg = pairedProbeCorrelationSweepDefaults()
% PAIREDPROBECORRELATIONSWEEPDEFAULTS  Defaults for the correlation sweep.

cfg=struct();
cfg.CorrelationLengthsDeg=[10 15 20 30 45 60 90];

paired=pairedProbeDefaults();
paired.OffsetsDeg=[10 25 45 60 75 90 120 135 170 180];
paired.NTrials=12000;
paired.NCalibrationTrials=8000;
paired.Visible='off';
paired.CreateFigures=false;
paired.SaveOutputs=false;
cfg.PairedConfig=paired;

cfg.Visible='on';
cfg.SaveOutputs=true;
cfg.OutputDir=fullfile(pwd,'PairedProbeCorrelationSweep_Output');
end
