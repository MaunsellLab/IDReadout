function savePairedProbeCorrelationSweepFigures(figures,outputDir)
% SAVEPAIREDPROBECORRELATIONSWEEPFIGURES  Export correlation-sweep figures.

names={'PairedProbe_CorrelationSweep_Gains.pdf', ...
    'PairedProbe_CorrelationSweep_OperatingPoints.pdf', ...
    'PairedProbe_CorrelationSweep_NoiseCalibration.pdf', ...
    'PairedProbe_CorrelationSweep_SelectedOffsets.pdf', ...
    'PairedProbe_CorrelationSweep_Winners.pdf'};
for k=1:min(numel(figures),numel(names))
    filename=fullfile(outputDir,names{k});
    try
        exportgraphics(figures(k),filename,'ContentType','vector');
    catch
        set(figures(k),'PaperPositionMode','auto');
        print(figures(k),filename,'-dpdf','-painters');
    end
end
end
