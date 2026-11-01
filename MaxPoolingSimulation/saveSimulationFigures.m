function saveSimulationFigures(figures, outputDir)
% SAVESIMULATIONFIGURES  Save figures with exportgraphics/print compatibility.

names = {'SparseNSweep_GainCurves.pdf', ...
         'SparseNSweep_OperatingPoints.pdf', ...
         'SparseNSweep_GainMaps.pdf', ...
         'SparseNSweep_CloseToMax.pdf'};
for k = 1:min(numel(figures),numel(names))
    filename = fullfile(outputDir,names{k});
    try
        exportgraphics(figures(k),filename,'ContentType','vector');
    catch
        set(figures(k),'PaperPositionMode','auto');
        print(figures(k),filename,'-dpdf','-painters');
    end
end
end

