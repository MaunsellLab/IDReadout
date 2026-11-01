function savePairedProbeFigures(figures,outputDir)
% SAVEPAIREDPROBEFIGURES  Export paired-probe figures without TeX labels.

names={'PairedProbe_Gains.pdf','PairedProbe_Winners.pdf', ...
    'PairedProbe_OperatingPoints.pdf','PairedProbe_Decomposition.pdf'};
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

