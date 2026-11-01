function scores = computeSplitSignalScoreSet(design,cfg)
% COMPUTESPLITSIGNALSCORESET  Smooth candidate summaries of deterministic MT.

scores=struct();
scores.max=zeros(design.nTrials,numel(cfg.SoftmaxTemperatures));
for k=1:numel(cfg.SoftmaxTemperatures)
    temperature=cfg.SoftmaxTemperatures(k);
    scores.max(:,k)=splitSignalSoftMax(design.detChanged,temperature)- ...
        splitSignalSoftMax(design.detNoChange,temperature);
end
scores.vector=poolSplitSignalResponse(design.detChanged,design.directionsDeg,'vector')- ...
    poolSplitSignalResponse(design.detNoChange,design.directionsDeg,'vector');
scores.energy=poolSplitSignalResponse(design.detChanged,design.directionsDeg,'energy')- ...
    poolSplitSignalResponse(design.detNoChange,design.directionsDeg,'energy');
end
