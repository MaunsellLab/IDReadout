function out = subsetSplitSignalTrialBank(bank,index)
% SUBSETSPLITSIGNALTRIALBANK  Select trials while preserving shared metadata.

out=bank; out.nTrials=numel(index);
fields={'externalChangedUnit','externalNoChangeUnit', ...
    'externalChangedFramesUnit','externalNoChangeFramesUnit', ...
    'internalChangedUnit','internalNoChangeUnit'};
for k=1:numel(fields)
    name=fields{k}; out.(name)=bank.(name)(index,:,:);
end
end
