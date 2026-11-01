function choice = simulateSplitSignalPoolingChoice(model,bank,modelName,noiseSD)
% SIMULATESPLITSIGNALPOOLINGCHOICE  One pre-pooling-noise choice realization.

changed=model.detChanged+single(noiseSD)*bank.internalChangedUnit;
noChange=model.detNoChange+single(noiseSD)*bank.internalNoChangeUnit;
changedValue=poolSplitSignalResponse(changed,bank.directionsDeg,modelName);
noChangeValue=poolSplitSignalResponse(noChange,bank.directionsDeg,modelName);
choice=changedValue>noChangeValue;
end
