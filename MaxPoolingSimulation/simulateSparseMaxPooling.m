function sim = simulateSparseMaxPooling(bank, stepCoherence, cfg)
% SIMULATESPARSEMAXPOOLING  Generate choices and recover probit stream gains.

responseChanged = bsxfun(@plus, bank.responseChangedBase, ...
    stepCoherence*bank.W(1,:));
[valueChanged, winnerChanged] = max(responseChanged, [], 2);
[valueNoChange, winnerNoChange] = max(bank.responseNoChange, [], 2);
correct = valueChanged > valueNoChange;

fit = fitProbitIRLS(bank.predictors, correct);
winProbabilityChanged = accumarray(winnerChanged, 1, [bank.n 1]) / bank.nTrials;
winProbabilityNoChange = accumarray(winnerNoChange, 1, [bank.n 1]) / bank.nTrials;

sim = struct();
sim.n = bank.n;
sim.directionsDeg = bank.directionsDeg;
sim.offsetDeg = bank.offsetDeg;
sim.stepCoherence = stepCoherence;
sim.performance = mean(correct);
sim.correct = correct;
sim.fit = fit;
sim.winProbabilityChanged = winProbabilityChanged;
sim.winProbabilityNoChange = winProbabilityNoChange;
sim.responseChanged = responseChanged;
sim.valueChanged = valueChanged;
sim.winnerChanged = winnerChanged;
sim.config = cfg;
end

