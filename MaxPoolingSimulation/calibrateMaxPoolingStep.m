function step = calibrateMaxPoolingStep(bank, targetPerformance, cfg)
% CALIBRATEMAXPOOLINGSTEP  Bisection calibration using common random numbers.

driftProfile = bank.W(1,:);
performanceAt = @(s) mean(max(bsxfun(@plus, bank.responseChangedBase, ...
    s*driftProfile), [], 2) > max(bank.responseNoChange, [], 2));

lower = 0;
upper = cfg.CalibrationInitialUpper;
while performanceAt(upper) < targetPerformance
    upper = 2*upper;
    if upper > 1e6
        error('Could not bracket the target performance.');
    end
end

for iter = 1:cfg.CalibrationMaxIterations
    mid = (lower + upper)/2;
    p = performanceAt(mid);
    if abs(p-targetPerformance) <= cfg.CalibrationTolerance
        step = mid;
        return;
    elseif p < targetPerformance
        lower = mid;
    else
        upper = mid;
    end
end
step = (lower + upper)/2;
end

