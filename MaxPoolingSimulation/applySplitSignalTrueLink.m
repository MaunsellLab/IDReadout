function probability = applySplitSignalTrueLink(link,scores,cfg)
% APPLYSPLITSIGNALTRUELINK  Apply an independently calibrated smooth link.

if strcmpi(link.modelName,'max')
    [~,index]=min(abs(cfg.SoftmaxTemperatures-link.bestTemperature));
    predictor=scores.max(:,index);
else
    predictor=scores.(link.modelName);
end
eta=link.fit.intercept+double(predictor)*link.fit.beta;
probability=1./(1+exp(-max(min(eta,35),-35)));
probability=cfg.TrueLapse/2+(1-cfg.TrueLapse)*probability;
end
