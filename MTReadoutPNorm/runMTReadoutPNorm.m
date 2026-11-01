function result = runMTReadoutPNorm(sourceFiles,outputFolder,varargin)
%RUNMTREADOUTPNORM Fit the agreed two-window MT forward model on IDS data.
% result=runMTReadoutPNorm(files,out,'PGrid',[5 Inf]);
% sourceFiles: corrected raw IDS *.mat files OR their ProbeSessions, not both.
% Set Config to mtpnorm.defaultConfig() with explicit overrides if needed.
p=inputParser;
p.addParameter('PGrid',[5 Inf],@(x)isnumeric(x)&&isvector(x)&&all(x>=1)&&all(~isnan(x)));
p.addParameter('Config',mtpnorm.defaultConfig(),@isstruct);
p.addParameter('Summaries',[],@(x)isempty(x)||isstruct(x)||iscell(x));
p.parse(varargin{:}); cfg=p.Results.Config;
if ~isfolder(outputFolder), mkdir(outputFolder); end
cache=mtpnorm.readIDS(sourceFiles,cfg);
result=mtpnorm.fit(cache,cfg,p.Results.PGrid);
save(fullfile(outputFolder,'MTReadoutPNormTrialCache.mat'),'cache','-v7.3');
writetable(cache.sourceAudit,fullfile(outputFolder,'MTReadoutPNormInputAudit.csv'));
writetable(result.fits,fullfile(outputFolder,'MTReadoutPNormFits.csv'));
if ~isempty(p.Results.Summaries)
  result.gainRecovery=recoverMTReadoutPNormGains(result,p.Results.Summaries);
  plotMTReadoutPNormGains(result.gainRecovery, ...
    fullfile(outputFolder,'MTReadoutPNormGainComparison.pdf'));
end
save(fullfile(outputFolder,'MTReadoutPNormResult.mat'),'result','-v7.3');
fprintf('Saved %d trials from %d sessions.\n',height(cache.trials), ...
  numel(unique(cache.trials.sessionID)));
disp(result.fits);
end
