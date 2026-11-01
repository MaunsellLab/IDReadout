function [result,cacheReport] = updateAndRunMaxPoolingAnalysis( ...
    projectRoot,cachePath,outputFolder,varargin)
%UPDATEANDRUNMAXPOOLINGANALYSIS Synchronize the trial cache, then refit.
%
% [result,report] = updateAndRunMaxPoolingAnalysis(root,cachePath,out)
% Optional name/value arguments:
%   'SyncOptions' - cell array forwarded to syncExperimentalMaxPoolingCache
%   'FitOptions'  - cell array forwarded to runMaxPoolingFit

% A synchronization failure leaves the old cache untouched and prevents the
% fit, so a saved result cannot accidentally use a partly updated cache.

p=inputParser;
p.addParameter('SyncOptions',{},@iscell);
p.addParameter('FitOptions',{},@iscell);
p.parse(varargin{:});

[~,cacheReport]=syncExperimentalMaxPoolingCache(projectRoot,cachePath, ...
  p.Results.SyncOptions{:});
result=runMaxPoolingFit(cachePath,outputFolder,p.Results.FitOptions{:});
end
