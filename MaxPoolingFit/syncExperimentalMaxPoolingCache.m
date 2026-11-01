function [cache,report,scanManifest] = syncExperimentalMaxPoolingCache(projectRoot,cachePath,varargin)
%SYNCEXPERIMENTALMAXPOOLINGCACHE Standard incremental experimental cache run.
%
% First run extracts every discovered ProbeSession. Later runs open only new
% or modified ProbeSessions, update their cached trials, recompute LOSO
% predictors from cached framewise values, and atomically replace the cache.

thisFolder=fileparts(mfilename('fullpath'));
addpath(fullfile(thisFolder,'adapters'));
p=inputParser;
p.KeepUnmatched=true;
p.addParameter('Domains',["IDR" "IDQ" "IDS"]);
p.addParameter('Animals',["Meetz" "Neesha"]);
p.addParameter('FileSelectionArgs',{},@iscell);
p.parse(varargin{:});
scanManifest=maxpoolfit.discoverProbeSessions(projectRoot, ...
  'Domains',p.Results.Domains,'Animals',p.Results.Animals, ...
  'FileSelectionArgs',p.Results.FileSelectionArgs);
args=unmatchedToCell(p.Unmatched);
[cache,report]=maxpoolfit.syncCache(scanManifest,@readMaxPoolingProbeSession, ...
  cachePath,args{:});
end

function args=unmatchedToCell(S)
names=fieldnames(S); args=cell(1,2*numel(names));
for k=1:numel(names), args{2*k-1}=names{k}; args{2*k}=S.(names{k}); end
end
