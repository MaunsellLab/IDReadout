function result = runMTReadoutPNormLocal(varargin)
%RUNMTREADOUTPNORMLOCAL Folder selection, incremental extraction, local fit.
% First run: runMTReadoutPNormLocal  (one folder chooser; folder is remembered)
% Short first check: runMTReadoutPNormLocal('MaxFiles',2)
% Select another root: runMTReadoutPNormLocal('ChooseFolder',true)
p=inputParser;
p.addParameter('DataRoot','',@(x)ischar(x)||isstring(x));
p.addParameter('OutputFolder','',@(x)ischar(x)||isstring(x));
p.addParameter('ChooseFolder',false,@islogical);
p.addParameter('Config',mtpnorm.defaultConfig(),@isstruct);
p.addParameter('PGrid',[5 Inf],@(x)isnumeric(x)&&isvector(x)&&all(x>=1)&&all(~isnan(x)));
p.addParameter('MaxFiles',Inf,@(x)isscalar(x)&&x>=1&&(isinf(x)||x==fix(x)));
p.addParameter('ExtractOnly',false,@islogical);
p.addParameter('RebuildCache',false,@islogical);
p.addParameter('RecoverGains',true,@islogical);
p.addParameter('Summaries',[],@(x)isempty(x)||isstruct(x)||iscell(x));
p.addParameter('ReportExcluded',false,@islogical);
p.parse(varargin{:}); O=p.Results; cfg=O.Config;
root=char(O.DataRoot);
if isempty(root) && ~O.ChooseFolder && ispref('MTReadoutPNorm','dataRoot')
  root=getpref('MTReadoutPNorm','dataRoot');
end
if O.ChooseFolder || isempty(root) || ~isfolder(root)
  start=pwd; if isfolder(root), start=root; end
  chosen=uigetdir(start,'Choose IDS, IDS/Data, or a ProbeSessions folder');
  assert(~isequal(chosen,0),'Folder selection cancelled.'); root=chosen;
end
[files,info,folders]=mtpnorm.selectLocalFiles(root,cfg.animal,O.ReportExcluded);
setpref('MTReadoutPNorm','dataRoot',root);
subset=numel(files)>O.MaxFiles;
if subset, files=files(1:O.MaxFiles); info=info(1:O.MaxFiles,:); end
out=char(O.OutputFolder);
if isempty(out), out=fullfile(root,'MTReadoutPNormOutput'); end
if ~isfolder(out), mkdir(out); end
cfg.checkpointFile=string(fullfile(out,'MTReadoutPNormFitCheckpoint.mat'));
fprintf('Selected %d %s IDS files through selectAnalysisFiles.\n',numel(files),cfg.animal);
logFile=fullfile(out,'MTReadoutPNormRunLog.txt');
diary(logFile); cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
selection=struct('fileInfo',info,'sourceFolders',{folders},'subsetRun',subset, ...
  'dataRoot',string(root),'selector',string(which('selectAnalysisFiles')));
save(fullfile(out,'MTReadoutPNormSelection.mat'),'selection');
safeNames=intersect({'fileName','animal','sourceDomain','probeDirDeg', ...
  'parentNProbeDirections','nTrials','nNoiseTrials'},info.Properties.VariableNames,'stable');
writetable(info(:,safeNames),fullfile(out,'MTReadoutPNormSelection.csv'));
[cache,manifest]=mtpnorm.prepareLocalCache(files,fullfile(out,'InputCache'),cfg,O.RebuildCache);
save(fullfile(out,'MTReadoutPNormTrialCache.mat'),'cache','-v7.3');
writetable(manifest,fullfile(out,'MTReadoutPNormCacheManifest.csv'));
writetable(cache.sourceAudit,fullfile(out,'MTReadoutPNormInputAudit.csv'));
if O.ExtractOnly
  result=struct('config',cfg,'trials',cache.trials,'fitMask', ...
    cache.trials.hasStepNoise & cache.trials.isIncrement,'fits',table(), ...
    'gainStatus',"extraction only",'gainNotes',strings(0,1),'selection',selection);
else
  result=mtpnorm.fit(cache,cfg,O.PGrid);
  writetable(result.fits,fullfile(out,'MTReadoutPNormFits.csv'));
  result.selection=selection; result.gainStatus="not requested"; result.gainNotes=strings(0,1);
  result.gainPlotAvailable=false;
  % Save the completed fit before optional downstream analysis.
  save(fullfile(out,'MTReadoutPNormResult.mat'),'result','-v7.3');
  if O.RecoverGains && ~subset
    try
      summaries=O.Summaries;
      if isempty(summaries)
        [summaries,result.gainNotes]=mtpnorm.findLocalSummaries(folders,cfg.animal);
      end
      if isempty(summaries)
        result.gainStatus="matching saved IDS summaries not found";
      else
        result.gainRecovery=recoverMTReadoutPNormGains(result,summaries);
        result.gainStatus="complete";
      end
    catch ME
      result.gainStatus="gain recovery failed; fit is saved";
      result.gainNotes(end+1,1)=string(getReport(ME,'extended','hyperlinks','off'));
      warning('MTReadoutPNorm:GainRecovery','%s',ME.message);
    end
    if strcmp(string(result.gainStatus),"complete")
      try
        plotMTReadoutPNormGains(result.gainRecovery,fullfile(out,'MTReadoutPNormGainComparison.pdf'));
        result.gainPlotAvailable=true;
      catch ME
        result.gainNotes(end+1,1)="Gain values recovered; plot failed: "+string(ME.message);
        warning('MTReadoutPNorm:GainPlot','%s',ME.message);
      end
    end
  elseif subset && O.RecoverGains
    result.gainStatus="subset check; full-summary gain recovery skipped";
  end
end
save(fullfile(out,'MTReadoutPNormResult.mat'),'result','-v7.3');
sharePath=exportMTReadoutPNormSummary(result,out);
fprintf('\nLocal outputs: %s\nShare this compact file: %s\n',out,sharePath);
if ~isempty(result.gainNotes), disp(result.gainNotes); end
end
