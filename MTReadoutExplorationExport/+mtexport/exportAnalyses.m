function [relativeFiles,notes] = exportAnalyses(root,out,options)
%EXPORTANALYSES Preserve existing predictors/calibrations as plain values.
relativeFiles={}; notes=cell(0,3); directories={};
roots={root}; requested=cellstr(string(options.Domains));
for k=1:numel(requested)
  p=fullfile(root,requested{k}); if isfolder(p), roots{end+1}=p; end %#ok<AGROW>
end
for k=1:numel(roots)
  base=fullfile(roots{k},'Data'); if ~isfolder(base), continue; end
  D=dir(fullfile(base,'Probe*')); D=D([D.isdir]);
  for j=1:numel(D)
    p=fullfile(D(j).folder,D(j).name,'AcrossSessionSummaries');
    if isfolder(p), directories{end+1}=p; end %#ok<AGROW>
  end
  for name={'AcrossSessions','AcrossSessionSummaries','AcrossOffsetSummaries','AcrossSessionAnalysis'}
    p=fullfile(base,name{1}); if isfolder(p), directories{end+1}=p; end %#ok<AGROW>
  end
  D=dir(fullfile(base,'*Across*')); D=D([D.isdir]);
  for j=1:numel(D), directories{end+1}=fullfile(D(j).folder,D(j).name); end %#ok<AGROW>
end
directories=unique(directories,'stable'); target=fullfile(out,'analyses');
if ~isfolder(target), mkdir(target); end
for k=1:numel(directories)
  D=dir(fullfile(directories{k},'*.mat'));
  for j=1:numel(D)
    path=fullfile(D(j).folder,D(j).name);
    try
      names=who('-file',path);
      wanted=intersect(names,{'acrossSummary','summary','result','gainTable','metadata','sessionTable'},'stable');
      if isempty(wanted), continue; end
      A=load(path,wanted{:});
      % Remove recursive animal copies/large resampling outputs; retain
      % fitted curves, actual predictors, kernels, and psychometric tables.
      if isfield(A,'result')&&isstruct(A.result)
        remove=intersect(fieldnames(A.result),{'animalResults','bootstrapProfileLL','recovery'});
        if ~isempty(remove), A.result=rmfield(A.result,remove); end
      end
      analysis=mtexport.plain(A); %#ok<NASGU>
      analysis.export_source_name=D(j).name;
      analysis.export_note='Saved analysis preserved as-is; animal/task selection and calibration must be checked against source shards.';
      relative=fullfile('analyses',['analysis_' mtexport.id(path) '.mat']);
      save(fullfile(out,relative),'analysis','-v7');
      relativeFiles{end+1,1}=relative; %#ok<AGROW>
    catch ME
      notes(end+1,:)={'analysis',path,ME.message}; %#ok<AGROW>
    end
  end
end
codeFolder=fullfile(out,'analysis_code'); if ~isfolder(codeFolder), mkdir(codeFolder); end
functions={'selectAnalysisFiles','fitIDSNoiseGain','plotIDSSummary', ...
  'idqWeibullP','idqWeibullAlphaForThreshold','fitIDQInitialSessionThresholds', ...
  'fitIDQAcrossAlignedWeibull','fitIDQNoiseGain','plotIDQSummary', ...
  'fitIDJSessionPsychometrics','idjWeibullP','idjWeibullAlphaForThreshold', ...
  'plotIDJAcrossSessionBehavior','fitIDJModelRecoveryAcrossSessions','pilotIDJPNormProfile', ...
  'idjHasTopUpExclusion','normalizeIDQFullSessionForIDS'};
for k=1:numel(functions)
  path=which(functions{k});
  if isempty(path)||~endsWith(path,'.m'), continue; end
  relative=fullfile('analysis_code',[functions{k} '.m']);
  copyfile(path,fullfile(out,relative)); relativeFiles{end+1,1}=relative; %#ok<AGROW>
end
end
