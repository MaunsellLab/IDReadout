function sharePath = exportMTReadoutPNormSummary(result,outputFolder)
%EXPORTMTREADOUTPNORMSUMMARY A small, explicit-list ZIP with no trial data.
% Accepts a result struct or a saved MTReadoutPNormResult.mat pathname.
if ischar(result)||isstring(result)
  A=load(result,'result'); result=A.result;
end
folder=fullfile(outputFolder,'Share'); if ~isfolder(folder), mkdir(folder); end
T=result.trials; mask=result.fitMask; offsets=unique(T.probeOffsetDeg);
rows=cell(numel(offsets),6);
for k=1:numel(offsets)
  q=T.probeOffsetDeg==offsets(k); m=q & mask;
  rows(k,:)={offsets(k),numel(unique(T.sessionID(q))),nnz(q),nnz(m), ...
    nnz(q & ~T.hasStepNoise),mean(T.correct(m))};
end
coverage=cell2table(rows,'VariableNames',{'probeOffsetDeg','nSessions', ...
  'nValidTrials','nFitTrials','nNoStepNoiseTrials','observedFitAccuracy'});
performance=table();
if ~isempty(result.fits)
  rows=cell(height(result.fits)*numel(offsets),7); n=0;
  for j=1:height(result.fits)
    pr=result.predictions{j}.pCorrect;
    p1=min(max(result.predictions{j}.pRDK1,1e-8),1-1e-8);
    y=double(T.choiceRDK==1);
    for k=1:numel(offsets)
      m=mask & T.probeOffsetDeg==offsets(k); n=n+1;
      ll=sum(y(m).*log(p1(m))+(1-y(m)).*log(1-p1(m)));
      rows(n,:)={result.fits.p(j),offsets(k),nnz(m),mean(T.correct(m)), ...
        mean(pr(m)),ll,result.fits.noiseSD(j)};
    end
  end
  performance=cell2table(rows,'VariableNames',{'p','probeOffsetDeg','nFitTrials', ...
    'observedAccuracy','predictedAccuracy','logLikelihood','noiseSD'});
end
summary=struct('config',result.config,'fits',result.fits,'coverage',coverage, ...
  'performanceByOffset',performance,'createdAt',datetime('now'));
summary.config.checkpointFile="";
if isfield(result,'noiseProfiles'), summary.noiseProfiles=result.noiseProfiles; end
if isfield(result,'responseAudit'), summary.responseAudit=result.responseAudit; end
summary.gainStatus="not supplied"; summary.gainNotes=strings(0,1);
if isfield(result,'gainStatus'), summary.gainStatus=result.gainStatus; end
if isfield(result,'gainNotes'), summary.gainNotes=result.gainNotes; end
if isfield(result,'selection')
  info=result.selection.fileInfo;
  names=intersect({'fileName','animal','sourceDomain','probeDirDeg', ...
    'parentNProbeDirections','nTrials','nNoiseTrials'},info.Properties.VariableNames,'stable');
  summary.selectedFiles=info(:,names);
  summary.subsetRun=result.selection.subsetRun;
end
files={'MTReadoutPNormSummary.mat','MTReadoutPNormCoverage.csv','README.txt'};
writetable(coverage,fullfile(folder,'MTReadoutPNormCoverage.csv'));
if ~isempty(result.fits)
  writetable(result.fits,fullfile(folder,'MTReadoutPNormFits.csv'));
  writetable(performance,fullfile(folder,'MTReadoutPNormPerformanceByOffset.csv'));
  files=[files {'MTReadoutPNormFits.csv','MTReadoutPNormPerformanceByOffset.csv'}];
end
if isfield(result,'gainRecovery') && strcmp(string(summary.gainStatus),"complete")
  summary.gains=gainTable(result.gainRecovery);
  writetable(summary.gains,fullfile(folder,'MTReadoutPNormGains.csv'));
  files{end+1}='MTReadoutPNormGains.csv';
  if isfield(result,'gainPlotAvailable') && result.gainPlotAvailable && ...
      isfile(fullfile(outputFolder,'MTReadoutPNormGainComparison.pdf'))
    copyfile(fullfile(outputFolder,'MTReadoutPNormGainComparison.pdf'), ...
      fullfile(folder,'MTReadoutPNormGainComparison.pdf'));
    files{end+1}='MTReadoutPNormGainComparison.pdf';
  end
end
% This struct intentionally omits trial tables, window inputs, individual
% predictions, raw traces, cache fragments, and the large source recordings.
save(fullfile(folder,'MTReadoutPNormSummary.mat'),'summary','-v7');
fid=fopen(fullfile(folder,'README.txt'),'w'); assert(fid>=0,'Cannot create share notes.');
closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'MTReadoutPNorm compact results\nModel: %s\n',result.config.modelVersion);
fprintf(fid,'Gain status: %s\n',summary.gainStatus);
fprintf(fid,'The ZIP contains aggregate summaries only; source and trial data stay local.\n');
fprintf(fid,'Likelihood calibration uses noisy increments. Gain predictors remain those of the existing IDS analysis.\n');
fprintf(fid,'A subset check or extraction-only run is not a full gain comparison.\n');
for k=1:numel(summary.gainNotes), fprintf(fid,'%s\n',summary.gainNotes(k)); end
clear closeFile
sharePath=fullfile(outputFolder,'MTReadoutPNormShare.zip');
zip(sharePath,files,folder);
end

function G=gainTable(recovered)
if isstruct(recovered), recovered=num2cell(recovered); end
rows=cell(0,11);
for k=1:numel(recovered)
  S=recovered{k}; fits=[{S.observed};S.predicted(:)];
  for j=1:numel(fits)
    F=fits{j}.driftProbe; N=fits{j}.normalized;
    if j==1
      source="observed"; p=NaN; ci=N.pooledCI95;
    else
      source="model"; p=S.p(j-1); ci=[NaN NaN];
    end
    rows(end+1,:)={S.probeOffsetDeg,source,p,S.nTrials,S.nSessions, ...
      F.gain(1),F.gain(2),N.pooledProbeOverDrift,ci(1),ci(2),string(S.noisePredictorType)}; %#ok<AGROW>
  end
end
G=cell2table(rows,'VariableNames',{'probeOffsetDeg','source','p','nTrials','nSessions', ...
  'driftGain','probeGain','probeOverDrift','observedRatioCILow','observedRatioCIHigh', ...
  'noisePredictorType'});
end
