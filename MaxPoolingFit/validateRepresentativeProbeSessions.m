function report = validateRepresentativeProbeSessions(filePaths)
%VALIDATEREPRESENTATIVEPROBESESSIONS Audit adapters before full extraction.
%
% report = validateRepresentativeProbeSessions([idrFile idqFile idsFile])

if ischar(filePaths), filePaths={filePaths}; end
if isstring(filePaths), filePaths=cellstr(filePaths); end
assert(iscell(filePaths)&&~isempty(filePaths),'Supply representative MAT paths.');
addpath(fullfile(fileparts(mfilename('fullpath')),'adapters'));

domain=strings(numel(filePaths),1); nCandidate=zeros(numel(filePaths),1);
nValid=zeros(numel(filePaths),1); nIncrement=zeros(numel(filePaths),1);
nDecrement=zeros(numel(filePaths),1); nStepNoise=zeros(numel(filePaths),1);
noiseMin=nan(numel(filePaths),1); noiseMax=noiseMin; noiseMean=noiseMin;
meanPreDriftPC=noiseMin; meanPreProbePC=noiseMin;
minPreStreamPC=noiseMin; maxPreStreamPC=noiseMin;
trialIdxContract=strings(numel(filePaths),1);
for k=1:numel(filePaths)
  F=readMaxPoolingProbeSession(filePaths{k});
  F.schemaVersion="MaxPoolingTrialCache-v2";
  maxpoolfit.validateCache(F);
  A=F.extractionAudit; T=F.trials;
  domain(k)=A.domain; nCandidate(k)=A.nCandidateTrials; nValid(k)=A.nValidTrials;
  nIncrement(k)=sum(T.isIncrement); nDecrement(k)=sum(~T.isIncrement);
  nStepNoise(k)=sum(T.hasStepNoise); noiseMin(k)=A.rawNoiseMinPC;
  noiseMax(k)=A.rawNoiseMaxPC; noiseMean(k)=A.rawNoiseMeanPC;
  [meanPreDriftPC(k),meanPreProbePC(k), ...
    minPreStreamPC(k),maxPreStreamPC(k)]=summarizePrestep(F.preCohPC);
  trialIdxContract(k)=A.trialIdxContract;
end
report=table(string(filePaths(:)),domain,nCandidate,nValid,nIncrement,nDecrement, ...
  nStepNoise,noiseMin,noiseMax,noiseMean,meanPreDriftPC,meanPreProbePC, ...
  minPreStreamPC,maxPreStreamPC,trialIdxContract,'VariableNames', ...
  {'file','domain','nCandidate','nValid','nIncrement','nDecrement', ...
   'nStepNoise','rawNoiseMinPC','rawNoiseMaxPC','rawNoiseMeanPC', ...
   'meanPreDriftPC','meanPreProbePC','minPreStreamPC','maxPreStreamPC', ...
   'trialIdxContract'});
disp(report)
assert(numel(unique(domain))==numel(domain), ...
  'Representative files contain a duplicate or misidentified domain.');
assert(all(noiseMin<0 & noiseMax>0), ...
  'At least one stored trace does not span negative and positive noise values.');
end

function [driftMean,probeMean,minValue,maxValue]=summarizePrestep(preCohPC)
driftValues=[]; probeValues=[]; allValues=[];
for k=1:numel(preCohPC)
  x=double(preCohPC{k}(:));
  driftValues(end+1,1)=x(1); %#ok<AGROW>
  if numel(x)>1, probeValues=[probeValues;x(2:end)]; end %#ok<AGROW>
  allValues=[allValues;x]; %#ok<AGROW>
end
driftMean=mean(driftValues,'omitnan');
probeMean=mean(probeValues,'omitnan');
minValue=min(allValues,[],'omitnan');
maxValue=max(allValues,[],'omitnan');
end
