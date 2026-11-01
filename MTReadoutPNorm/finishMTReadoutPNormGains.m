function result = finishMTReadoutPNormGains(varargin)
%FINISHMTREADOUTPNORMGAINS Recover gains/export without refitting the model.
% Run plotIDSSummary normally, then finishMTReadoutPNormGains('Summaries',S).
% With no arguments, choose the existing result MAT file in a file picker.
p=inputParser;
p.addParameter('ResultFile','',@(x)ischar(x)||isstring(x));
p.addParameter('Summaries',[],@(x)isempty(x)||isstruct(x)||iscell(x));
p.parse(varargin{:}); path=char(p.Results.ResultFile);
if isempty(path)
  start=pwd;
  if ispref('MTReadoutPNorm','dataRoot')
    start=fullfile(getpref('MTReadoutPNorm','dataRoot'),'MTReadoutPNormOutput');
  end
  [name,folder]=uigetfile('MTReadoutPNormResult.mat','Choose completed forward-model result',start);
  assert(~isequal(name,0),'Selection cancelled.'); path=fullfile(folder,name);
end
A=load(path,'result'); result=A.result; out=fileparts(path);
assert(~isempty(result.fits),'This is an extraction-only result; fit the model first.');
summaries=p.Results.Summaries;
if isempty(summaries)
  assert(isfield(result,'selection'),'Supply the current IDS summary structs.');
  [summaries,notes]=mtpnorm.findLocalSummaries(result.selection.sourceFolders,result.config.animal);
else
  notes=strings(0,1);
end
assert(~isempty(summaries),'Create Meetz-only IDS summaries, then rerun this function.');
result.gainRecovery=recoverMTReadoutPNormGains(result,summaries);
result.gainStatus="complete"; result.gainNotes=notes;
result.gainPlotAvailable=false;
try
  plotMTReadoutPNormGains(result.gainRecovery,fullfile(out,'MTReadoutPNormGainComparison.pdf'));
  result.gainPlotAvailable=true;
catch ME
  result.gainNotes(end+1,1)="Gain values recovered; plot failed: "+string(ME.message);
  warning('MTReadoutPNorm:GainPlot','%s',ME.message);
end
save(path,'result','-v7.3');
sharePath=exportMTReadoutPNormSummary(result,out);
fprintf('Share this compact file: %s\n',sharePath);
end
