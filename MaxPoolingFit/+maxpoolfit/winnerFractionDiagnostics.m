function W = winnerFractionDiagnostics(cache,prediction,cfg,trialGroup)
%WINNERFRACTIONDIAGNOSTICS Summarize within-RDK max winners.
% Winners are assigned to the nearest presented stream direction (a
% circular Voronoi assignment), then classified as drift or probe.

trialGroup=string(trialGroup(:));
assert(numel(trialGroup)==height(cache.trials), ...
  'trialGroup must contain one label per cached trial.');
rows=cell(0,8);
for g=1:numel(prediction.winnerDetails)
  item=prediction.winnerDetails{g};
  if isempty(item)||isempty(item.trialIdx), continue; end
  idx=item.trialIdx(:); D=item.details;
  for k=1:numel(idx)
    t=idx(k); changed=double(cache.trials.changeRDK(t));
    unchanged=3-changed;
    if changed==1
      wiChanged=D.winnerMechanismIndexRDK1(k,:);
      wiUnchanged=D.winnerMechanismIndexRDK2(k,:);
    else
      wiChanged=D.winnerMechanismIndexRDK2(k,:);
      wiUnchanged=D.winnerMechanismIndexRDK1(k,:);
    end
    changedProbe=probeWinnerFraction(cache,t,changed,wiChanged,cfg);
    unchangedProbe=probeWinnerFraction(cache,t,unchanged,wiUnchanged,cfg);
    rows(end+1,:)={t,trialGroup(t),string(cache.trials.animal(t)), ... %#ok<AGROW>
      double(cache.trials.probeOffsetDeg(t)), ...
      changedProbe,1-changedProbe,unchangedProbe,1-unchangedProbe};
  end
end
trialTable=cell2table(rows,'VariableNames',{'trialIdx','analysisGroup', ...
  'animal','probeOffsetDeg','changedProbeWinnerFraction', ...
  'changedDriftWinnerFraction','unchangedProbeWinnerFraction', ...
  'unchangedDriftWinnerFraction'});
if isempty(trialTable), W=trialTable; return; end
trialTable.analysisGroup=string(trialTable.analysisGroup);
trialTable.animal=string(trialTable.animal);
[G,analysisGroup,animal,probeOffsetDeg]=findgroups( ...
  trialTable.analysisGroup,trialTable.animal,trialTable.probeOffsetDeg);
nTrials=splitapply(@numel,trialTable.trialIdx,G);
nSessions=splitapply(@(x)numel(unique(string(cache.trials.sessionID(x)))), ...
  trialTable.trialIdx,G);
changedProbeWinnerFraction=splitapply(@mean, ...
  trialTable.changedProbeWinnerFraction,G);
changedDriftWinnerFraction=splitapply(@mean, ...
  trialTable.changedDriftWinnerFraction,G);
unchangedProbeWinnerFraction=splitapply(@mean, ...
  trialTable.unchangedProbeWinnerFraction,G);
unchangedDriftWinnerFraction=splitapply(@mean, ...
  trialTable.unchangedDriftWinnerFraction,G);
W=table(analysisGroup,animal,probeOffsetDeg,nTrials,nSessions, ...
  changedProbeWinnerFraction,changedDriftWinnerFraction, ...
  unchangedProbeWinnerFraction,unchangedDriftWinnerFraction);
W.winnerAssignment=repmat("nearest presented stream direction",height(W),1);
end

function f=probeWinnerFraction(cache,t,rdk,winnerIndex,cfg)
winnerDir=double(cfg.directionGridDeg(double(winnerIndex)));
drift=double(cache.driftDirDeg(t,rdk));
streams=double(cache.streamDirDeg{t,rdk}(:));
probe=streams(circularDistance(streams,drift)>1e-6);
if isempty(probe), f=0; return; end
driftDistance=circularDistance(winnerDir,drift);
probeDistance=inf(size(winnerDir));
for j=1:numel(probe)
  probeDistance=min(probeDistance,circularDistance(winnerDir,probe(j)));
end
f=mean(probeDistance<driftDistance);
end

function d=circularDistance(a,b)
d=abs(mod(double(a)-double(b)+180,360)-180);
end
