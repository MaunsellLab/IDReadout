function validateSplitSignalSimulation()
% VALIDATESPLITSIGNALSIMULATION  Fast structural and qualitative checks.

cfg=splitSignalDefaults();
cfg.SeparationsDeg=90; cfg.AllocationFractions=[0.5 0.75];
cfg.TotalStepCoherences=[16 28]; cfg.StepNoiseHalfAmplitudes=7.071;
cfg.PrestepPedestalFractions=1; cfg.IncludeBalancedPrestep=[false true];
cfg.AdaptationFractions=0; cfg.NDesignTrials=2500;
cfg.NPowerReplicates=10; cfg.MaxSessions=2; cfg.TrialsPerSession=100;
cfg.CreateFigures=false; cfg.SaveOutputs=false;
pairs=structToPairs(cfg); summary=runSplitSignalDesignSweep(pairs{:});
assert(~isempty(summary.conditions));
assert(all(isfinite([summary.conditions.performance])));
assert(all([summary.conditions.performance]>=0 & [summary.conditions.performance]<=1));

bank=makeSplitSignalTrialBank(200,cfg,cfg.Seed+1);
m=buildSplitSignalDesign(bank,90,0.75,20,7.071,1,false,0,cfg);
assert(all(size(m.noiseDifference)==[200 2]));
assert(all(abs(double(m.meanStreams)-[15 5])<1e-6));
mb=buildSplitSignalDesign(bank,90,0.75,20,7.071,1,true,0,cfg);
assert(abs(mb.balanceDirection-225)<1e-6);
assert(abs(mb.balanceCoherence-sqrt(2)*7.071)<1e-3);
fprintf('validateSplitSignalSimulation passed.\n');
end

function pairs=structToPairs(s)
names=fieldnames(s); pairs=cell(1,2*numel(names));
for k=1:numel(names), pairs{2*k-1}=names{k}; pairs{2*k}=s.(names{k}); end
end
