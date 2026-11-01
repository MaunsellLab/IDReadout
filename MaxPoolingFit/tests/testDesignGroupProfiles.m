function tests=testDesignGroupProfiles
tests=functiontests(localfunctions);
end

function testFourGroupsAndNestedComparisons(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',2);
T=cache.trials;
idq=string(T.domain)=="IDQ";
T.driftDirectionDeg(idq)=30;
idqSessions=unique(string(T.sessionID(idq)),'stable');
q=find(string(T.sessionID)==idqSessions(1));
T.driftDirectionDeg(q)=[10;130];
cache.trials=T;
[D,G]=maxpoolfit.driftRegimeMetadata(cache,'ReadHeaders',false);
verifyEqual(testCase,sort(unique(G)), ...
  sort(["IDR";"IDS";"IDQ-randomized";"IDQ-fixed"]));
verifyTrue(testCase,all(D.analysisGroup(D.domain=="IDQ" & ...
  D.realizedRandomBaseDir)=="IDQ-randomized"));
verifyTrue(testCase,all(D.analysisGroup(D.domain=="IDQ" & ...
  ~D.realizedRandomBaseDir)=="IDQ-fixed"));

cfg=maxpoolfit.defaultConfig(); cfg.verbose=false; cfg.batchSize=8;
R=maxpoolfit.profileDesignGroups(cache,cfg,G,'NDraws',8, ...
  'AdaptationGrid',[0 0.5]);
verifyEqual(testCase,R.groupSummary.analysisGroup,R.groupOrder);
verifyEqual(testCase,R.groupOrder, ...
  ["IDQ-fixed";"IDQ-randomized";"IDS";"IDR"]);
verifyEqual(testCase,height(R.groupSummary),4);
verifyEqual(testCase,height(R.groupProfile),8);
verifyEqual(testCase,height(R.modelComparison),5);
verifyTrue(testCase,all(isfinite(R.groupSummary.logLikelihood)));
verifyTrue(testCase,all(R.modelComparison.estimable));
verifyEqual(testCase,R.modelComparison.deltaLogLFromPooled([1 3]),[0;0], ...
  'AbsTol',1e-12);
end

function testMissingHeaderFallsBackToRealizedDirections(testCase)
% Use several trials so the selected synthetic session can contain more
% than one realized drift direction.
cache=maxpoolfit.syntheticCache('TrialsPerCondition',3);
q=find(string(cache.trials.domain)=="IDQ");
session=string(cache.trials.sessionID(q(1)));
q=find(string(cache.trials.sessionID)==session);
cache.trials.driftDirectionDeg(q)=[10;130;250];
[D,G]=maxpoolfit.driftRegimeMetadata(cache,'ReadHeaders',true);
row=D.sessionID==session;
verifyFalse(testCase,D.headerAvailable(row));
verifyTrue(testCase,D.realizedRandomBaseDir(row));
verifyEqual(testCase,D.analysisGroup(row),"IDQ-randomized");
verifyTrue(testCase,all(G(q)=="IDQ-randomized"));
end

function testFixedTripletWithMultipleDriftsRemainsFixed(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',3);
q=find(string(cache.trials.domain)=="IDQ");
session=string(cache.trials.sessionID(q(1)));
q=find(string(cache.trials.sessionID)==session);
cache.trials.driftDirectionDeg(q)=[30;150;270];
[D,G]=maxpoolfit.driftRegimeMetadata(cache,'ReadHeaders',false);
row=D.sessionID==session;
verifyEqual(testCase,D.nUniqueDriftDirections(row),3);
verifyFalse(testCase,D.realizedRandomBaseDir(row));
verifyEqual(testCase,D.analysisGroup(row),"IDQ-fixed");
verifyTrue(testCase,all(G(q)=="IDQ-fixed"));
end

function testHeaderIsAuthoritativeAndAudited(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
cache.trials.sourceFile=strings(height(cache.trials),1);
cache.trials.driftDirectionDeg(string(cache.trials.domain)=="IDQ")=30;
folder=tempname; mkdir(folder);
cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
idqSessions=unique(string(cache.trials.sessionID( ...
  string(cache.trials.domain)=="IDQ")),'stable');
for k=1:numel(idqSessions)
  path=fullfile(folder,sprintf('idq_%02d.mat',k));
  sessionHeader=struct('randomBaseDir',k==1); %#ok<NASGU>
  save(path,'sessionHeader');
  cache.trials.sourceFile(string(cache.trials.sessionID)==idqSessions(k))=path;
end
warningState=warning('off', ...
  'maxpoolfit:driftRegimeMetadata:HeaderDirectionMismatch');
warningCleanup=onCleanup(@()warning(warningState)); %#ok<NASGU>
[D,G]=maxpoolfit.driftRegimeMetadata(cache,'ReadHeaders',true);
row=D.sessionID==idqSessions(1);
verifyTrue(testCase,D.headerAvailable(row));
verifyEqual(testCase,D.randomBaseDir(row),1);
verifyTrue(testCase,D.auditMismatch(row));
verifyEqual(testCase,D.analysisGroup(row),"IDQ-randomized");
verifyTrue(testCase,all(G(string(cache.trials.sessionID)==idqSessions(1))== ...
  "IDQ-randomized"));
end

function testMissingIDQFixedIsReportedNotRejected(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',2);
idqSessions=unique(string(cache.trials.sessionID( ...
  string(cache.trials.domain)=="IDQ")),'stable');
for k=1:numel(idqSessions)
  q=find(string(cache.trials.sessionID)==idqSessions(k));
  cache.trials.driftDirectionDeg(q)=[10;130];
end
[~,G]=maxpoolfit.driftRegimeMetadata(cache,'ReadHeaders',false);
warningState=warning('off','maxpoolfit:profileDesignGroups:MissingGroups');
warningCleanup=onCleanup(@()warning(warningState)); %#ok<NASGU>
cfg=maxpoolfit.defaultConfig(); cfg.verbose=false; cfg.batchSize=8;
R=maxpoolfit.profileDesignGroups(cache,cfg,G,'NDraws',8, ...
  'AdaptationGrid',[0 0.5]);
verifyFalse(testCase,any(R.groupOrder=="IDQ-fixed"));
verifyFalse(testCase,any(R.modelComparison.estimable(1:2)));
verifyTrue(testCase,all(R.modelComparison.estimable(3:5)));
verifyTrue(testCase,all(isnan(R.modelComparison.logLikelihood(1:2))));
verifyTrue(testCase,all(contains(R.modelComparison.note(1:2),"IDQ-fixed")));
end

function testSelectedGroupedFitContract(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',2);
idq=string(cache.trials.domain)=="IDQ";
cache.trials.driftDirectionDeg(idq)=30;
idqSessions=unique(string(cache.trials.sessionID(idq)),'stable');
q=find(string(cache.trials.sessionID)==idqSessions(1));
cache.trials.driftDirectionDeg(q)=[10;130];
[~,trialGroup]=maxpoolfit.driftRegimeMetadata(cache,'ReadHeaders',false);
adaptationByGroup=table( ...
  ["IDQ-fixed";"IDQ-randomized";"IDS";"IDR"], ...
  [0.1;0.4;0.4;0.6], ...
  'VariableNames',{'analysisGroup','adaptation'});
cfg=maxpoolfit.defaultConfig(); cfg.verbose=false; cfg.batchSize=8;
fit=maxpoolfit.fitDesignGroupModel(cache,cfg,trialGroup, ...
  adaptationByGroup,'NDraws',8);
prediction=maxpoolfit.predictDesignGroupModel(cache,fit,cfg,'NDraws',8);
verifyEqual(testCase,height(fit.groupFit),4);
verifyEqual(testCase,fit.adaptationByGroup,adaptationByGroup);
verifyEqual(testCase,fit.groupFit.adaptation( ...
  fit.groupFit.analysisGroup=="IDQ-randomized"),0.4);
verifyEqual(testCase,fit.groupFit.adaptation( ...
  fit.groupFit.analysisGroup=="IDS"),0.4);
verifyTrue(testCase,all(isfinite( ...
  prediction.pRDK1(prediction.trialIdx))));
verifyEqual(testCase,prediction.trialIdx,find(fit.trialMask));
end
