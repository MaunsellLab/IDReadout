function G = groupDiagnostics(cache,prediction)
%GROUPDIAGNOSTICS Observed choice rate, calibration, and likelihood by group.
T=cache.trials(prediction.trialIdx,:);
p1=prediction.pRDK1(prediction.trialIdx);
y=double(T.choiceRDK==1);
pChosen=y.*p1+(1-y).*(1-p1);
groupVars={'animal','domain'};
optional={'sessionID','probeOffsetDeg','baseCohPC','stepCohPC'};
groupVars=[groupVars optional(ismember(optional,T.Properties.VariableNames))];
[groupID,G]=findgroups(T(:,groupVars));
G.nTrials=splitapply(@numel,y,groupID);
G.observedPChooseRDK1=splitapply(@mean,y,groupID);
G.predictedPChooseRDK1=splitapply(@mean,p1,groupID);
G.calibrationResidual=G.observedPChooseRDK1-G.predictedPChooseRDK1;
G.logLikelihood=splitapply(@(x)sum(log(x)),pChosen,groupID);
end
