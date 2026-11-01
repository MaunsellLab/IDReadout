function calibration = calibrationDiagnostics(cache,prediction,varargin)
%CALIBRATIONDIAGNOSTICS Calibration summaries for fitted trial predictions.
p=inputParser;
p.addParameter('NBins',10,@(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
p.parse(varargin{:});
nBins=p.Results.NBins;

idx=prediction.trialIdx(:);
assert(~isempty(idx),'Prediction contains no trials.');
assert(all(idx>=1 & idx<=height(cache.trials)),'Prediction trial indices are invalid.');
T=cache.trials(idx,:);
p1=double(prediction.pRDK1(idx));
assert(all(isfinite(p1)) && all(p1>0 & p1<1), ...
  'All selected predictions must be finite probabilities strictly between zero and one.');
y=double(T.choiceRDK==1);
observedCorrect=double(T.choiceRDK==T.changeRDK);
change1=double(T.changeRDK==1);
predictedCorrect=change1.*p1+(1-change1).*(1-p1);

groupVars={'domain','animal'};
[groupID,groupKeys]=findgroups(T(:,groupVars));
groupSummary=summarizeGroups(groupKeys,groupID,y,p1,observedCorrect,predictedCorrect,true);

reliability=groupKeys([],:);
reliability.bin=zeros(0,1);
reliability.nTrials=zeros(0,1);
reliability.meanPredictedPChooseRDK1=zeros(0,1);
reliability.observedPChooseRDK1=zeros(0,1);
reliability.calibrationResidual=zeros(0,1);
reliability.ciLow=zeros(0,1);
reliability.ciHigh=zeros(0,1);
groupSummary.expectedCalibrationError=nan(height(groupSummary),1);
for g=1:height(groupKeys)
  q=find(groupID==g);
  [~,order]=sort(p1(q));
  bin=zeros(numel(q),1);
  bin(order)=min(nBins,ceil((1:numel(q))'.*nBins./numel(q)));
  for b=1:nBins
    qb=q(bin==b);
    if isempty(qb), continue; end
    observed=mean(y(qb));
    se=sqrt(max(observed*(1-observed),0)/numel(qb));
    row=groupKeys(g,:);
    row.bin=b;
    row.nTrials=numel(qb);
    row.meanPredictedPChooseRDK1=mean(p1(qb));
    row.observedPChooseRDK1=observed;
    row.calibrationResidual=observed-row.meanPredictedPChooseRDK1;
    row.ciLow=max(0,observed-1.96*se);
    row.ciHigh=min(1,observed+1.96*se);
    reliability=[reliability;row]; %#ok<AGROW>
  end
  R=reliability(string(reliability.domain)==string(groupKeys.domain(g)) & ...
    string(reliability.animal)==string(groupKeys.animal(g)),:);
  groupSummary.expectedCalibrationError(g)=sum( ...
    R.nTrials.*abs(R.calibrationResidual))/sum(R.nTrials);
end

sessionVars={'domain','animal','sessionID'};
[sessionID,sessionKeys]=findgroups(T(:,sessionVars));
sessionSummary=summarizeGroups(sessionKeys,sessionID,y,p1,observedCorrect,predictedCorrect,false);

conditionVars={'domain','animal'};
optional={'probeOffsetDeg','stepMagnitudePC'};
conditionVars=[conditionVars optional(ismember(optional,T.Properties.VariableNames))];
[conditionID,conditionKeys]=findgroups(T(:,conditionVars));
conditionSummary=summarizeGroups(conditionKeys,conditionID,y,p1,observedCorrect,predictedCorrect,false);

calibration=struct();
calibration.groupSummary=groupSummary;
calibration.reliability=reliability;
calibration.sessionSummary=sessionSummary;
calibration.conditionSummary=conditionSummary;
calibration.nTrials=numel(idx);
calibration.nBins=nBins;
calibration.trialIdx=idx;
calibration.createdAt=datetime('now');
end

function S=summarizeGroups(keys,groupID,y,p1,observedCorrect,predictedCorrect,fitCalibration)
n=height(keys);
S=keys;
S.nTrials=zeros(n,1);
S.observedPChooseRDK1=nan(n,1);
S.predictedPChooseRDK1=nan(n,1);
S.choiceResidual=nan(n,1);
S.observedAccuracy=nan(n,1);
S.predictedAccuracy=nan(n,1);
S.accuracyResidual=nan(n,1);
S.logLoss=nan(n,1);
S.brierScore=nan(n,1);
S.calibrationIntercept=nan(n,1);
S.calibrationSlope=nan(n,1);
for g=1:n
  q=groupID==g;
  S.nTrials(g)=sum(q);
  S.observedPChooseRDK1(g)=mean(y(q));
  S.predictedPChooseRDK1(g)=mean(p1(q));
  S.choiceResidual(g)=S.observedPChooseRDK1(g)-S.predictedPChooseRDK1(g);
  S.observedAccuracy(g)=mean(observedCorrect(q));
  S.predictedAccuracy(g)=mean(predictedCorrect(q));
  S.accuracyResidual(g)=S.observedAccuracy(g)-S.predictedAccuracy(g);
  pChosen=y(q).*p1(q)+(1-y(q)).*(1-p1(q));
  S.logLoss(g)=-mean(log(pChosen));
  S.brierScore(g)=mean((y(q)-p1(q)).^2);
  if fitCalibration
    [S.calibrationIntercept(g),S.calibrationSlope(g)]=calibrationCoefficients(y(q),p1(q));
  end
end
end

function [intercept,slope]=calibrationCoefficients(y,p)
% Recalibrate logit(p) with an intercept and slope; ideal values are 0 and 1.
p=min(max(p,1e-6),1-1e-6);
x=log(p./(1-p));
objective=@(b) bernoulliNLL(y,b(1)+b(2).*x);
options=optimset('Display','off','TolX',1e-7,'TolFun',1e-7,'MaxIter',1000);
b=fminsearch(objective,[0 1],options);
intercept=b(1); slope=b(2);
end

function nll=bernoulliNLL(y,eta)
q=1./(1+exp(-eta));
q=min(max(q,1e-10),1-1e-10);
nll=-sum(y.*log(q)+(1-y).*log(1-q));
end
