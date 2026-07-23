function sessions=makeIDRIDQMTCalibrationSessions(idqSummaryPath,idrSummaryPath)
% makeIDRIDQMTCalibrationSessions  Session thresholds and no-noise baselines.

Q=load(idqSummaryPath,'acrossSummary');
R=load(idrSummaryPath,'inventory');
if ~isfield(Q,'acrossSummary')||~isfield(R,'inventory')
  error('makeIDRIDQMTCalibrationSessions:MissingSummary', ...
    'Expected acrossSummary and inventory in the source files.');
end

qFits=Q.acrossSummary.psychFit.sessionFits;
if ~istable(qFits)||~ismember('threshold',qFits.Properties.VariableNames)
  error('makeIDRIDQMTCalibrationSessions:MissingIDQThreshold', ...
    'IDQ psychFit.sessionFits must contain threshold.');
end
nQ=height(qFits);
qSession=(1:nQ)';
if ismember('sessionIndex',qFits.Properties.VariableNames)
  qSession=double(qFits.sessionIndex);
end
Qtable=table(repmat("IDQ",nQ,1),qSession,double(qFits.threshold), ...
  zeros(nQ,1),'VariableNames', ...
  {'dataset','sessionIndex','thresholdPC','baseCohPC'});

if ~isfield(R.inventory,'sessionPsychFits')|| ...
    ~isfield(R.inventory,'trialTable')
  error('makeIDRIDQMTCalibrationSessions:MissingIDRFields', ...
    'IDR inventory lacks sessionPsychFits or trialTable.');
end
rFits=R.inventory.sessionPsychFits;
requiredFit={'sessionIndex','threshold'};
requiredTrial={'sessionIndex','baseCohPC'};
requireTableVariables(rFits,requiredFit,'IDR sessionPsychFits');
T=R.inventory.trialTable;
requireTableVariables(T,requiredTrial,'IDR trialTable');
nR=height(rFits); base=nan(nR,1);
for k=1:nR
  use=double(T.sessionIndex)==double(rFits.sessionIndex(k)) & ...
    isfinite(double(T.baseCohPC));
  if ~any(use)
    error('makeIDRIDQMTCalibrationSessions:MissingIDRBaseline', ...
      'No finite baseCohPC for IDR session %g.',rFits.sessionIndex(k));
  end
  base(k)=median(double(T.baseCohPC(use)),'omitnan');
end
Rtable=table(repmat("IDR",nR,1),double(rFits.sessionIndex), ...
  double(rFits.threshold),base,'VariableNames', ...
  {'dataset','sessionIndex','thresholdPC','baseCohPC'});

sessions=[Qtable;Rtable];
valid=isfinite(sessions.thresholdPC)&sessions.thresholdPC>0& ...
  isfinite(sessions.baseCohPC)&sessions.baseCohPC>=0;
if ~all(valid)
  error('makeIDRIDQMTCalibrationSessions:InvalidSessionValues', ...
    '%d sessions have invalid thresholds or baselines.',sum(~valid));
end
sessions.baseToThresholdRatio=sessions.baseCohPC./sessions.thresholdPC;
end

function requireTableVariables(T,names,label)
if ~istable(T),error('makeIDRIDQMTCalibrationSessions:ExpectedTable', ...
  '%s must be a table.',label);end
missing=setdiff(names,T.Properties.VariableNames);
if ~isempty(missing)
  error('makeIDRIDQMTCalibrationSessions:MissingVariables', ...
    '%s is missing: %s',label,strjoin(missing,', '));
end
end
