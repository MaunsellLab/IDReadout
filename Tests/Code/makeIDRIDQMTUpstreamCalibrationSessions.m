function sessions=makeIDRIDQMTUpstreamCalibrationSessions( ...
  idqSummaryPath,idrSummaryPath,idqPreStepPedestalPC)
% makeIDRIDQMTUpstreamCalibrationSessions  Threshold-epoch specifications.

Q=load(idqSummaryPath,'acrossSummary');
R=load(idrSummaryPath,'inventory');
if ~isfield(Q,'acrossSummary')||~isfield(R,'inventory')
  error('makeIDRIDQMTUpstreamCalibrationSessions:MissingSummary', ...
    'Expected acrossSummary and inventory in the source files.');
end
if ~(isscalar(idqPreStepPedestalPC)&&isfinite(idqPreStepPedestalPC)&& ...
    idqPreStepPedestalPC>=0)
  error('makeIDRIDQMTUpstreamCalibrationSessions:BadIDQPedestal', ...
    'idqPreStepPedestalPC must be a finite nonnegative scalar.');
end

qFits=Q.acrossSummary.psychFit.sessionFits;
requireTableVariables(qFits,{'threshold'},'IDQ sessionFits');
nQ=height(qFits);qSession=(1:nQ)';
if ismember('sessionIndex',qFits.Properties.VariableNames)
  qSession=double(qFits.sessionIndex);
end
Qtable=table(repmat("IDQ",nQ,1),qSession,double(qFits.threshold), ...
  zeros(nQ,1),zeros(nQ,1),repmat(idqPreStepPedestalPC,nQ,1), ...
  'VariableNames',{'dataset','sessionIndex','thresholdPC', ...
  'decisionBaseCohPC','preStepPrimaryCohPC','preStepThreeDirPedestalPC'});

if ~isfield(R.inventory,'sessionPsychFits')||~isfield(R.inventory,'trialTable')
  error('makeIDRIDQMTUpstreamCalibrationSessions:MissingIDRFields', ...
    'IDR inventory lacks sessionPsychFits or trialTable.');
end
rFits=R.inventory.sessionPsychFits;T=R.inventory.trialTable;
requireTableVariables(rFits,{'sessionIndex','threshold'},'IDR sessionPsychFits');
requireTableVariables(T,{'sessionIndex','baseCohPC'},'IDR trialTable');
nR=height(rFits);base=nan(nR,1);
for k=1:nR
  use=double(T.sessionIndex)==double(rFits.sessionIndex(k)) & ...
    isfinite(double(T.baseCohPC));
  if ~any(use)
    error('makeIDRIDQMTUpstreamCalibrationSessions:MissingIDRBaseline', ...
      'No finite baseCohPC for IDR session %g.',rFits.sessionIndex(k));
  end
  base(k)=median(double(T.baseCohPC(use)),'omitnan');
end
Rtable=table(repmat("IDR",nR,1),double(rFits.sessionIndex), ...
  double(rFits.threshold),base,base,zeros(nR,1), ...
  'VariableNames',{'dataset','sessionIndex','thresholdPC', ...
  'decisionBaseCohPC','preStepPrimaryCohPC','preStepThreeDirPedestalPC'});
sessions=[Qtable;Rtable];
valid=isfinite(sessions.thresholdPC)&sessions.thresholdPC>0& ...
  isfinite(sessions.decisionBaseCohPC)&sessions.decisionBaseCohPC>=0;
if ~all(valid)
  error('makeIDRIDQMTUpstreamCalibrationSessions:InvalidSessionValues', ...
    '%d sessions have invalid thresholds or decision baselines.',sum(~valid));
end
sessions.baseToThresholdRatio= ...
  sessions.decisionBaseCohPC./sessions.thresholdPC;
preStepText=sprintf([ ...
  'IDQ: three %.6g%% pedestals on every trial; ' ...
  'IDR: preferred base coherence'],idqPreStepPedestalPC);
sessions.preStepDefinition=repmat(string(preStepText),height(sessions),1);
sessions.decisionEpochDefinition=repmat( ...
  "No-noise threshold calibration after IDQ pedestals terminate at step onset", ...
  height(sessions),1);
end

function requireTableVariables(T,names,label)
if ~istable(T),error('makeIDRIDQMTUpstreamCalibrationSessions:ExpectedTable', ...
  '%s must be a table.',label);end
missing=setdiff(names,T.Properties.VariableNames);
if ~isempty(missing)
  error('makeIDRIDQMTUpstreamCalibrationSessions:MissingVariables', ...
    '%s is missing: %s',label,strjoin(missing,', '));
end
end
