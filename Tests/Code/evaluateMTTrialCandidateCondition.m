function result = evaluateMTTrialCandidateCondition(condition,mtModel, ...
  sigmaReadoutDeg,rectifyCandidates)
% evaluateMTTrialCandidateCondition  Deterministic candidate activations.
%
% Candidate directions comprise every physical component direction and its
% opponent. Coincident directions are collapsed before readout; roleLabels
% retain all roles represented by a collapsed mechanism.

arguments
  condition (1,1) struct
  mtModel (1,1) struct
  sigmaReadoutDeg (1,1) double {mustBeFinite,mustBePositive}
  rectifyCandidates (1,1) logical = false
end

[candidateDirections,roleLabels,isPhysical,isOpponentOnly] = ...
  candidateDefinitions(condition.componentDirectionsDeg,condition.componentRoles);
bank = makeGaussianReadoutBank(candidateDirections,sigmaReadoutDeg,mtModel);

change = evaluatePatch(condition.changeComponentCoherence, ...
  condition.changeSignalCoherence,condition.componentDirectionsDeg, ...
  mtModel,bank,rectifyCandidates);
noChange = evaluatePatch(condition.noChangeComponentCoherence, ...
  condition.noChangeSignalCoherence,condition.componentDirectionsDeg, ...
  mtModel,bank,rectifyCandidates);

result = struct();
result.dataset = condition.dataset;
result.offsetDeg = condition.offsetDeg;
result.sigmaMTDeg = mtModel.sigmaMTDeg;
result.sigmaReadoutDeg = sigmaReadoutDeg;
result.rectifyCandidates = rectifyCandidates;
result.sessionIndex = condition.sessionIndex;
result.trialIndex = condition.trialIndex;
result.nExcludedNonfinite = condition.nExcludedNonfinite;
result.pedestalPC = condition.pedestalPC;
result.candidateDirectionsDeg = candidateDirections;
result.roleLabels = roleLabels;
result.isPhysicalCandidate = isPhysical;
result.isOpponentOnlyCandidate = isOpponentOnly;
result.change = change;
result.noChange = noChange;
result.patchMargin = change.maximum-noChange.maximum;
result.changePatchSelected = result.patchMargin>0;
end

function patch = evaluatePatch(coherence,signal,directions,mtModel,bank,rectify)
M = evaluateMTCoherenceDistribution(directions,coherence,mtModel);
M0 = evaluateMTCoherenceDistribution(directions,signal,mtModel);
A = evaluateGaussianDirectionCandidates(M,bank, ...
  'RectifyCandidates',rectify,'PoolingMode','max');
A0 = evaluateGaussianDirectionCandidates(M0,bank, ...
  'RectifyCandidates',rectify,'PoolingMode','max');
patch = struct('evidence',A.evidence,'activation',A.activation, ...
  'maximum',A.pooled,'winner',A.winner, ...
  'signalEvidence',A0.evidence,'signalActivation',A0.activation, ...
  'signalMaximum',A0.pooled,'signalWinner',A0.winner, ...
  'signalWinnerUnique',uniqueMaximum(A0.activation), ...
  'winnerChangedByNoise',double(A.winner~=A0.winner));
patch.winnerChangedByNoise(~patch.signalWinnerUnique)=nan;
end

function tf=uniqueMaximum(x)
if size(x,2)<2,tf=true(size(x,1),1);return,end
s=sort(x,2,'descend');
scale=max(1,max(abs(s(:,1:2)),[],2));
tf=abs(s(:,1)-s(:,2))>1e-12.*scale;
end

function [directions,labels,isPhysical,isOpponentOnly] = ...
    candidateDefinitions(physicalDirections,physicalRoles)
physicalDirections = wrap(double(physicalDirections(:)));
physicalRoles = string(physicalRoles(:));
rawDirections = [physicalDirections;wrap(physicalDirections+180)];
rawRoles = [physicalRoles;"opponentOf_"+physicalRoles];

directions = zeros(0,1); labels = strings(0,1);
isPhysical = false(0,1); hasOpponent = false(0,1);
for k=1:numel(rawDirections)
  row = find(abs(wrap(directions-rawDirections(k)))<1e-9,1,'first');
  if isempty(row)
    directions(end+1,1)=rawDirections(k); %#ok<AGROW>
    labels(end+1,1)=rawRoles(k); %#ok<AGROW>
    isPhysical(end+1,1)=k<=numel(physicalDirections); %#ok<AGROW>
    hasOpponent(end+1,1)=k>numel(physicalDirections); %#ok<AGROW>
  else
    labels(row)=labels(row)+"|"+rawRoles(k);
    isPhysical(row)=isPhysical(row)||(k<=numel(physicalDirections));
    hasOpponent(row)=hasOpponent(row)||(k>numel(physicalDirections));
  end
end
isOpponentOnly = hasOpponent & ~isPhysical;
end

function x=wrap(x)
x=mod(x+180,360)-180;
% Use +180 as the canonical representation of the null direction.
x(abs(x+180)<1e-12)=180;
end
