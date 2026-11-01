function d = computeWinnerDirectionDiagnostic(winnerChanged,winnerNoChange, ...
    directionsDeg,W,binWidthDeg)
% COMPUTEWINNERDIRECTIONDIAGNOSTIC  Winner mass and signed local effects.
%
% The changed- and no-change-patch step winners are binned by preferred
% direction.  W contains the drift, +probe, and -probe response profiles.
% The probe diagnostic is the per-stream average 0.5*(W_plus+W_minus), which
% is also half of the identifiable pair derivative in the yoked condition.

directions=mod(double(directionsDeg(:))+180,360)-180;
nMechanisms=numel(directions); nTrials=numel(winnerChanged);
nBins=round(360/binWidthDeg);
edges=linspace(-180,180,nBins+1);
centers=0.5*(edges(1:end-1)+edges(2:end));
bin=discretize(directions,edges);

pChanged=accumarray(winnerChanged(:),1,[nMechanisms 1])/nTrials;
pNoChange=accumarray(winnerNoChange(:),1,[nMechanisms 1])/nTrials;
pAverage=0.5*(pChanged+pNoChange);
driftProfile=double(W(1,:))';
probeProfile=0.5*double(W(2,:)+W(3,:))';

d=struct();
d.directionDeg=centers(:);
d.changedWinnerFraction=accumarray(bin,pChanged,[nBins 1],@sum,0);
d.noChangeWinnerFraction=accumarray(bin,pNoChange,[nBins 1],@sum,0);
d.patchAverageWinnerFraction=accumarray(bin,pAverage,[nBins 1],@sum,0);
d.driftLocalSensitivity=accumarray(bin,driftProfile,[nBins 1],@mean,NaN);
d.probeLocalSensitivityPerStream=accumarray(bin,probeProfile,[nBins 1],@mean,NaN);
d.driftWinnerWeightedContribution=accumarray(bin,pAverage.*driftProfile, ...
    [nBins 1],@sum,0);
d.probeWinnerWeightedContributionPerStream=accumarray(bin, ...
    pAverage.*probeProfile,[nBins 1],@sum,0);
d.driftFixedWinnerTotal=sum(d.driftWinnerWeightedContribution);
d.probeFixedWinnerTotalPerStream=sum( ...
    d.probeWinnerWeightedContributionPerStream);
end
