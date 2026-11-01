function winners = classifyDenseWinners(changedIndex,noChangeIndex,bankDirections,streamDirections,radius)
% CLASSIFYDENSEWINNERS  Fractions near drift, either probe, or neither stream.

winners=classifyOne(changedIndex,noChangeIndex,bankDirections,streamDirections,radius);
end

function out=classifyOne(changedIndex,noChangeIndex,bankDirections,streamDirections,radius)
changedDirections=bankDirections(changedIndex(:));
noChangeDirections=bankDirections(noChangeIndex(:));
[out.changedDrift,out.changedProbe,out.changedOther]=fractions(changedDirections,streamDirections,radius);
[out.noChangeDrift,out.noChangeProbe,out.noChangeOther]=fractions(noChangeDirections,streamDirections,radius);
end

function [drift,probe,other]=fractions(winnerDirections,streamDirections,radius)
driftDistance=circularDistance(winnerDirections,streamDirections(1));
probeDistance=min(circularDistance(winnerDirections,streamDirections(2)), ...
    circularDistance(winnerDirections,streamDirections(3)));
isDrift=driftDistance<=radius & driftDistance<=probeDistance;
isProbe=~isDrift & probeDistance<=radius;
drift=mean(isDrift); probe=mean(isProbe); other=1-drift-probe;
end

function d=circularDistance(a,b)
d=abs(mod(a-b+180,360)-180);
end

