function [cache,kernels] = crossfitKernels(cache, varargin)
%CROSSFITKERNELS Leave-one-session-out temporal kernels and stream integrals.
%
% Required optional cache fields:
%   raw.kernelDriverByFrame{t}: step-period zero-centered scalar driver used
%     only to estimate temporal weighting. Its sign must be aligned so that a
%     positive value favors the correct/change RDK.
%   raw.stepStreamCohPC{t,j}: nFrames x nStreams total coherence, including
%     pedestals, with columns matching streamDirDeg{t,j}.
%
% All trials in a kernel family must have the same number of step frames.

p=inputParser;
p.addParameter('MinIntegral',1e-6,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
assert(isfield(cache,'raw'),'Cache has no raw framewise fields.');
T=cache.trials; n=height(T); families=unique(string(T.kernelFamily),'stable');
kernels=struct('family',{},'heldOutSession',{},'weights',{},'nTrainingTrials',{});
for f=1:numel(families)
  famMask=string(T.kernelFamily)==families(f) & T.valid & T.hasStepNoise & T.isIncrement;
  familyIdx=find(famMask);
  sessions=unique(string(T.sessionID(familyIdx)),'stable');
  assert(numel(sessions)>=2,'LOSO kernel family %s has fewer than two sessions.',families(f));
  X=stackVectors(cache.raw.kernelDriverByFrame(familyIdx));
  assert(all(isfinite(X),'all'),'Kernel drivers must be finite.');
  y=logical(T.correct(familyIdx));
  sessionForTrial=string(T.sessionID(familyIdx));
  correctTotal=sum(X(:,y),2); errorTotal=sum(X(:,~y),2);
  nCorrectTotal=sum(y); nErrorTotal=sum(~y);
  for s=1:numel(sessions)
    heldOut=sessionForTrial==sessions(s);
    nCorrect=nCorrectTotal-sum(y&heldOut);
    nError=nErrorTotal-sum(~y&heldOut);
    assert(nCorrect>0&&nError>0, ...
      'LOSO fold %s excluding %s lacks correct or error training trials.',families(f),sessions(s));
    correctSum=correctTotal-sum(X(:,y&heldOut),2);
    errorSum=errorTotal-sum(X(:,~y&heldOut),2);
    raw=correctSum./nCorrect-errorSum./nError;
    integral=sum(raw);
    assert(isfinite(integral)&&abs(integral)>p.Results.MinIntegral, ...
      'Kernel integral is unstable for %s excluding %s.',families(f),sessions(s));
    if integral<0, raw=-raw; integral=-integral; end
    w=raw./integral;
    apply=find(string(T.kernelFamily)==families(f) & ...
      string(T.sessionID)==sessions(s) & T.valid);
    for q=apply(:)'
      for j=1:2
        C=cache.raw.stepStreamCohPC{q,j};
        assert(size(C,1)==numel(w),'Frame count differs within kernel family.');
        assert(size(C,2)==numel(cache.streamDirDeg{q,j}), ...
          'Framewise stream count does not match stream directions.');
        cache.stepCohPC{q,j}=(w(:)'*double(C))';
      end
    end
    kernels(end+1)=struct('family',families(f),'heldOutSession',sessions(s), ...
      'weights',w,'nTrainingTrials',nCorrect+nError); %#ok<AGROW>
  end
end
assert(all(cellfun(@(x)all(isfinite(x)),cache.stepCohPC(:))), ...
  'One or more trials did not receive cross-fitted stream integrals.');
cache.kernelMethod="leave-one-session-out";
cache.kernelCreatedAt=datetime('now');
cache.kernelStatus="current";
end

function X=stackVectors(c)
lengths=cellfun(@numel,c);
assert(numel(unique(lengths))==1,'Kernel drivers have inconsistent frame counts.');
X=nan(lengths(1),numel(c));
for k=1:numel(c), X(:,k)=double(c{k}(:)); end
end
