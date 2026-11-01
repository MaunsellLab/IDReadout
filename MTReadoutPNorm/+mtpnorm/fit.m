function result = fit(cache,cfg,pGrid)
%FIT Fit a shared noise SD per p, using bounded-memory trial batches.
T=cache.trials; mask=T.hasStepNoise & T.isIncrement;
assert(any(mask),'No noisy increment trials.');
if ~isfield(cfg,'trialBatchSize'), cfg.trialBatchSize=128; end
choice=double(T.choiceRDK(mask)==1);
Z=mtpnorm.noiseFields(cfg,cfg.fitDraws);
Zfinal=mtpnorm.noiseFields(cfg,cfg.finalDraws);
grid=unique([cfg.noiseBounds cfg.noiseGrid]);
grid=grid(grid>=cfg.noiseBounds(1) & grid<=cfg.noiseBounds(2));
profile=cell(numel(pGrid),1); predictions=cell(numel(pGrid),1);
rows=cell(numel(pGrid),8); responseAudit=struct();
for k=1:numel(pGrid)
  p=pGrid(k); fprintf('Pooling p = %g (%d/%d; %d fit trials)\n',p,k,numel(pGrid),nnz(mask));
  f=@(sd)objective(sd,cache,cfg,p,Z,mask,choice);
  values=arrayfun(f,grid); [~,best]=min(values);
  lo=grid(max(1,best-1)); hi=grid(min(numel(grid),best+1));
  candidates=grid(best);
  if hi>lo
    candidates(end+1)=fminbnd(f,lo,hi,optimset('Display','off','TolX',.01));
  end
  v=arrayfun(f,candidates); [~,q]=min(v); sd=candidates(q);
  % Re-evaluate with higher integration precision; retain both scores.
  [prob,audit,responseAudit]=mtpnorm.predictCache(cache,cfg,p,sd,Zfinal);
  prob=min(max(prob,cfg.probabilityFloor),1-cfg.probabilityFloor);
  y=double(T.choiceRDK==1);
  finalLL=sum(y(mask).*log(prob(mask))+(1-y(mask)).*log(1-prob(mask)));
  pCorrect=prob; c=T.changeRDK==2; pCorrect(c)=1-prob(c);
  rows(k,:)={p,sd,-min(v),finalLL,mean(T.correct(mask)), ...
    mean(pCorrect(mask)),audit.negativeRateFraction,abs(sd-cfg.noiseBounds(2))<.05};
  profile{k}=table(grid(:),values(:),'VariableNames',{'noiseSD','negLogLikelihood'});
  predictions{k}=struct('p',p,'noiseSD',sd,'pRDK1',prob, ...
    'pCorrect',pCorrect,'rateAudit',audit);
  if isfield(cfg,'checkpointFile') && strlength(string(cfg.checkpointFile))>0
    checkpoint=struct('config',cfg,'pGrid',pGrid,'completedThrough',k, ...
      'rows',{rows(1:k,:)},'profiles',{profile(1:k)}, ...
      'predictions',{predictions(1:k)},'sourceAudit',cache.sourceAudit); %#ok<NASGU>
    save(cfg.checkpointFile,'checkpoint','-v7.3');
  end
end
fits=cell2table(rows,'VariableNames',{'p','noiseSD','fitLogLikelihood', ...
  'finalLogLikelihood','observedAccuracy','predictedAccuracy', ...
  'negativeRateFraction','noiseAtUpperBound'});
[~,bestIndex]=max(fits.finalLogLikelihood);
result=struct('config',cfg,'trials',T,'fitMask',mask,'fits',fits, ...
  'noiseProfiles',{profile},'predictions',{predictions}, ...
  'bestIndex',bestIndex,'responseAudit',responseAudit, ...
  'fitDraws',cfg.fitDraws,'finalDraws',cfg.finalDraws,'createdAt',datetime('now'));
end

function value=objective(sd,cache,cfg,p,Z,mask,y)
started=tic;
prob=mtpnorm.predictCache(cache,cfg,p,sd,Z,find(mask));
prob=min(max(prob,cfg.probabilityFloor),1-cfg.probabilityFloor);
value=-sum(y.*log(prob)+(1-y).*log(1-prob));
fprintf('  SD=%7.3f  NLL=%12.4f  %.1fs\n',sd,value,toc(started));
end
