function cache = syntheticCache(varargin)
%SYNTHETICCACHE Small two-animal, three-domain cache for tests and examples.
p=inputParser;
p.addParameter('TrialsPerCondition',60,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('NoNoiseFraction',0, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
p.parse(varargin{:});
animals=["Meetz" "Neesha"]; domains=["IDR" "IDQ" "IDS"];
offsets={ [25 45 90 179],120,[10 25 45 120 170] };
rows={}; dirs={}; cohs={}; preCohs={}; drift=[]; pre=[]; trial=0;
old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>; rng(913,'twister');
for ia=1:2
  for id=1:3
    for off=offsets{id}
      for r=1:p.Results.TrialsPerCondition
        trial=trial+1; change=uint8(1+(rand>.5));
        if domains(id)=="IDR", base=95; else, base=7.071; end
        if domains(id)=="IDR", family="IDR"; else, family="IDQIDS"; end
        step=12+6*rand; pedestal=7.071;
        hasNoise=rand>=p.Results.NoNoiseFraction;
        rows(end+1,:)={animals(ia),domains(id),sprintf('%s_%s_%03d',animals(ia),domains(id),off), ...
          trial,uint8(1),change,false,true,hasNoise,uint8(2),true,animals(ia)+"_"+family, ...
          base,base+step,step,off,0}; %#ok<AGROW>
        for j=1:2
          d=[0 off -off];
          if hasNoise, noise=[2.5*randn 3*randn 3*randn]; else, noise=zeros(1,3); end
          total=[base pedestal pedestal]+noise;
          if j==change, total(1)=total(1)+step; end
          dirs{trial,j}=d; cohs{trial,j}=total; %#ok<AGROW>
          preCohs{trial,j}=[base pedestal pedestal]; %#ok<AGROW>
        end
        drift(trial,:)=[0 0]; pre(trial,:)=[base base]; %#ok<AGROW>
      end
    end
  end
end
T=cell2table(rows,'VariableNames',{'animal','domain','sessionID','trialID', ...
  'choiceRDK','changeRDK','correct','valid','hasStepNoise','stepSignIndex', ...
  'isIncrement','kernelFamily','baseCohPC','stepEndpointCohPC', ...
  'stepMagnitudePC','probeOffsetDeg','driftDirectionDeg'});
T.animal=string(T.animal); T.domain=string(T.domain); T.sessionID=string(T.sessionID);
T.kernelFamily=string(T.kernelFamily); T.correct=logical(T.correct); T.valid=logical(T.valid);
T.hasStepNoise=logical(T.hasStepNoise);
T.isIncrement=logical(T.isIncrement);
cache=struct('schemaVersion',"MaxPoolingTrialCache-v2",'trials',T, ...
  'streamDirDeg',{dirs},'stepCohPC',{cohs},'preCohPC',{preCohs},'driftDirDeg',drift, ...
  'preDriftCohPC',pre,'createdBy',mfilename,'createdAt',datetime('now'));
maxpoolfit.validateCache(cache);
end
