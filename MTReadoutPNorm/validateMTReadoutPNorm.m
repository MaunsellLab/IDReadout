function validateMTReadoutPNorm()
%VALIDATEMTREADOUTPNORM Analytic checks for pooling and the 2AFC decision.
assert(abs(mtpnorm.pool([3;4],2)-5)<1e-12);
assert(mtpnorm.pool([3;4],Inf)==4);
assert(all(mtpnorm.pool(zeros(3,2),5)==zeros(1,2)),'Zero pool failed.');
% Positive firing-rate increases and pooled changes are distinct operations.
a=[20;10]; b=[10;20];
assert(abs(mtpnorm.pool(b,5)-mtpnorm.pool(a,5))<1e-12);
assert(abs(mtpnorm.pool(max(b-a,0),5)-10)<1e-12);
% Keep signed deltas: the smaller decrease must win, not become a zero tie.
R=zeros(2,2,2,1); R(:,1,1,1)=[20;10]; R(:,1,2,1)=[19;9];
R(:,2,1,1)=[20;10]; R(:,2,2,1)=[15;5]; Z=zeros(2,4,2,2);
assert(mtpnorm.predict(R,5,0,Z)==1,'Signed negative deltas were lost.');
R(:,2,:,:)=R(:,1,:,:);
assert(mtpnorm.predict(R,5,0,Z)==.5,'Equal RDKs must tie.');
% Isotropic fixed-dot input produces its full positive response.
cfg=mtpnorm.defaultConfig(); T=table(true,'VariableNames',{'hasStepNoise'});
C=struct('trials',T,'streamDirDeg',{{[0 120 240]}}, ...
  'preCohPC',{{[0;0;0],[0;0;0]}},'stepCohPC',{{[0;0;0],[0;0;0]}});
[F,A]=mtpnorm.fullResponses(C,cfg);
assert(max(abs(F(:)-100*A.randomDirectionMean))<1e-10);
% Densifying a circular grid does not create independent noise dimensions.
Z=mtpnorm.noiseFields(cfg,64);
assert(isequal(size(Z),[360 64 2 2]));
fprintf('MTReadoutPNorm analytic checks passed.\n');
end
