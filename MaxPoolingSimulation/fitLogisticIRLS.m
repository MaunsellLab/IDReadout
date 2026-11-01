function fit = fitLogisticIRLS(predictors,outcome)
% FITLOGISTICIRLS  Toolbox-independent logistic regression.

y=double(outcome(:));
X=[ones(size(predictors,1),1),double(predictors)];
p=size(X,2); beta=zeros(p,1); ridge=1e-9*diag([0;ones(p-1,1)]);
converged=false;
for iter=1:100
    eta=max(min(X*beta,35),-35);
    mu=1./(1+exp(-eta));
    w=max(mu.*(1-mu),1e-10);
    z=eta+(y-mu)./w;
    info=X'*(w.*X)+ridge;
    betaNew=info\(X'*(w.*z));
    if norm(betaNew-beta)<=1e-8*(1+norm(beta))
        beta=betaNew; converged=true; break;
    end
    beta=betaNew;
end
eta=max(min(X*beta,35),-35); mu=1./(1+exp(-eta));
ll=sum(y.*log(max(mu,1e-12))+(1-y).*log(max(1-mu,1e-12)));
w=max(mu.*(1-mu),1e-10); covAll=pinv(X'*(w.*X)+ridge);
fit=struct('intercept',beta(1),'beta',beta(2:end), ...
    'covariance',covAll(2:end,2:end),'logLikelihood',ll, ...
    'probability',mu,'pointLogLikelihood', ...
    y.*log(max(mu,1e-12))+(1-y).*log(max(1-mu,1e-12)), ...
    'converged',converged,'iterations',iter);
end
