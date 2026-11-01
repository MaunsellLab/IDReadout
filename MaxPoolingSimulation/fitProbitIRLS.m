function fit = fitProbitIRLS(predictors, outcome)
% FITPROBITIRLS  Toolbox-independent probit GLM using Fisher scoring/IRLS.
%
% Returns stream coefficients separately from the intercept.  A very small
% ridge stabilizes rare nearly separated simulations but does not penalize the
% intercept.

outcome = double(outcome(:));
X = [ones(size(predictors,1),1), double(predictors)];
p = size(X,2);
mu0 = min(max(mean(outcome), 1e-6), 1-1e-6);
beta = zeros(p,1);
beta(1) = -sqrt(2)*erfcinv(2*mu0);
ridge = 1e-9 * diag([0; ones(p-1,1)]);
converged = false;

for iter = 1:100
    eta = X*beta;
    mu = 0.5*erfc(-eta/sqrt(2));
    mu = min(max(mu,1e-9),1-1e-9);
    density = exp(-0.5*eta.^2)/sqrt(2*pi);
    density = max(density,1e-12);
    weights = density.^2 ./ (mu.*(1-mu));
    adjusted = eta + (outcome-mu)./density;
    information = X'*(weights.*X) + ridge;
    betaNew = information \ (X'*(weights.*adjusted));
    if norm(betaNew-beta) <= 1e-8*(1+norm(beta))
        beta = betaNew;
        converged = true;
        break;
    end
    beta = betaNew;
end

eta = X*beta;
mu = min(max(0.5*erfc(-eta/sqrt(2)),1e-12),1-1e-12);
density = max(exp(-0.5*eta.^2)/sqrt(2*pi),1e-12);
weights = density.^2 ./ (mu.*(1-mu));
information = X'*(weights.*X) + ridge;
covarianceAll = pinv(information);

fit = struct();
fit.intercept = beta(1);
fit.beta = beta(2:end);
fit.se = sqrt(max(diag(covarianceAll(2:end,2:end)),0));
fit.covariance = covarianceAll(2:end,2:end);
fit.logLikelihood = sum(outcome.*log(mu) + (1-outcome).*log(1-mu));
fit.converged = converged;
fit.iterations = iter;
end

