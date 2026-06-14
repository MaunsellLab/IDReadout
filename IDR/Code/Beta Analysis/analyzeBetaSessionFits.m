function sessionBetaAnalysis = analyzeBetaSessionFits()
% analyzeBetaSessionFits
% Summarize session beta distribution, heterogeneity, and dependence on
% trial count, error count, performance, and beta uncertainty.

% cleanupObj = initProjectPath(); %#ok<NASGU>

baseFolder = domainFolder(mfilename('fullpath'));
acrossFolder = fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'AcrossSessions');
inputPath = fullfile(acrossFolder, 'BetaNoiseSessionFits.mat');

S = load(inputPath, 'sessionFitTable');
T = S.sessionFitTable;

valid = isfinite(T.beta) & isfinite(T.betaSE) & T.betaSE > 0 & T.exitflag > 0;
T = T(valid,:);

beta = T.beta;
betaSE = T.betaSE;
v = betaSE.^2;
n = height(T);
referenceBeta = unique(T.referenceBeta);

% Ordinary summaries.
ordinary.mean = mean(beta);
ordinary.median = median(beta);
ordinary.sd = std(beta);
ordinary.sem = ordinary.sd / sqrt(n);
ordinary.ci95 = ordinary.mean + tinv([0.025 0.975], n-1)*ordinary.sem;
ordinary.min = min(beta);
ordinary.max = max(beta);

% Fixed-effect meta-analysis.
w = 1 ./ v;
fixedMean = sum(w.*beta) / sum(w);
fixedSE = sqrt(1/sum(w));
fixedCI95 = fixedMean + 1.96*fixedSE*[-1 1];

% Heterogeneity and DerSimonian-Laird random effects.
Q = sum(w .* (beta-fixedMean).^2);
df = n-1;
C = sum(w) - sum(w.^2)/sum(w);
tau2 = max(0, (Q-df)/C);
I2 = max(0, (Q-df)/Q);

wr = 1 ./ (v + tau2);
randomMean = sum(wr.*beta) / sum(wr);
randomSE = sqrt(1/sum(wr));
randomCI95 = randomMean + 1.96*randomSE*[-1 1];

reliability = tau2 ./ (tau2 + v);
shrunkenBeta = reliability.*beta + (1-reliability).*randomMean;

T.reliability = reliability;
T.shrunkenBeta = shrunkenBeta;
T.betaDifferenceFromReference = beta-referenceBeta;
T.betaRatioToReference = beta/referenceBeta;

A.beta_vs_nTrials = corrSummary(beta,T.nTrials);
A.beta_vs_nError = corrSummary(beta,T.nError);
A.beta_vs_fractionCorrect = corrSummary(beta,T.fractionCorrect);
A.beta_vs_betaSE = corrSummary(beta,betaSE);
A.absDeviation_vs_betaSE = corrSummary(abs(beta-referenceBeta),betaSE);

sessionBetaAnalysis = struct();
sessionBetaAnalysis.version = 1;
sessionBetaAnalysis.referenceBeta = referenceBeta;
sessionBetaAnalysis.ordinary = ordinary;
sessionBetaAnalysis.fixedEffect = struct('mean',fixedMean,'se',fixedSE,'ci95',fixedCI95);
sessionBetaAnalysis.randomEffects = struct( ...
  'mean',randomMean,'se',randomSE,'ci95',randomCI95, ...
  'tau2',tau2,'tau',sqrt(tau2),'Q',Q,'df',df,'I2',I2);
sessionBetaAnalysis.association = A;
sessionBetaAnalysis.sessionTable = T;
sessionBetaAnalysis.createdBy = mfilename;
sessionBetaAnalysis.createdDate = datetime('now');

outputPath = fullfile(acrossFolder, 'BetaSessionFitAnalysis.mat');
save(outputPath, 'sessionBetaAnalysis', '-v7.3');

fprintf('Saved %s\n', outputPath);
fprintf('Reference beta: %.6g\n', referenceBeta);
fprintf('Unweighted mean: %.6g (95%% CI %.6g to %.6g)\n', ...
  ordinary.mean, ordinary.ci95(1), ordinary.ci95(2));
fprintf('Median: %.6g\n', ordinary.median);
fprintf('Fixed-effect mean: %.6g (95%% CI %.6g to %.6g)\n', ...
  fixedMean, fixedCI95(1), fixedCI95(2));
fprintf('Random-effects mean: %.6g (95%% CI %.6g to %.6g)\n', ...
  randomMean, randomCI95(1), randomCI95(2));
fprintf('Between-session SD tau: %.6g\n', sqrt(tau2));
fprintf('I^2: %.3f\n', I2);
fprintf('Median session reliability: %.3f\n', median(reliability));

fprintf('\nSpearman correlations with beta:\n');
fprintf('  nTrials:         rho %.3f, p %.4g\n', A.beta_vs_nTrials.rho, A.beta_vs_nTrials.p);
fprintf('  nError:          rho %.3f, p %.4g\n', A.beta_vs_nError.rho, A.beta_vs_nError.p);
fprintf('  fractionCorrect: rho %.3f, p %.4g\n', ...
  A.beta_vs_fractionCorrect.rho, A.beta_vs_fractionCorrect.p);
fprintf('  betaSE:          rho %.3f, p %.4g\n', A.beta_vs_betaSE.rho, A.beta_vs_betaSE.p);

plotDiagnostics(T, referenceBeta, randomMean);
end

function R = corrSummary(x,y)
use = isfinite(x) & isfinite(y);
[rho,p] = corr(x(use),y(use),'Type','Spearman','Rows','complete');
R = struct('rho',rho,'p',p,'n',sum(use));
end

function plotDiagnostics(T, referenceBeta, randomMean)

figure;
histogram(T.beta);
hold on
xline(referenceBeta,'r--','Reference');
xline(randomMean,'k-','Random-effects mean');
xlabel('Session noise beta');
ylabel('Number of sessions');
title('Distribution of session noise betas');
box off

figure;
scatter(T.betaSE,T.beta,'filled');
hold on
yline(referenceBeta,'r--');
xlabel('Session beta SE');
ylabel('Session beta');
title('Session beta versus uncertainty');
box off

figure;
scatter(T.nError,T.beta,'filled');
hold on
yline(referenceBeta,'r--');
xlabel('Number of error trials');
ylabel('Session beta');
title('Session beta versus error count');
box off

figure;
scatter(T.fractionCorrect,T.beta,'filled');
hold on
yline(referenceBeta,'r--');
xlabel('Fraction correct');
ylabel('Session beta');
title('Session beta versus operating performance');
box off

figure;
scatter(T.beta,T.shrunkenBeta,'filled');
hold on
lims = [min([T.beta;T.shrunkenBeta]), max([T.beta;T.shrunkenBeta])];
plot(lims,lims,'k:');
xline(referenceBeta,'r--');
yline(referenceBeta,'r--');
xlabel('Raw session beta');
ylabel('Shrunken session beta');
title('Random-effects shrinkage');
axis square
box off
end
