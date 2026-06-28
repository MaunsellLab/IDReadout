% betaFullAnalysis
% Rebuild and run the preferred-direction beta-analysis workflow.
%
% This script assumes that makeBetaSessionData(true) uses
% selectAnalysisFiles and purges only BetaAnalysis-derived outputs.
%
% Set these before running if desired.
nBootstrap = 1000;
bootstrapSeed = 1;
convergenceTolerance = 1e-4;
maxFitPasses = 6;

fprintf('\n=== 1. Rebuild selected session data ===\n');
makeBetaSessionData(true);

fprintf('\n=== 2. Build all-session kernel ===\n');
kernelData = makeBetaKernel(); %#ok<NASGU>

fprintf('\n=== 4. Compute trialwise effective coherence ===\n');
updateBetaSessionEffectiveCoherence(true);

fprintf('\n=== 5. Fit pooled psychometric beta ===\n');
for pass = 1:maxFitPasses
  psychometricFit = fitBetaPsychometric();
  if psychometricFit.gradientInfNorm <= convergenceTolerance
    break;
  end
end
if psychometricFit.gradientInfNorm > convergenceTolerance
  warning('betaFullAnalysis:PsychometricConvergence', ...
    'Psychometric gradient infinity norm remains %.6g.', ...
    psychometricFit.gradientInfNorm);
end

fprintf('\n=== 6. Fit pooled noise beta ===\n');
for pass = 1:maxFitPasses
  noiseFit = fitBetaNoiseRegression();
  if noiseFit.gradientInfNorm <= convergenceTolerance
    break;
  end
end
if noiseFit.gradientInfNorm > convergenceTolerance
  warning('betaFullAnalysis:NoiseConvergence', ...
    'Noise-fit gradient infinity norm remains %.6g.', ...
    noiseFit.gradientInfNorm);
end

fprintf('\n=== 7. Paired hierarchical bootstrap ===\n');
agreement = bootstrapBetaAgreement(nBootstrap, bootstrapSeed); %#ok<NASGU>

fprintf('\n=== 8. Fit session-by-session noise betas ===\n');
sessionFits = fitBetaNoiseBySession(); %#ok<NASGU>

fprintf('\n=== Beta analysis complete ===\n');
fprintf('Psychometric beta: %.6g\n', psychometricFit.beta);
fprintf('Noise beta:        %.6g\n', noiseFit.betaNoise);
fprintf('Difference:        %.6g\n', agreement.difference);
fprintf('Ratio:             %.6g\n', agreement.ratio);
