function results = simulateThreeDirectionReadout()
% simulateThreeDirectionReadout
%
% First-pass design simulation for 3-direction unpredictable coherence
% increment task.
%
% Each trial has:
%   C0      = change-side noise aligned to actual drift direction
%   Cp120   = change-side noise at +120 deg
%   Cm120   = change-side noise at -120 deg
%
% Models:
%   1) randomReadout:
%        animal randomly selects one of the three ordinary readouts.
%        If selected readout matches actual drift: signal = S.
%        If selected readout is +/-120: signal = r120*S.
%
%   2) flatReadout:
%        animal uses all three channels equally.
%        evidence = S + gamma*(C0 + Cp120 + Cm120)
%
%   3) ordinaryScale:
%        ordinary aligned scale-function model.
%        evidence = S + gamma*(C0 + r120*Cp120 + r120*Cm120)
%
% Analysis fit:
%   correct ~ C0 + (Cp120 + Cm120)
%
% With this convention:
%   flatReadout predicts beta120sum / beta0 ~= 1
%   ordinaryScale predicts beta120sum / beta0 ~= r120
%
% The randomReadout model can produce different patterns because only one
% channel controls a trial.

rng(1);

% -----------------------------
% User-adjustable parameters
% -----------------------------

r120List       = [0.05 0.10 0.20 0.30];
targetPerfList = [0.60 0.65 0.70 0.75];

% Include the 15% trial-loss comparison explicitly.
nTrialsList = [10000 20000 40000 80000];

modelList = {'randomReadout', 'flatReadout', 'ordinaryScale'};

nRep = 200;              % increase later for final power estimates
noiseSD = 1;             % arbitrary units
gamma = 0.20;            % noise-to-log-odds gain; sweep later if desired

% Saturation thresholds for random-readout diagnostic.
satWarn = 0.95;
satBad  = 0.97;

% Monte Carlo samples used to tune S for target performance.
nTune = 300000;

% -----------------------------
% Run simulation
% -----------------------------

row = 0;
results = table();

for iR = 1:numel(r120List)
    r120 = r120List(iR);

    for iP = 1:numel(targetPerfList)
        targetPerf = targetPerfList(iP);

        for iM = 1:numel(modelList)
            modelName = modelList{iM};

            % Tune signal strength for this model/target/r120.
            [S, tuneInfo] = tuneSignalForTarget( ...
                modelName, r120, targetPerf, gamma, noiseSD, nTune);

            for iN = 1:numel(nTrialsList)
                nTrials = nTrialsList(iN);

                beta0Hat      = nan(nRep,1);
                beta120Hat    = nan(nRep,1);
                ratioHat      = nan(nRep,1);
                pCorrectHat   = nan(nRep,1);
                glmP120       = nan(nRep,1);

                for iRep = 1:nRep
                    sim = simulateOneDataset( ...
                        modelName, nTrials, r120, S, gamma, noiseSD);

                    C120sum = sim.Cp120 + sim.Cm120;

                    % GLM: correct ~ C0 + C120sum
                    X = [sim.C0(:), C120sum(:)];
                    y = sim.correct(:);

                    b = glmfit(X, y, 'binomial', 'link', 'logit');
                    beta0Hat(iRep)   = b(2);
                    beta120Hat(iRep) = b(3);
                    ratioHat(iRep)   = b(3) / b(2);

                    % Get p-value for C120 term from fitglm for convenience.
                    % If this is too slow, remove this block.
                    try
                        tbl = table(y, sim.C0(:), C120sum(:), ...
                            'VariableNames', {'correct','C0','C120sum'});
                        mdl = fitglm(tbl, ...
                            'correct ~ C0 + C120sum', ...
                            'Distribution', 'binomial', ...
                            'Link', 'logit');
                        glmP120(iRep) = mdl.Coefficients.pValue(3);
                    catch
                        glmP120(iRep) = NaN;
                    end

                    pCorrectHat(iRep) = mean(y);
                end

                row = row + 1;

                results.model{row,1}       = modelName;
                results.r120(row,1)        = r120;
                results.targetPerf(row,1)  = targetPerf;
                results.nTrials(row,1)     = nTrials;
                results.S(row,1)           = S;
                results.gamma(row,1)       = gamma;
                results.noiseSD(row,1)     = noiseSD;

                results.tunedPerf(row,1)   = tuneInfo.meanPerf;
                results.pMatch(row,1)      = tuneInfo.pMatch;
                results.pWrong(row,1)      = tuneInfo.pWrong;

                results.meanPerf(row,1)    = mean(pCorrectHat);
                results.sdPerf(row,1)      = std(pCorrectHat);

                results.meanBeta0(row,1)   = mean(beta0Hat, 'omitnan');
                results.sdBeta0(row,1)     = std(beta0Hat, 'omitnan');

                results.meanBeta120(row,1) = mean(beta120Hat, 'omitnan');
                results.sdBeta120(row,1)   = std(beta120Hat, 'omitnan');

                results.meanRatio(row,1)   = mean(ratioHat, 'omitnan');
                results.sdRatio(row,1)     = std(ratioHat, 'omitnan');

                results.powerP120(row,1)   = mean(glmP120 < 0.05, 'omitnan');

                if strcmp(modelName, 'randomReadout')
                    if tuneInfo.pMatch >= satBad
                        results.saturationFlag{row,1} = 'bad';
                    elseif tuneInfo.pMatch >= satWarn
                        results.saturationFlag{row,1} = 'warning';
                    else
                        results.saturationFlag{row,1} = 'ok';
                    end
                else
                    results.saturationFlag{row,1} = 'n/a';
                end
            end
        end
    end
end

disp(results);

makeSummaryPlots(results);

end

function sim = simulateOneDataset(modelName, nTrials, r120, S, gamma, noiseSD)

C0    = noiseSD * randn(nTrials,1);
Cp120 = noiseSD * randn(nTrials,1);
Cm120 = noiseSD * randn(nTrials,1);

switch modelName

  case 'randomReadout'
    % 1 = actual drift direction
    % 2 = +120 readout
    % 3 = -120 readout
    R = randi(3, nTrials, 1);

    signal = zeros(nTrials,1);
    noiseTerm = zeros(nTrials,1);

    is0 = R == 1;
    isP = R == 2;
    isM = R == 3;

    signal(is0) = S;
    signal(isP | isM) = r120 * S;

    noiseTerm(is0) = gamma * C0(is0);
    noiseTerm(isP) = gamma * Cp120(isP);
    noiseTerm(isM) = gamma * Cm120(isM);

    eta = signal + noiseTerm;

  case 'flatReadout'
    eta = S + gamma * (C0 + Cp120 + Cm120);

  case 'ordinaryScale'
    eta = S + gamma * (C0 + r120*Cp120 + r120*Cm120);

  otherwise
    error('Unknown modelName: %s', modelName);
end

pCorrect = logisticLocal(eta);
correct = rand(nTrials,1) < pCorrect;

sim.C0 = C0;
sim.Cp120 = Cp120;
sim.Cm120 = Cm120;
sim.correct = correct;
sim.pCorrect = pCorrect;

end

function [S, info] = tuneSignalForTarget(modelName, r120, targetPerf, gamma, noiseSD, nTune)

% Use fixed noise samples for stable tuning.
C0    = noiseSD * randn(nTune,1);
Cp120 = noiseSD * randn(nTune,1);
Cm120 = noiseSD * randn(nTune,1);

% Fixed latent readout states for randomReadout.
R = randi(3, nTune, 1);

perfAtS = @(S) meanPerfGivenS(modelName, S, r120, gamma, C0, Cp120, Cm120, R);

% Bracket S.
lo = -10;
hi = 30;

plo = perfAtS(lo);
phi = perfAtS(hi);

if targetPerf < plo || targetPerf > phi
  warning('Target %.3f not bracketed for %s, r120=%.3f. Range %.3f-%.3f.', ...
    targetPerf, modelName, r120, plo, phi);
end

% Bisection.
for i = 1:80
  mid = (lo + hi)/2;
  pmid = perfAtS(mid);

  if pmid < targetPerf
    lo = mid;
  else
    hi = mid;
  end
end

S = (lo + hi)/2;

[meanPerf, pMatch, pWrong] = meanPerfGivenS( ...
  modelName, S, r120, gamma, C0, Cp120, Cm120, R);

info.meanPerf = meanPerf;
info.pMatch = pMatch;
info.pWrong = pWrong;

end

function [meanPerf, pMatch, pWrong] = meanPerfGivenS( ...
  modelName, S, r120, gamma, C0, Cp120, Cm120, R)

switch modelName

  case 'randomReadout'
    signal = zeros(size(C0));
    noiseTerm = zeros(size(C0));

    is0 = R == 1;
    isP = R == 2;
    isM = R == 3;

    signal(is0) = S;
    signal(isP | isM) = r120 * S;

    noiseTerm(is0) = gamma * C0(is0);
    noiseTerm(isP) = gamma * Cp120(isP);
    noiseTerm(isM) = gamma * Cm120(isM);

    p = logisticLocal(signal + noiseTerm);
    meanPerf = mean(p);

    % Saturation diagnostics.
    pMatch = mean(logisticLocal(S + gamma*C0));
    pWrong = mean(logisticLocal(r120*S + gamma*C0));

  case 'flatReadout'
    p = logisticLocal(S + gamma*(C0 + Cp120 + Cm120));
    meanPerf = mean(p);
    pMatch = NaN;
    pWrong = NaN;

  case 'ordinaryScale'
    p = logisticLocal(S + gamma*(C0 + r120*Cp120 + r120*Cm120));
    meanPerf = mean(p);
    pMatch = NaN;
    pWrong = NaN;

  otherwise
    error('Unknown modelName: %s', modelName);
end

end

function y = logisticLocal(x)

% Numerically safe enough for this range.
y = 1 ./ (1 + exp(-x));

end

function makeSummaryPlots(results)

% Plot beta120/beta0 as a function of r120 for each target performance
% and model. One figure per nTrials.

nTrialsList = unique(results.nTrials);
targetPerfList = unique(results.targetPerf);
modelList = unique(results.model, 'stable');

for iN = 1:numel(nTrialsList)
  nTrials = nTrialsList(iN);

  figure('Color','w');
  tiledlayout(numel(targetPerfList), 1, 'TileSpacing','compact');

  for iP = 1:numel(targetPerfList)
    targetPerf = targetPerfList(iP);
    nexttile;
    hold on;

    for iM = 1:numel(modelList)
      modelName = modelList{iM};

      idx = results.nTrials == nTrials & ...
        results.targetPerf == targetPerf & ...
        strcmp(results.model, modelName);

      sub = results(idx,:);
      [rSorted, order] = sort(sub.r120);

      errorbar(rSorted, ...
        sub.meanRatio(order), ...
        sub.sdRatio(order), ...
        '-o', ...
        'DisplayName', modelName);
    end

    yline(1, ':', 'flat ratio = 1', 'HandleVisibility','off');
    xlabel('r_{120}');
    ylabel('\beta_{120sum} / \beta_0');
    title(sprintf('nTrials = %d, target = %.2f', nTrials, targetPerf));
    legend('Location','best');
    box off;
  end
end

% Saturation diagnostic for random-readout model.
idx = strcmp(results.model, 'randomReadout');
sub = results(idx,:);

figure('Color','w');
hold on;

for iP = 1:numel(targetPerfList)
  targetPerf = targetPerfList(iP);
  idxP = sub.targetPerf == targetPerf & sub.nTrials == min(sub.nTrials);
  this = sub(idxP,:);
  [rSorted, order] = sort(this.r120);

  plot(rSorted, this.pMatch(order), '-o', ...
    'DisplayName', sprintf('target %.2f', targetPerf));
end

yline(0.95, ':', 'warning');
yline(0.97, '--', 'bad');
xlabel('r_{120}');
ylabel('P(correct | matched readout)');
title('Random-readout saturation diagnostic');
legend('Location','best');
box off;

end
