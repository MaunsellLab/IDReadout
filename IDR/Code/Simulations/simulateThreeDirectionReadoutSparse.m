function results = simulateThreeDirectionReadoutSparse()
% simulateThreeDirectionReadoutSparse
%
% Sparse first-pass simulation for a 3-direction unpredictable coherence
% increment task.
%
% Assumptions:
%   - Three possible drift directions separated by 120 deg.
%   - Trials are aligned to the actual drift direction.
%   - Independent change-side noise streams: C0, Cp120, Cm120.
%   - +120 and -120 are symmetric.
%   - Chance performance is 50%.
%   - Mean signal and coherence-noise effects are separated:
%         eta = signal term + gamma * noise term
%
% Models:
%   1) randomReadout
%        One of the three ordinary readouts is selected randomly each trial.
%        If the selected readout matches the drift:
%             signal = S(c), noise = C0
%        Otherwise:
%             signal = r120*S(c), noise = corresponding +/-120 stream
%
%   2) flatReadout
%        Signal is read out equally well for any of the three possible
%        directions. All three noise streams contribute equally:
%             eta = S(c) + gamma*(C0 + Cp120 + Cm120)
%
%   3) ordinaryScale
%        Signal is available, but noise weighting retains the ordinary
%        scale-function shape:
%             eta = S(c) + gamma*(C0 + r120*(Cp120 + Cm120))
%
% Provisional signal calibration:
%   A 20% coherence step gives 75% correct in the single-direction task:
%       S(c) = logit(0.75) * c/20
%
% Analysis:
%   correct ~ C0 + C120sum
%   where C120sum = Cp120 + Cm120
%
% With this convention:
%   flatReadout predicts beta120 ~= beta0
%   ordinaryScale predicts beta120 ~= r120*beta0
%
% The beta120/beta0 ratio is deliberately not used as a primary output.

rng(1);

r120List       = [0.05 0.10 0.20];
coherenceList  = [20 35 50 65];
gammaList      = [0.10 0.20 0.30];
nTrialsList    = [10000 40000];
modelList      = {'randomReadout', 'flatReadout', 'ordinaryScale'};

nRep = 100;
noiseSD = 1;

singleDirectionCoherence75 = 20;
signalAt75 = log(0.75 / 0.25);

row = 0;
results = table();

for iR = 1:numel(r120List)
    r120 = r120List(iR);

    for iC = 1:numel(coherenceList)
        coherence = coherenceList(iC);
        S = signalAt75 * coherence / singleDirectionCoherence75;

        for iG = 1:numel(gammaList)
            gamma = gammaList(iG);

            for iM = 1:numel(modelList)
                modelName = modelList{iM};

                for iN = 1:numel(nTrialsList)
                    nTrials = nTrialsList(iN);

                    perfHat       = nan(nRep,1);
                    beta0Hat      = nan(nRep,1);
                    beta120Hat    = nan(nRep,1);
                    seBeta0Hat    = nan(nRep,1);
                    seBeta120Hat  = nan(nRep,1);
                    pBeta0Hat     = nan(nRep,1);
                    pBeta120Hat   = nan(nRep,1);
                    pMatchHat     = nan(nRep,1);
                    pWrongHat     = nan(nRep,1);

                    for iRep = 1:nRep
                        sim = simulateOneDataset( ...
                            modelName, nTrials, r120, S, gamma, noiseSD);

                        C120sum = sim.Cp120 + sim.Cm120;
                        X = [sim.C0, C120sum];

                        [b,~,stats] = glmfit( ...
                            X, sim.correct, 'binomial', 'link', 'logit');

                        perfHat(iRep)      = mean(sim.correct);
                        beta0Hat(iRep)     = b(2);
                        beta120Hat(iRep)   = b(3);
                        seBeta0Hat(iRep)   = stats.se(2);
                        seBeta120Hat(iRep) = stats.se(3);
                        pBeta0Hat(iRep)    = stats.p(2);
                        pBeta120Hat(iRep)  = stats.p(3);

                        if strcmp(modelName, 'randomReadout')
                            pMatchHat(iRep) = mean(sim.pCorrect(sim.readoutState == 1));
                            pWrongHat(iRep) = mean(sim.pCorrect(sim.readoutState ~= 1));
                        end
                    end

                    row = row + 1;
                    results.model{row,1}      = modelName;
                    results.r120(row,1)       = r120;
                    results.coherence(row,1)  = coherence;
                    results.gamma(row,1)      = gamma;
                    results.nTrials(row,1)    = nTrials;
                    results.S(row,1)          = S;
                    results.meanPerf(row,1)   = mean(perfHat);
                    results.sdPerf(row,1)     = std(perfHat);
                    results.meanBeta0(row,1)  = mean(beta0Hat);
                    results.sdBeta0(row,1)    = std(beta0Hat);
                    results.meanSEBeta0(row,1)= mean(seBeta0Hat);
                    results.powerBeta0(row,1) = mean(pBeta0Hat < 0.05);
                    results.meanBeta120(row,1)   = mean(beta120Hat);
                    results.sdBeta120(row,1)     = std(beta120Hat);
                    results.meanSEBeta120(row,1) = mean(seBeta120Hat);
                    results.powerBeta120(row,1)  = mean(pBeta120Hat < 0.05);
                    results.meanPMatch(row,1) = mean(pMatchHat, 'omitnan');
                    results.meanPWrong(row,1) = mean(pWrongHat, 'omitnan');
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
readoutState = nan(nTrials,1);

switch modelName
    case 'randomReadout'
        readoutState = randi(3, nTrials, 1);
        eta = zeros(nTrials,1);
        is0 = readoutState == 1;
        isP = readoutState == 2;
        isM = readoutState == 3;
        eta(is0) = S + gamma*C0(is0);
        eta(isP) = r120*S + gamma*Cp120(isP);
        eta(isM) = r120*S + gamma*Cm120(isM);

    case 'flatReadout'
        eta = S + gamma*(C0 + Cp120 + Cm120);

    case 'ordinaryScale'
        eta = S + gamma*(C0 + r120*(Cp120 + Cm120));

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
sim.readoutState = readoutState;

end


function makeSummaryPlots(results)

modelList = unique(results.model, 'stable');
r120List = unique(results.r120);
gammaList = unique(results.gamma);
nTrialsList = unique(results.nTrials);
nTrialsToPlot = max(nTrialsList);

for iG = 1:numel(gammaList)
    gamma = gammaList(iG);

    figure('Color','w');
    tiledlayout(1, numel(r120List), 'TileSpacing','compact');
    for iR = 1:numel(r120List)
        r120 = r120List(iR);
        nexttile;
        hold on;
        for iM = 1:numel(modelList)
            modelName = modelList{iM};
            idx = strcmp(results.model, modelName) & ...
                  results.r120 == r120 & ...
                  results.gamma == gamma & ...
                  results.nTrials == nTrialsToPlot;
            sub = sortrows(results(idx,:), 'coherence');
            errorbar(sub.coherence, sub.meanPerf, sub.sdPerf, '-o', ...
                'DisplayName', modelName);
        end
        yline(0.50, ':', 'chance', 'HandleVisibility','off');
        yline(0.75, '--', '75%', 'HandleVisibility','off');
        xlabel('Coherence step (%)');
        ylabel('P(correct)');
        title(sprintf('r_{120}=%.2f', r120));
        ylim([0.45 1]);
        box off;
    end
    legend('Location','best');
    sgtitle(sprintf('Performance versus coherence, \\gamma = %.2f, N = %d', ...
        gamma, nTrialsToPlot));

    figure('Color','w');
    tiledlayout(1, numel(r120List), 'TileSpacing','compact');
    for iR = 1:numel(r120List)
        r120 = r120List(iR);
        nexttile;
        hold on;
        for iM = 1:numel(modelList)
            modelName = modelList{iM};
            idx = strcmp(results.model, modelName) & ...
                  results.r120 == r120 & ...
                  results.gamma == gamma & ...
                  results.nTrials == nTrialsToPlot;
            sub = sortrows(results(idx,:), 'coherence');
            errorbar(sub.coherence, sub.meanBeta0, sub.sdBeta0, '-o', ...
                'DisplayName', modelName);
        end
        xlabel('Coherence step (%)');
        ylabel('\\beta_0');
        title(sprintf('r_{120}=%.2f', r120));
        box off;
    end
    legend('Location','best');
    sgtitle(sprintf('\\beta_0 versus coherence, \\gamma = %.2f, N = %d', ...
        gamma, nTrialsToPlot));

    figure('Color','w');
    tiledlayout(1, numel(r120List), 'TileSpacing','compact');
    for iR = 1:numel(r120List)
        r120 = r120List(iR);
        nexttile;
        hold on;
        for iM = 1:numel(modelList)
            modelName = modelList{iM};
            idx = strcmp(results.model, modelName) & ...
                  results.r120 == r120 & ...
                  results.gamma == gamma & ...
                  results.nTrials == nTrialsToPlot;
            sub = sortrows(results(idx,:), 'coherence');
            errorbar(sub.coherence, sub.meanBeta120, sub.sdBeta120, '-o', ...
                'DisplayName', modelName);
        end
        xlabel('Coherence step (%)');
        ylabel('\\beta_{120 sum}');
        title(sprintf('r_{120}=%.2f', r120));
        box off;
    end
    legend('Location','best');
    sgtitle(sprintf('\\beta_{120} versus coherence, \\gamma = %.2f, N = %d', ...
        gamma, nTrialsToPlot));

    figure('Color','w');
    tiledlayout(1, numel(r120List), 'TileSpacing','compact');
    for iR = 1:numel(r120List)
        r120 = r120List(iR);
        nexttile;
        hold on;
        for iM = 1:numel(modelList)
            modelName = modelList{iM};
            idx = strcmp(results.model, modelName) & ...
                  results.r120 == r120 & ...
                  results.gamma == gamma & ...
                  results.nTrials == nTrialsToPlot;
            sub = sortrows(results(idx,:), 'coherence');
            plot(sub.meanBeta0, sub.meanBeta120, '-o', ...
                'DisplayName', modelName);
            for i = 1:height(sub)
                text(sub.meanBeta0(i), sub.meanBeta120(i), ...
                    sprintf(' %d', sub.coherence(i)), 'FontSize', 8);
            end
        end
        xlabel('\\beta_0');
        ylabel('\\beta_{120 sum}');
        title(sprintf('r_{120}=%.2f', r120));
        axis square;
        box off;
    end
    legend('Location','best');
    sgtitle(sprintf('Coefficient plane, \\gamma = %.2f, N = %d', ...
        gamma, nTrialsToPlot));

    figure('Color','w');
    tiledlayout(1, numel(r120List), 'TileSpacing','compact');
    for iR = 1:numel(r120List)
        r120 = r120List(iR);
        nexttile;
        hold on;
        idx = strcmp(results.model, 'randomReadout') & ...
              results.r120 == r120 & ...
              results.gamma == gamma & ...
              results.nTrials == nTrialsToPlot;
        sub = sortrows(results(idx,:), 'coherence');
        plot(sub.coherence, sub.meanPMatch, '-o', ...
            'DisplayName', 'matched readout');
        plot(sub.coherence, sub.meanPWrong, '-o', ...
            'DisplayName', 'wrong readout');
        plot(sub.coherence, sub.meanPerf, '-o', ...
            'DisplayName', 'overall');
        yline(0.50, ':', 'chance', 'HandleVisibility','off');
        yline(0.75, '--', '75%', 'HandleVisibility','off');
        xlabel('Coherence step (%)');
        ylabel('P(correct)');
        title(sprintf('r_{120}=%.2f', r120));
        ylim([0.45 1]);
        box off;
    end
    legend('Location','best');
    sgtitle(sprintf('Random-readout latent states, \\gamma = %.2f, N = %d', ...
        gamma, nTrialsToPlot));
end

end


function y = logisticLocal(x)
y = 1 ./ (1 + exp(-x));
end
