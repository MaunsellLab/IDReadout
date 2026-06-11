function results = simulateThreeDirectionReadoutMatched()
% simulateThreeDirectionReadoutMatched
%
% Matched-performance simulation for the 3-direction unpredictable
% coherence-increment task.
%
% For each model, r120, gamma, and target performance, the program finds
% the coherence step that gives the requested mean performance. It then
% simulates analyzed trials at that coherence and estimates:
%
%   correct ~ C0 + C120sum
%
% where C120sum = Cp120 + Cm120. Because C120sum is used rather than the
% mean, a flat three-direction readout predicts beta120 = beta0.
%
% Models:
%   randomReadout  One ordinary readout is selected randomly per trial.
%   flatReadout    All three candidate directions are weighted equally.
%   ordinaryScale  The ordinary directional scale is retained.
%
% Chance performance is 50%. Signal and noise effects are separated:
%   eta = signal term + gamma * noise term
%
% The provisional signal calibration is linear in log odds:
%   S(c) = logit(0.75) * c/20
% so a 20% step gives 75% correct in the single-direction condition.
%
% Top-up trials are not simulated. completedTrialsEquivalent reports how
% many total trials would be needed if 15% are easy top-up trials excluded
% from analysis.

rng(1);

r120List      = [0.10 0.20 0.25];
targetPerfList = [0.70 0.75];
gammaList     = [0.10 0.20 0.30];
nTrialsList   = [10000 40000];
modelList     = {'randomReadout', 'flatReadout', 'ordinaryScale'};

nRep = 100;
noiseSD = 1;
topUpFraction = 0.15;

singleDirectionCoherence75 = 20;
signalAt75 = log(0.75 / 0.25);

% Stable Monte Carlo samples for matching performance.
nTune = 500000;
tuneNoise.C0    = noiseSD * randn(nTune,1);
tuneNoise.Cp120 = noiseSD * randn(nTune,1);
tuneNoise.Cm120 = noiseSD * randn(nTune,1);
tuneNoise.R     = randi(3,nTune,1);

row = 0;
results = table();

for iR = 1:numel(r120List)
    r120 = r120List(iR);

    for iP = 1:numel(targetPerfList)
        targetPerf = targetPerfList(iP);

        for iG = 1:numel(gammaList)
            gamma = gammaList(iG);

            for iM = 1:numel(modelList)
                modelName = modelList{iM};

                [coherence,S,predPerf,pMatch,pWrong] = findMatchedCoherence( ...
                    modelName,r120,targetPerf,gamma,tuneNoise, ...
                    singleDirectionCoherence75,signalAt75);

                for iN = 1:numel(nTrialsList)
                    nTrials = nTrialsList(iN);

                    perfHat = nan(nRep,1);
                    beta0Hat = nan(nRep,1);
                    beta120Hat = nan(nRep,1);
                    contrastFlatHat = nan(nRep,1);
                    contrastOrdinaryHat = nan(nRep,1);
                    pFlatHat = nan(nRep,1);
                    pOrdinaryHat = nan(nRep,1);
                    pMatchHat = nan(nRep,1);
                    pWrongHat = nan(nRep,1);

                    for iRep = 1:nRep
                        sim = simulateOneDataset( ...
                            modelName,nTrials,r120,S,gamma,noiseSD);

                        C120sum = sim.Cp120 + sim.Cm120;
                        X = [sim.C0,C120sum];
                        [b,~,stats] = glmfit( ...
                            X,sim.correct,'binomial','link','logit');

                        beta0 = b(2);
                        beta120 = b(3);
                        covSlopes = stats.covb(2:3,2:3);

                        % Flat-readout contrast: beta0 - beta120 = 0.
                        Lflat = [1 -1];
                        contrastFlat = Lflat * [beta0;beta120];
                        seFlat = sqrt(Lflat * covSlopes * Lflat');
                        zFlat = contrastFlat / seFlat;

                        % Ordinary-scale contrast: beta120-r120*beta0 = 0.
                        Lordinary = [-r120 1];
                        contrastOrdinary = Lordinary * [beta0;beta120];
                        seOrdinary = sqrt(Lordinary * covSlopes * Lordinary');
                        zOrdinary = contrastOrdinary / seOrdinary;

                        perfHat(iRep) = mean(sim.correct);
                        beta0Hat(iRep) = beta0;
                        beta120Hat(iRep) = beta120;
                        contrastFlatHat(iRep) = contrastFlat;
                        contrastOrdinaryHat(iRep) = contrastOrdinary;
                        pFlatHat(iRep) = 2*normcdf(-abs(zFlat));
                        pOrdinaryHat(iRep) = 2*normcdf(-abs(zOrdinary));

                        if strcmp(modelName,'randomReadout')
                            pMatchHat(iRep) = mean(sim.pCorrect(sim.readoutState==1));
                            pWrongHat(iRep) = mean(sim.pCorrect(sim.readoutState~=1));
                        end
                    end

                    row = row + 1;
                    results.model{row,1} = modelName;
                    results.r120(row,1) = r120;
                    results.targetPerf(row,1) = targetPerf;
                    results.gamma(row,1) = gamma;
                    results.nTrials(row,1) = nTrials;
                    results.completedTrialsEquivalent(row,1) = ...
                        nTrials/(1-topUpFraction);
                    results.coherence(row,1) = coherence;
                    results.S(row,1) = S;
                    results.predictedPerf(row,1) = predPerf;
                    results.predictedPMatch(row,1) = pMatch;
                    results.predictedPWrong(row,1) = pWrong;
                    results.meanPerf(row,1) = mean(perfHat);
                    results.sdPerf(row,1) = std(perfHat);
                    results.meanBeta0(row,1) = mean(beta0Hat);
                    results.sdBeta0(row,1) = std(beta0Hat);
                    results.meanBeta120(row,1) = mean(beta120Hat);
                    results.sdBeta120(row,1) = std(beta120Hat);
                    results.meanFlatContrast(row,1) = mean(contrastFlatHat);
                    results.powerRejectFlat(row,1) = mean(pFlatHat<0.05);
                    results.meanOrdinaryContrast(row,1) = ...
                        mean(contrastOrdinaryHat);
                    results.powerRejectOrdinary(row,1) = ...
                        mean(pOrdinaryHat<0.05);
                    results.meanPMatch(row,1) = mean(pMatchHat,'omitnan');
                    results.meanPWrong(row,1) = mean(pWrongHat,'omitnan');
                end
            end
        end
    end
end

disp(results);
makeSummaryPlots(results);

end


function [coherence,S,meanPerf,pMatch,pWrong] = findMatchedCoherence( ...
    modelName,r120,targetPerf,gamma,tuneNoise,coh75,signalAt75)

perfAtC = @(c) predictedPerformance( ...
    modelName,r120,c,gamma,tuneNoise,coh75,signalAt75);

lo = 0;
hi = 200;
pLo = perfAtC(lo);
pHi = perfAtC(hi);

if targetPerf < pLo || targetPerf > pHi
    error(['Target performance %.3f is not bracketed for %s, ' ...
        'r120 %.2f, gamma %.2f. Range is %.3f to %.3f.'], ...
        targetPerf,modelName,r120,gamma,pLo,pHi);
end

for i = 1:70
    mid = (lo+hi)/2;
    if perfAtC(mid) < targetPerf
        lo = mid;
    else
        hi = mid;
    end
end

coherence = (lo+hi)/2;
S = signalAt75 * coherence/coh75;
[meanPerf,pMatch,pWrong] = predictedPerformance( ...
    modelName,r120,coherence,gamma,tuneNoise,coh75,signalAt75);

end


function [meanPerf,pMatch,pWrong] = predictedPerformance( ...
    modelName,r120,coherence,gamma,tuneNoise,coh75,signalAt75)

S = signalAt75 * coherence/coh75;
C0 = tuneNoise.C0;
Cp120 = tuneNoise.Cp120;
Cm120 = tuneNoise.Cm120;
R = tuneNoise.R;

switch modelName
    case 'randomReadout'
        eta = zeros(size(C0));
        is0 = R==1;
        isP = R==2;
        isM = R==3;
        eta(is0) = S + gamma*C0(is0);
        eta(isP) = r120*S + gamma*Cp120(isP);
        eta(isM) = r120*S + gamma*Cm120(isM);
        p = logisticLocal(eta);
        pMatch = mean(p(is0));
        pWrong = mean(p(~is0));

    case 'flatReadout'
        p = logisticLocal(S + gamma*(C0+Cp120+Cm120));
        pMatch = NaN;
        pWrong = NaN;

    case 'ordinaryScale'
        p = logisticLocal(S + gamma*(C0+r120*(Cp120+Cm120)));
        pMatch = NaN;
        pWrong = NaN;

    otherwise
        error('Unknown modelName: %s',modelName);
end

meanPerf = mean(p);

end


function sim = simulateOneDataset(modelName,nTrials,r120,S,gamma,noiseSD)

C0 = noiseSD*randn(nTrials,1);
Cp120 = noiseSD*randn(nTrials,1);
Cm120 = noiseSD*randn(nTrials,1);
readoutState = nan(nTrials,1);

switch modelName
    case 'randomReadout'
        readoutState = randi(3,nTrials,1);
        eta = zeros(nTrials,1);
        is0 = readoutState==1;
        isP = readoutState==2;
        isM = readoutState==3;
        eta(is0) = S + gamma*C0(is0);
        eta(isP) = r120*S + gamma*Cp120(isP);
        eta(isM) = r120*S + gamma*Cm120(isM);

    case 'flatReadout'
        eta = S + gamma*(C0+Cp120+Cm120);

    case 'ordinaryScale'
        eta = S + gamma*(C0+r120*(Cp120+Cm120));

    otherwise
        error('Unknown modelName: %s',modelName);
end

pCorrect = logisticLocal(eta);
correct = rand(nTrials,1)<pCorrect;

sim.C0 = C0;
sim.Cp120 = Cp120;
sim.Cm120 = Cm120;
sim.pCorrect = pCorrect;
sim.correct = correct;
sim.readoutState = readoutState;

end


function makeSummaryPlots(results)

modelList = unique(results.model,'stable');
r120List = unique(results.r120);
targetPerfList = unique(results.targetPerf);
gammaList = unique(results.gamma);
nTrialsToPlot = max(results.nTrials);

% Matched coherence required by each model.
for iP = 1:numel(targetPerfList)
    targetPerf = targetPerfList(iP);
    figure('Color','w');
    tiledlayout(1,numel(gammaList),'TileSpacing','compact');
    for iG = 1:numel(gammaList)
        gamma = gammaList(iG);
        nexttile;
        hold on;
        for iM = 1:numel(modelList)
            modelName = modelList{iM};
            idx = strcmp(results.model,modelName) & ...
                results.targetPerf==targetPerf & ...
                results.gamma==gamma & ...
                results.nTrials==nTrialsToPlot;
            sub = sortrows(results(idx,:),'r120');
            plot(sub.r120,sub.coherence,'-o','DisplayName',modelName);
        end
        xlabel('r_{120}');
        ylabel('Matched coherence (%)');
        title(sprintf('\\gamma = %.2f',gamma));
        box off;
    end
    legend('Location','best');
    sgtitle(sprintf('Coherence required for %.0f%% correct',100*targetPerf));
end

% Coefficient plane at matched performance.
for iP = 1:numel(targetPerfList)
    targetPerf = targetPerfList(iP);
    for iG = 1:numel(gammaList)
        gamma = gammaList(iG);
        figure('Color','w');
        tiledlayout(1,numel(r120List),'TileSpacing','compact');
        for iR = 1:numel(r120List)
            r120 = r120List(iR);
            nexttile;
            hold on;
            for iM = 1:numel(modelList)
                modelName = modelList{iM};
                idx = strcmp(results.model,modelName) & ...
                    results.r120==r120 & ...
                    results.targetPerf==targetPerf & ...
                    results.gamma==gamma & ...
                    results.nTrials==nTrialsToPlot;
                sub = results(idx,:);
                errorbar(sub.meanBeta0,sub.meanBeta120, ...
                    sub.sdBeta120,sub.sdBeta120, ...
                    sub.sdBeta0,sub.sdBeta0,'o', ...
                    'DisplayName',modelName);
            end
            xlabel('\\beta_0');
            ylabel('\\beta_{120 sum}');
            title(sprintf('r_{120}=%.2f',r120));
            axis square;
            box off;
        end
        legend('Location','best');
        sgtitle(sprintf(['Matched %.0f%% correct, \\gamma=%.2f, ' ...
            'N=%d'],100*targetPerf,gamma,nTrialsToPlot));
    end
end

% Power to reject flat and ordinary-scale coefficient constraints.
for iP = 1:numel(targetPerfList)
    targetPerf = targetPerfList(iP);
    for iG = 1:numel(gammaList)
        gamma = gammaList(iG);
        figure('Color','w');
        tiledlayout(1,2,'TileSpacing','compact');

        nexttile;
        hold on;
        for iM = 1:numel(modelList)
            modelName = modelList{iM};
            idx = strcmp(results.model,modelName) & ...
                results.targetPerf==targetPerf & ...
                results.gamma==gamma & ...
                results.nTrials==nTrialsToPlot;
            sub = sortrows(results(idx,:),'r120');
            plot(sub.r120,sub.powerRejectFlat,'-o', ...
                'DisplayName',modelName);
        end
        xlabel('r_{120}');
        ylabel('P(reject flat constraint)');
        ylim([0 1]);
        box off;

        nexttile;
        hold on;
        for iM = 1:numel(modelList)
            modelName = modelList{iM};
            idx = strcmp(results.model,modelName) & ...
                results.targetPerf==targetPerf & ...
                results.gamma==gamma & ...
                results.nTrials==nTrialsToPlot;
            sub = sortrows(results(idx,:),'r120');
            plot(sub.r120,sub.powerRejectOrdinary,'-o', ...
                'DisplayName',modelName);
        end
        xlabel('r_{120}');
        ylabel('P(reject ordinary constraint)');
        ylim([0 1]);
        box off;

        legend('Location','best');
        sgtitle(sprintf(['Constraint tests, %.0f%% correct, ' ...
            '\\gamma=%.2f, N=%d'], ...
            100*targetPerf,gamma,nTrialsToPlot));
    end
end

% Latent-state performance for random readout.
figure('Color','w');
tiledlayout(numel(targetPerfList),numel(gammaList),'TileSpacing','compact');
for iP = 1:numel(targetPerfList)
    targetPerf = targetPerfList(iP);
    for iG = 1:numel(gammaList)
        gamma = gammaList(iG);
        nexttile;
        idx = strcmp(results.model,'randomReadout') & ...
            results.targetPerf==targetPerf & ...
            results.gamma==gamma & ...
            results.nTrials==nTrialsToPlot;
        sub = sortrows(results(idx,:),'r120');
        hold on;
        plot(sub.r120,sub.meanPMatch,'-o','DisplayName','matched');
        plot(sub.r120,sub.meanPWrong,'-o','DisplayName','wrong');
        yline(targetPerf,'--','overall','HandleVisibility','off');
        xlabel('r_{120}');
        ylabel('P(correct)');
        ylim([0.45 1]);
        title(sprintf('target %.2f, \\gamma %.2f',targetPerf,gamma));
        box off;
    end
end
legend('Location','best');
sgtitle('Random-readout latent-state performance');

end


function y = logisticLocal(x)
y = 1./(1+exp(-x));
end
