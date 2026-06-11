function results = simulateRotationVsFlat2AFC()
% simulateRotationVsFlat2AFC
%
% Sparse simulation comparing:
%   1) random rotation among three ordinary direction-selective readouts
%   2) flat simultaneous readout across the same three directions
%
% This version uses the correct 2AFC psychometric form with a 50% floor:
%
%   P(correct) = 0.5 + 0.5*logistic(eta)
%
% Therefore:
%   - eta -> -Inf gives chance performance (0.5)
%   - eta = 0 gives 0.75 correct
%   - eta -> +Inf gives 1.0
%
% The main question is whether random rotation pushes:
%   - matched-readout trials toward the upper asymptote
%   - mismatched-readout trials toward the lower asymptote
% thereby reducing both marginal betas, while a flat readout keeps trials
% nearer the high-sensitivity midrange around 75%.
%
% Both models are tuned to approximately the same overall performance.
%
% Trial variables:
%   C0      = change-side noise aligned to the actual drift direction
%   Cp120   = change-side noise at +120 deg
%   Cm120   = change-side noise at -120 deg
%
% Regression:
%   correct ~ C0 + (Cp120 + Cm120)
%
% The beta ratio is not used.

rng(1);

% Sparse first-pass grid
r120List       = [0.10 0.20 0.25];
targetPerfList = 0.75;
gammaList      = [0.25 0.40 0.55];
nTrialsList    = [10000 40000];
modelList      = {'randomRotation','flatReadout'};

nRep = 100;
noiseSD = 1;
topUpFraction = 0.15;

% Calibration assumption:
% In the single-direction task, 20% coherence gives 75% correct.
% With the 2AFC-floor psychometric, 75% corresponds to eta = 0.
% Therefore we need an offset and a coherence gain. We parameterize:
%
%   eta = k*(coherence - threshold75)
%
% so that threshold75 = 20% coherence.
%
% k controls psychometric steepness. This is uncertain, so sweep sparsely.
kList = [0.10 0.14 0.18];   % logit-evidence units per % coherence
threshold75 = 20;

% Calibration sample for performance matching
nCal = 300000;
cal.C0 = noiseSD*randn(nCal,1);
cal.Cp = noiseSD*randn(nCal,1);
cal.Cm = noiseSD*randn(nCal,1);
cal.R  = randi(3,nCal,1);

row = 0;
results = table();

for iR = 1:numel(r120List)
    r120 = r120List(iR);

    for iP = 1:numel(targetPerfList)
        targetPerf = targetPerfList(iP);

        for iK = 1:numel(kList)
            k = kList(iK);
            fprintf('r120 %d of %d; perf %d of %d; k %d of %d\n', ...
              iR, numel(r120List), iP, numel(targetPerfList), iK, numel(kList));

            for iG = 1:numel(gammaList)
                gamma = gammaList(iG);

                % Find coherence separately for each model so that both
                % operate at approximately the same overall performance.
                modelParams = struct();

                for iM = 1:numel(modelList)
                    modelName = modelList{iM};

                    coherence = findCoherenceForTarget( ...
                        modelName,targetPerf,r120,k,threshold75, ...
                        gamma,cal);

                    modelParams.(modelName).coherence = coherence;
                end

                for iM = 1:numel(modelList)
                    modelName = modelList{iM};
                    coherence = modelParams.(modelName).coherence;

                    [predPerf,pMatch,pWrong] = predictedPerformance( ...
                        modelName,coherence,r120,k,threshold75, ...
                        gamma,cal);

                    for iN = 1:numel(nTrialsList)
                        nTrials = nTrialsList(iN);

                        perfHat = nan(nRep,1);
                        beta0Hat = nan(nRep,1);
                        beta120Hat = nan(nRep,1);
                        pMatchHat = nan(nRep,1);
                        pWrongHat = nan(nRep,1);

                        for iRep = 1:nRep
                            sim = simulateDataset( ...
                                modelName,nTrials,coherence,r120,k, ...
                                threshold75,gamma,noiseSD);

                            b = glmfit([sim.C0,sim.Cp+sim.Cm], ...
                                sim.correct,'binomial','link','logit');

                            perfHat(iRep) = mean(sim.correct);
                            beta0Hat(iRep) = b(2);
                            beta120Hat(iRep) = b(3);

                            if strcmp(modelName,'randomRotation')
                                pMatchHat(iRep) = ...
                                    mean(sim.pCorrect(sim.R==1));
                                pWrongHat(iRep) = ...
                                    mean(sim.pCorrect(sim.R~=1));
                            end
                        end

                        row = row+1;

                        results.model{row,1} = modelName;
                        results.r120(row,1) = r120;
                        results.targetPerf(row,1) = targetPerf;
                        results.k(row,1) = k;
                        results.gamma(row,1) = gamma;
                        results.nTrials(row,1) = nTrials;
                        results.completedTrialsEquivalent(row,1) = ...
                            nTrials/(1-topUpFraction);
                        results.coherence(row,1) = coherence;
                        results.predictedPerf(row,1) = predPerf;
                        results.predictedPMatch(row,1) = pMatch;
                        results.predictedPWrong(row,1) = pWrong;

                        results.meanPerf(row,1) = mean(perfHat);
                        results.sdPerf(row,1) = std(perfHat);

                        results.meanBeta0(row,1) = mean(beta0Hat);
                        results.sdBeta0(row,1) = std(beta0Hat);

                        results.meanBeta120(row,1) = mean(beta120Hat);
                        results.sdBeta120(row,1) = std(beta120Hat);

                        results.meanPMatch(row,1) = ...
                            mean(pMatchHat,'omitnan');
                        results.meanPWrong(row,1) = ...
                            mean(pWrongHat,'omitnan');
                    end
                end
            end
        end
    end
end

disp(results);
makeSummaryPlots(results);

end


function coherence = findCoherenceForTarget( ...
    modelName,targetPerf,r120,k,threshold75,gamma,cal)

perfAtC = @(c) predictedPerformanceScalar( ...
    modelName,c,r120,k,threshold75,gamma,cal);

lo = 0;
hi = 200;

% Expand upper bound if needed.
while perfAtC(hi) < targetPerf && hi < 1000
    hi = hi*2;
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

end


function m = predictedPerformanceScalar( ...
    modelName,coherence,r120,k,threshold75,gamma,cal)

[m,~,~] = predictedPerformance( ...
    modelName,coherence,r120,k,threshold75,gamma,cal);

end


function [meanPerf,pMatch,pWrong] = predictedPerformance( ...
    modelName,coherence,r120,k,threshold75,gamma,cal)

switch modelName
    case 'randomRotation'
        eta = randomEta(cal.C0,cal.Cp,cal.Cm,cal.R, ...
            coherence,r120,k,threshold75,gamma);
        p = psych2AFC(eta);
        meanPerf = mean(p);
        pMatch = mean(p(cal.R==1));
        pWrong = mean(p(cal.R~=1));

    case 'flatReadout'
        % Flat readout uses the same mean signal for all directions.
        % Noise streams contribute simultaneously with equal gain.
        etaSignal = k*(coherence-threshold75);
        eta = etaSignal + gamma*(cal.C0+cal.Cp+cal.Cm);
        p = psych2AFC(eta);
        meanPerf = mean(p);
        pMatch = NaN;
        pWrong = NaN;

    otherwise
        error('Unknown modelName: %s',modelName);
end

end


function sim = simulateDataset( ...
    modelName,nTrials,coherence,r120,k,threshold75,gamma,noiseSD)

C0 = noiseSD*randn(nTrials,1);
Cp = noiseSD*randn(nTrials,1);
Cm = noiseSD*randn(nTrials,1);

switch modelName
    case 'randomRotation'
        R = randi(3,nTrials,1);
        eta = randomEta(C0,Cp,Cm,R, ...
            coherence,r120,k,threshold75,gamma);

    case 'flatReadout'
        R = nan(nTrials,1);
        etaSignal = k*(coherence-threshold75);
        eta = etaSignal + gamma*(C0+Cp+Cm);

    otherwise
        error('Unknown modelName: %s',modelName);
end

pCorrect = psych2AFC(eta);
correct = rand(nTrials,1) < pCorrect;

sim.C0 = C0;
sim.Cp = Cp;
sim.Cm = Cm;
sim.R = R;
sim.pCorrect = pCorrect;
sim.correct = correct;

end


function eta = randomEta( ...
    C0,Cp,Cm,R,coherence,r120,k,threshold75,gamma)

% Signal evidence for the matched readout.
etaMatched = k*(coherence-threshold75);

% The mismatched readout receives a scaled coherence signal.
% Scale coherence around zero before applying the common threshold offset.
%
% This is a provisional model:
%   matched:   k*(c - threshold75)
%   wrong:     k*(r120*c - threshold75)
%
% Thus, when r120*c is far below the ordinary 75% threshold, wrong-readout
% trials approach the 50% floor.

etaWrong = k*(r120*coherence-threshold75);

eta = zeros(size(C0));

is0 = R==1;
isP = R==2;
isM = R==3;

eta(is0) = etaMatched + gamma*C0(is0);
eta(isP) = etaWrong + gamma*Cp(isP);
eta(isM) = etaWrong + gamma*Cm(isM);

end


function p = psych2AFC(eta)

% Correct 2AFC psychometric with 50% floor and no lapse term.
p = 0.5 + 0.5./(1+exp(-eta));

end


function makeSummaryPlots(results)

r120List = unique(results.r120);
targetPerfList = unique(results.targetPerf);
kList = unique(results.k);
gammaList = unique(results.gamma);
nTrialsToPlot = max(results.nTrials);

% 1. Required coherence
for iP = 1:numel(targetPerfList)
    targetPerf = targetPerfList(iP);

    figure('Color','w');
    tiledlayout(1,numel(kList),'TileSpacing','compact');

    for iK = 1:numel(kList)
        k = kList(iK);
        nexttile;
        hold on;

        for iG = 1:numel(gammaList)
            gamma = gammaList(iG);

            idx = strcmp(results.model,'randomRotation') & ...
                results.targetPerf==targetPerf & ...
                results.k==k & ...
                results.gamma==gamma & ...
                results.nTrials==nTrialsToPlot;

            sub = sortrows(results(idx,:),'r120');
            plot(sub.r120,sub.coherence,'-o', ...
                'DisplayName',sprintf('\\gamma=%.2f',gamma));
        end

        xlabel('r_{120}');
        ylabel('Coherence for target performance (%)');
        title(sprintf('k=%.2f',k));
        box off;
    end

    legend('Location','best');
    sgtitle(sprintf('Random rotation coherence requirement, target=%.2f', ...
        targetPerf));
end

% 2. Beta plane
for iP = 1:numel(targetPerfList)
    targetPerf = targetPerfList(iP);

    for iK = 1:numel(kList)
        k = kList(iK);

        figure('Color','w');
        tiledlayout(1,numel(gammaList),'TileSpacing','compact');

        for iG = 1:numel(gammaList)
            gamma = gammaList(iG);
            nexttile;
            hold on;

            for iM = 1:2
                if iM==1
                    modelName = 'randomRotation';
                else
                    modelName = 'flatReadout';
                end

                idx = strcmp(results.model,modelName) & ...
                    results.targetPerf==targetPerf & ...
                    results.k==k & ...
                    results.gamma==gamma & ...
                    results.nTrials==nTrialsToPlot;

                sub = sortrows(results(idx,:),'r120');

                plot(sub.meanBeta0,sub.meanBeta120,'-o', ...
                    'DisplayName',modelName);

                for i = 1:height(sub)
                    text(sub.meanBeta0(i),sub.meanBeta120(i), ...
                        sprintf(' %.2f',sub.r120(i)), ...
                        'FontSize',8);
                end
            end

            xlabel('\beta_0');
            ylabel('\beta_{120 sum}');
            title(sprintf('\\gamma=%.2f',gamma));
            axis square;
            box off;
        end

        legend('Location','best');
        sgtitle(sprintf('Beta plane, target=%.2f, k=%.2f, N=%d', ...
            targetPerf,k,nTrialsToPlot));
    end
end

% 3. Matched and wrong latent performance for random rotation
for iP = 1:numel(targetPerfList)
    targetPerf = targetPerfList(iP);

    for iK = 1:numel(kList)
        k = kList(iK);

        figure('Color','w');
        tiledlayout(1,numel(gammaList),'TileSpacing','compact');

        for iG = 1:numel(gammaList)
            gamma = gammaList(iG);
            nexttile;
            hold on;

            idx = strcmp(results.model,'randomRotation') & ...
                results.targetPerf==targetPerf & ...
                results.k==k & ...
                results.gamma==gamma & ...
                results.nTrials==nTrialsToPlot;

            sub = sortrows(results(idx,:),'r120');

            plot(sub.r120,sub.meanPMatch,'-o', ...
                'DisplayName','matched');
            plot(sub.r120,sub.meanPWrong,'-o', ...
                'DisplayName','wrong');
            yline(targetPerf,'--','overall target', ...
                'HandleVisibility','off');

            xlabel('r_{120}');
            ylabel('P(correct)');
            ylim([0.48 1]);
            title(sprintf('\\gamma=%.2f',gamma));
            box off;
        end

        legend('Location','best');
        sgtitle(sprintf('Latent-state performance, target=%.2f, k=%.2f', ...
            targetPerf,k));
    end
end

end
