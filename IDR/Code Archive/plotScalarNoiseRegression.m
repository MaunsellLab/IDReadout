function fig = plotScalarNoiseRegression(reg)
% plotScalarNoiseRegression  Plot binned performance and logistic fits.

T = reg.trialTable;

fig = figure('Color', 'w', 'Position', [100 100 1150 380]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

plotOne(T.DchangePref, T.correct, reg.bins.changePref, reg.models.changePref, ...
    'DchangePref', 'Change-side preferred noise');

plotOne(T.DnoChangePref, T.correct, reg.bins.noChangePref, reg.models.noChangePref, ...
    'DnoChangePref', 'No-change-side preferred noise');

plotOne(T.DdiffPref, T.correct, reg.bins.diffPref, reg.models.diffPref, ...
    'DdiffPref', 'Symmetric evidence Dchange + DnoChange');

if isfield(reg, 'sessionProbeHeader') && isfield(reg.sessionProbeHeader, 'sessionID')
    sessionID = reg.sessionProbeHeader.sessionID;
else
    sessionID = '';
end

if isfield(reg, 'sessionProbeHeader') && isfield(reg.sessionProbeHeader, 'probeTag')
    probeTag = reg.sessionProbeHeader.probeTag;
else
    probeTag = '';
end

sgtitle(sprintf('%s %s scalar noise regression, n=%d, window=%s', ...
    sessionID, probeTag, reg.nTrials, reg.params.noiseWindowName), ...
    'Interpreter', 'none');

end

% =========================================================================
function plotOne(x, y, bins, glm, predictorName, titleText)

nexttile;
hold on;

errorbar(bins.xMean, bins.pCorrect, bins.pSem, 'o', ...
    'LineStyle', 'none', 'MarkerFaceColor', 'auto');

xGrid = linspace(min(x), max(x), 200)';
TT = table(xGrid, 'VariableNames', {predictorName});
pHat = predict(glm, TT);
plot(xGrid, pHat, '-', 'LineWidth', 1.5);

yline(mean(y), ':');

xlabel(sprintf('%s\nmean %% coherence, favor-correct signed', predictorName), ...
    'Interpreter', 'none');
ylabel('P(correct)');
title(titleText, 'Interpreter', 'none');

ylim([0.4 1.0]);
box off;

coef = glm.Coefficients;
row = strcmp(coef.Properties.RowNames, predictorName);
if any(row)
    beta = coef.Estimate(row);
    se = coef.SE(row);
    p = coef.pValue(row);
    txt = sprintf('beta = %.3g\nSE = %.3g\np = %.3g', beta, se, p);

    xl = xlim;
    yl = ylim;
    text(xl(1) + 0.05 * range(xl), yl(2) - 0.08 * range(yl), txt, ...
        'VerticalAlignment', 'top', 'FontSize', 9);
end

end
