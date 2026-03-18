function out = PvIN_Plot(figureNum, dn_signed, isCorrect, nBins, titleStr)

if nargin < 4 || isempty(nBins)
    nBins = 10;
end

isCorrect = isCorrect(:);
S = sum(dn_signed, 1)';

if numel(S) ~= numel(isCorrect)
    error('PvIN_Plot:SizeMismatch', ...
        'dn_signed and isCorrect must contain the same number of trials.');
end

valid = isfinite(S) & isfinite(isCorrect);
S = S(valid);
isCorrect = isCorrect(valid);

N = numel(S);
if N < nBins
    error('PvIN_Plot:TooFewTrials', ...
        'Number of trials (%d) is smaller than number of bins (%d).', N, nBins);
end

[Ss, sortIdx] = sort(S);
correctSorted = isCorrect(sortIdx);

binEdges = round(linspace(0, N, nBins + 1));

x = nan(nBins, 1);
p = nan(nBins, 1);
sem = nan(nBins, 1);
nPerBin = nan(nBins, 1);

for b = 1:nBins
    ii = (binEdges(b) + 1):binEdges(b + 1);
    x(b) = mean(Ss(ii));
    p(b) = mean(correctSorted(ii));
    nPerBin(b) = numel(ii);
    sem(b) = sqrt(p(b) * (1 - p(b)) / nPerBin(b));
end

[B, ~, stats] = glmfit(S, isCorrect, 'binomial', 'link', 'logit');
xFit = linspace(min(S), max(S), 400)';
pFit = glmval(B, xFit, 'logit');

figure(figureNum);
clf;
hold on;
errorbar(x, p, sem, 'o', 'MarkerSize', 6, 'LineWidth', 1.2, 'CapSize', 0);
plot(x, p, '-', 'LineWidth', 1.0);
plot(xFit, pFit, '--', 'LineWidth', 1.5);
xlabel('Integrated noise evidence');
ylabel('P(correct)');
title(sprintf('%s: P(correct) vs integrated noise (%d bins, N = %d)', titleStr, nBins, N));
grid on;
box off;
hold off;

out = struct();
out.S = S;
out.isCorrect = isCorrect;
out.xBin = x;
out.pBin = p;
out.semBin = sem;
out.nPerBin = nPerBin;
out.beta = B;
out.betaSE = stats.se;
out.xFit = xFit;
out.pFit = pFit;
end