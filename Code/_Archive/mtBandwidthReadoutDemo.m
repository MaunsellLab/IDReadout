function mtBandwidthReadoutDemo(varargin)
% mtBandwidthReadoutDemo
%
% Explore how pooled precision (SNR) varies as a function of the bandwidth
% of MT direction preferences sampled around the drift direction.
%
% Assumptions:
%   - MT neurons are independent Poisson units
%   - Direction tuning is Gaussian in preferred-direction offset
%   - Peak mean response is normalized to 1 for all baseline conditions
%   - Pooling is over neurons with preferred directions within +/- bandwidth
%
% Precision metric:
%   SNR = sum(w .* mu) / sqrt(sum((w.^2) .* mu))
%
% Pooling models:
%   1) Equal weights within band
%   2) Graded positive weights within band (same Gaussian form as tuning)
%   3) Optimal weights within band for estimating the scalar "motion strength"
%      represented by the Gaussian tuning profile under Poisson noise
%
% Usage:
%   mtBandwidthReadoutDemo
%   mtBandwidthReadoutDemo('sigmaDeg', 40, 'dPhiDeg', 1, 'maxBandwidthDeg', 180)
%
% Notes:
%   For the optimal model here, with independent Poisson noise and scalar
%   signal profile s(phi) = g(phi), the optimal linear weights are
%       w_opt(phi) proportional to s(phi) / mu(phi)
%   within the included band, and 0 outside.
%
% John Maunsell discussion version

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
p.addParameter('sigmaDeg', 40, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('dPhiDeg', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('maxBandwidthDeg', 180, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('r0List', [0 0.25], @(x) isnumeric(x) && isvector(x));
p.addParameter('gradedSigmaDeg', [], @(x) isempty(x) || (isscalar(x) && x > 0));
p.parse(varargin{:});

sigmaDeg        = p.Results.sigmaDeg;
dPhiDeg         = p.Results.dPhiDeg;
maxBandwidthDeg = p.Results.maxBandwidthDeg;
r0List          = p.Results.r0List;
gradedSigmaDeg  = p.Results.gradedSigmaDeg;

if isempty(gradedSigmaDeg)
    gradedSigmaDeg = sigmaDeg;  % default heuristic graded weights
end

% -------------------------------------------------------------------------
% Direction preference grid and bandwidths
% -------------------------------------------------------------------------
phiDeg = -maxBandwidthDeg:dPhiDeg:maxBandwidthDeg;
phiRad = deg2rad(phiDeg); %#ok<NASGU> % retained in case you later want cosine/projection models
bandwidthDeg = 0:dPhiDeg:maxBandwidthDeg;

nPhi = numel(phiDeg);
nBW  = numel(bandwidthDeg);
nR0  = numel(r0List);

% Gaussian tuning profile centered on drift direction (offset = 0)
g = exp(-(phiDeg.^2) / (2 * sigmaDeg^2));

% Heuristic graded positive weighting profile
gW = exp(-(phiDeg.^2) / (2 * gradedSigmaDeg^2));

% Density of neurons per degree is constant, so sums approximate integrals.
% Since dPhi is constant, it cancels in normalized comparisons; no need to
% multiply explicitly by dPhiDeg unless absolute scaling is needed.

modelNames = {'Equal weights', 'Graded positive', 'Optimal'};
nModels = numel(modelNames);

snrAbs  = nan(nR0, nModels, nBW);
snrNorm = nan(nR0, nModels, nBW);

% -------------------------------------------------------------------------
% Compute SNR curves
% -------------------------------------------------------------------------
for iR = 1:nR0
    r0 = r0List(iR);

    % Mean response profile, peak-normalized to 1
    mu = r0 + (1 - r0) * g;

    % Signal profile for estimating motion strength in the drift direction.
    % Here we take it as the tuning-shaped component.
    s = g;

    for iB = 1:nBW
        B = bandwidthDeg(iB);
        inBand = abs(phiDeg) <= B;

        % Ensure at least one sampled bin
        if ~any(inBand)
            continue;
        end

        % -----------------------------
        % Model 1: equal weights
        % -----------------------------
        w1 = double(inBand);

        % -----------------------------
        % Model 2: graded positive weights
        % -----------------------------
        w2 = gW .* inBand;

        % -----------------------------
        % Model 3: optimal weights
        % For independent Poisson noise and scalar signal profile s(phi),
        % optimal linear weights are proportional to s(phi)/mu(phi).
        % -----------------------------
        w3 = (s ./ mu) .* inBand;

        W = {w1, w2, w3};

        for iM = 1:nModels
            w = W{iM};

            % pooled mean signal and Poisson variance
            num = sum(w .* s);
            den = sqrt(sum((w .^ 2) .* mu));

            if den > 0
                snrAbs(iR, iM, iB) = num / den;
            else
                snrAbs(iR, iM, iB) = NaN;
            end
        end
    end

    % Normalize each model by its zero-bandwidth value
    for iM = 1:nModels
        ref = snrAbs(iR, iM, 1);
        snrNorm(iR, iM, :) = snrAbs(iR, iM, :) / ref;
    end
end

% -------------------------------------------------------------------------
% Plot normalized SNR
% -------------------------------------------------------------------------
figure(1); clf;

tiledlayout(1, nR0, 'TileSpacing', 'compact', 'Padding', 'compact');

for iR = 1:nR0
    nexttile;
    hold on;

    for iM = 1:nModels
        plot(bandwidthDeg, squeeze(snrNorm(iR, iM, :)), 'LineWidth', 2);
    end

    yline(1, '--k');
    xlabel('Bandwidth around drift direction (deg)');
    ylabel('Normalized pooled SNR');
    title(sprintf('r_0 = %.2f, \\sigma = %.1f^\\circ', r0List(iR), sigmaDeg));
    legend(modelNames, 'Location', 'best');
    box off;
end

% -------------------------------------------------------------------------
% Plot absolute SNR
% -------------------------------------------------------------------------
figure(2); clf;

tiledlayout(1, nR0, 'TileSpacing', 'compact', 'Padding', 'compact');

for iR = 1:nR0
    nexttile;
    hold on;

    for iM = 1:nModels
        plot(bandwidthDeg, squeeze(snrAbs(iR, iM, :)), 'LineWidth', 2);
    end

    xlabel('Bandwidth around drift direction (deg)');
    ylabel('Absolute pooled SNR');
    title(sprintf('r_0 = %.2f, \\sigma = %.1f^\\circ', r0List(iR), sigmaDeg));
    legend(modelNames, 'Location', 'best');
    box off;
end

% -------------------------------------------------------------------------
% Report bandwidth of peak normalized SNR for each condition/model
% -------------------------------------------------------------------------
fprintf('\nPeak normalized SNR bandwidths:\n');
for iR = 1:nR0
    fprintf('  r0 = %.2f\n', r0List(iR));
    for iM = 1:nModels
        y = squeeze(snrNorm(iR, iM, :));
        [peakVal, idx] = max(y);
        fprintf('    %-16s peak = %.4f at %3d deg\n', ...
            modelNames{iM}, peakVal, bandwidthDeg(idx));
    end
end

% -------------------------------------------------------------------------
% Optional third figure: show example tuning and weight profiles
% at a representative bandwidth
% -------------------------------------------------------------------------
repB = min(60, maxBandwidthDeg);
inBand = abs(phiDeg) <= repB;

figure(3); clf;
tiledlayout(1, nR0, 'TileSpacing', 'compact', 'Padding', 'compact');

for iR = 1:nR0
    r0 = r0List(iR);
    mu = r0 + (1 - r0) * g;
    s  = g;

    w1 = double(inBand);
    w2 = gW .* inBand;
    w3 = (s ./ mu) .* inBand;

    nexttile;
    yyaxis left
    plot(phiDeg, mu, 'k', 'LineWidth', 2);
    ylabel('Mean response');

    yyaxis right
    plot(phiDeg, w1, 'LineWidth', 1.5); hold on;
    plot(phiDeg, w2, 'LineWidth', 1.5);
    plot(phiDeg, w3, 'LineWidth', 1.5);
    ylabel('Weight');
    xlabel('Preferred direction offset (deg)');
    title(sprintf('Profiles at B = %d^\\circ, r_0 = %.2f', repB, r0));
    legend({'\mu(\phi)', 'Equal', 'Graded', 'Optimal'}, 'Location', 'best');
    box off;
end
end