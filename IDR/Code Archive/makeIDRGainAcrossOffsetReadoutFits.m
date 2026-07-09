function fitSummary = makeIDRGainAcrossOffsetReadoutFits(varargin)
% makeIDRGainAcrossOffsetReadoutFits
%
% Adapter from IDR Weibull probe-gain offset summary to the generic
% fitAcrossOffsetReadout.m readout-function fitter.
%
% Required table variables:
%   probeDirDeg
%   medianGProbeOverGPref
%
% Preferred uncertainty variable:
%   bootstrapVarGProbeOverGPref
%
% Also accepted, if bootstrap variance is not yet present:
%   varGProbeOverGPref
%   semGProbeOverGPref
%   ci95LowGProbeOverGPref / ci95HighGProbeOverGPref
%   q25GProbeOverGPref / q75GProbeOverGPref  [fallback only]
%
% The 0-deg normalization anchor is added if absent:
%   measurement = 1
%   variance    = 0
%   nSessions   = Inf
% The 0-deg point is retained in the output but is ignored by
% fitAcrossOffsetReadout.m during fitting.

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'Bounds', struct(), @isstruct);
addParameter(p, 'SigmaMTDeg', 37.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'PhiDeg', -180:1:179, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)));


%% REMOVE THIS KLUDGE ONCE WE ARE EXCLUDING 1° PROBES
addParameter(p, 'MinProbeSessionsForFit', 2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);

parse(p, varargin{:});
opts = p.Results;

% -------------------------------------------------------------------------
% Load summary and find offsetSummary table
% -------------------------------------------------------------------------
gainSummaryFile = fullfile(domainFolder(mfilename('fullpath')), 'Data', 'AcrossOffsetSummaries', 'GainAnalysis', ...
  sprintf('IDRGainProbeOffsetSummary_%s.mat', opts.Animal));
L = load(gainSummaryFile);
offsetSummary = L.summaryData.offsetSummary;
vars = offsetSummary.Properties.VariableNames;

% mustHave = {'probeDirDeg', 'medianGProbeOverGPref'};
% for i = 1:numel(mustHave)
%     assert(ismember(mustHave{i}, vars), 'makeIDRGainAcrossOffsetReadoutFits:MissingVariable', ...
%         'offsetSummary is missing required variable: %s.', mustHave{i});
% end


% -------------------------------------------------------------------------
% Extract measurements
% -------------------------------------------------------------------------
offsetsDeg   = offsetSummary.probeDirDeg(:)';
measurements = offsetSummary.medianGProbeOverGPref(:)';
nSessions = offsetSummary.nProbeSessions(:)';
variances = offsetSummary.bootstrapVarGProbeOverGPref(:)';

% Number of contributing probe sessions, if available.
% if ismember('nProbeSessions', vars)
    % nSessions = offsetSummary.nProbeSessions(:)';
% elseif ismember('nSessions', vars)
%     nSessions = offsetSummary.nSessions(:)';
% else
%     nSessions = ones(size(offsetsDeg));
% end

% -------------------------------------------------------------------------
% Extract or infer variance
% -------------------------------------------------------------------------
% if ismember('bootstrapVarGProbeOverGPref', vars)
%     variances = offsetSummary.bootstrapVarGProbeOverGPref(:)';
% 
% elseif ismember('varGProbeOverGPref', vars)
%     variances = offsetSummary.varGProbeOverGPref(:)';
% 
% elseif ismember('semGProbeOverGPref', vars)
%     variances = offsetSummary.semGProbeOverGPref(:)'.^2;
% 
% elseif all(ismember({'ci95LowGProbeOverGPref', ...
%                      'ci95HighGProbeOverGPref'}, vars))
%     ciLow  = offsetSummary.ci95LowGProbeOverGPref(:)';
%     ciHigh = offsetSummary.ci95HighGProbeOverGPref(:)';
%     variances = ((ciHigh - ciLow) ./ (2 * 1.96)).^2;
% 
% elseif all(ismember({'q25GProbeOverGPref', ...
%                      'q75GProbeOverGPref'}, vars))
%     % Fallback only. This treats the interquartile range as if it came from
%     % an approximately Gaussian bootstrap distribution.
%     q25 = offsetSummary.q25GProbeOverGPref(:)';
%     q75 = offsetSummary.q75GProbeOverGPref(:)';
%     bootstrapSDApprox = (q75 - q25) ./ 1.349;
%     variances = bootstrapSDApprox.^2;
% 
%     warning('makeIDRGainAcrossOffsetReadoutFits:ApproxVariance', ...
%         ['Using q25/q75 to approximate bootstrap variance. ' ...
%          'Prefer saving bootstrapVarGProbeOverGPref directly.']);
% 
% else
%     error('makeIDRGainAcrossOffsetReadoutFits:MissingVariance', ...
%         ['Could not find bootstrap variance or a recognized uncertainty ' ...
%          'field for medianGProbeOverGPref.']);
% end

% -------------------------------------------------------------------------
% Add or validate 0-deg normalization anchor
% -------------------------------------------------------------------------
anchorIdx = find(abs(offsetsDeg) < 1e-9);

if isempty(anchorIdx)
    offsetsDeg   = [0, offsetsDeg];
    measurements = [1, measurements];
    variances    = [0, variances];
    nSessions    = [Inf, nSessions];

elseif isscalar(anchorIdx)
    % Treat the anchor as analytical. The fitter excludes it, but keeping it
    % exactly normalized avoids ambiguity in downstream plots/summaries.
    measurements(anchorIdx) = 1;
    variances(anchorIdx)    = 0;

else
    error('makeIDRGainAcrossOffsetReadoutFits:DuplicateAnchor', ...
        'More than one 0-deg row found in offsetSummary.');
end

% Sort for clean summaries and plots.
[offsetsDeg, sortIdx] = sort(offsetsDeg);
measurements = measurements(sortIdx);
variances    = variances(sortIdx);
nSessions    = nSessions(sortIdx);

% Basic validation. The 0-deg variance may be zero because it is not fit.
isAnchor = abs(offsetsDeg) < 1e-9;
badNonAnchor = ~isAnchor & ...
    (~isfinite(offsetsDeg) | ~isfinite(measurements) | ...
     ~isfinite(variances) | variances <= 0 | nSessions <= 0);

assert(~any(badNonAnchor), ...
    'makeIDRGainAcrossOffsetReadoutFits:BadInputRows', ...
    'One or more non-anchor offset rows have invalid measurement, variance, or nSessions.');

% Exclude offsets with too few contributing probe sessions from the readout fit,
% while preserving their measured values in the saved summary.
fitMeasurements = measurements;
fitVariances = variances;
excludeFromFit = ~isAnchor & nSessions < opts.MinProbeSessionsForFit;
fitMeasurements(excludeFromFit) = NaN;
fitVariances(excludeFromFit) = NaN;

% Call generic readout fitter
fitSummary = fitAcrossOffsetReadout(fitMeasurements, fitVariances, offsetsDeg, ...
    'NSessions', nSessions, 'Bounds', opts.Bounds, 'SigmaMTDeg', opts.SigmaMTDeg, ...
    'PhiDeg', opts.PhiDeg, 'SourceMeasureType', 'WeibullGain_gProbeOverGPref', ...
    'SourceSideType', 'probe', 'SourceStepType', 'INC', 'SourceMode', 'sessionBootstrapMedian');

fitSummary.type = 'Beta';
fitSummary.source = struct();
fitSummary.source.gainSummaryFile = gainSummaryFile;
fitSummary.source.rawMeasurements = struct( ...               %% this one entry is a kludge to be removed
  'offsetsDeg', offsetsDeg, ...
  'values', measurements, ...
  'variances', variances, ...
  'nSessions', nSessions, ...
  'excludeFromFit', excludeFromFit, ...
  'exclusionReason', sprintf('nSessions < %g', opts.MinProbeSessionsForFit));
fitSummary.source.adapterFunction = mfilename;
fitSummary.source.offsetSummaryVariable = 'medianGProbeOverGPref';
fitSummary.source.varianceInterpretation = ...
  'Variance of bootstrap distribution for median gProbe/gPrefParent by offset.';
fitSummary.source.anchorNote = ...
  '0-deg anchor set to measurement=1 and excluded from DOG fitting by fitAcrossOffsetReadout.';

% -------------------------------------------------------------------------
% Add compatibility fields for makeAcrossOffsetPlots.
% These fields mimic the historical acrossOffsetSummary layout used by
% updateAcrossOffsetSummaries, but they are derived from the gain-summary
% offset table rather than from kernel bootstraps.
% -------------------------------------------------------------------------

raw = fitSummary.source.rawMeasurements;

plotMask = ~isAnchor & ~excludeFromFit & isfinite(raw.values) & isfinite(raw.variances) & raw.variances > 0;
plotOffsets = raw.offsetsDeg(plotMask);
plotValues  = raw.values(plotMask);
plotVars    = raw.variances(plotMask);
plotNSess   = raw.nSessions(plotMask);

% Use bootstrap percentile CI fields from offsetSummary if present.
ciLow = nan(size(plotOffsets));
ciHigh = nan(size(plotOffsets));

if all(ismember({'ciLowGProbeOverGPref','ciHighGProbeOverGPref'}, vars))
  [tfCI, locCI] = ismember(plotOffsets, offsetSummary.probeDirDeg(:)');
  ciLow(tfCI)  = offsetSummary.ciLowGProbeOverGPref(locCI(tfCI));
  ciHigh(tfCI) = offsetSummary.ciHighGProbeOverGPref(locCI(tfCI));
else
  % Fallback: approximate 95% CI from bootstrap variance.
  ciLow  = plotValues - 1.96 .* sqrt(plotVars);
  ciHigh = plotValues + 1.96 .* sqrt(plotVars);
end

% Trial counts are not used deeply by the plots except for annotation.
% Preserve total INC trials where available, placing them in column 2 because
% stepType is INC.
nTrials = nan(size(plotOffsets));
if ismember('nIncTrials', vars)
  [tfN, locN] = ismember(plotOffsets, offsetSummary.probeDirDeg(:)');
  nTrials(tfN) = offsetSummary.nIncTrials(locN(tfN));
end

empirical = repmat(struct( ...
  'probeOffsetDeg', NaN, ...
  'stepType', 2, ...
  'pooledScale', NaN, ...
  'meanSessionScale', NaN, ...
  'medianSessionScale', NaN, ...
  'semSessionScale', NaN, ...
  'boot95', [NaN NaN], ...
  'nSessions', 0, ...
  'nTrials', NaN, ...
  'nTrialsByStep', [NaN NaN], ...
  'sessionWeighting', 'session_bootstrap_median_gain_ratio', ...
  'distributionSkew', NaN, ...
  'meanScale', NaN, ...
  'medianScale', NaN, ...
  'semScale', NaN), numel(plotOffsets), 1);

offsetBootstrap = repmat(struct( ...
  'probeOffsetDeg', NaN, ...
  'bootScale', [], ...
  'bootMean', NaN, ...
  'bootMedian', NaN, ...
  'boot68', [NaN NaN], ...
  'boot95', [NaN NaN]), numel(plotOffsets), 1);

for k = 1:numel(plotOffsets)
  empirical(k).probeOffsetDeg = plotOffsets(k);
  empirical(k).stepType = 2;  % INC
  empirical(k).pooledScale = plotValues(k);
  empirical(k).meanSessionScale = NaN;
  empirical(k).medianSessionScale = plotValues(k);
  empirical(k).semSessionScale = sqrt(plotVars(k));
  empirical(k).boot95 = [ciLow(k), ciHigh(k)];
  empirical(k).nSessions = plotNSess(k);
  empirical(k).nTrials = nTrials(k);
  empirical(k).nTrialsByStep = [NaN, nTrials(k)];

  % Backward-compatible aliases.
  empirical(k).meanScale = empirical(k).pooledScale;
  empirical(k).medianScale = empirical(k).medianSessionScale;
  empirical(k).semScale = empirical(k).semSessionScale;

  offsetBootstrap(k).probeOffsetDeg = plotOffsets(k);
  offsetBootstrap(k).bootScale = [];
  offsetBootstrap(k).bootMean = NaN;
  offsetBootstrap(k).bootMedian = plotValues(k);
  offsetBootstrap(k).boot68 = [NaN NaN];
  offsetBootstrap(k).boot95 = [ciLow(k), ciHigh(k)];
end

fitSummary.empirical = empirical;

fitSummary.bootstrap = struct();
fitSummary.bootstrap.offsetBootstrap = offsetBootstrap;
fitSummary.bootstrap.fitBootstrap = struct( ...
  'modelName', 'gain_offset_variance_only', ...
  'offsetsDeg', plotOffsets, ...
  'offsetFitVar', plotVars, ...
  'offsetFitWeights', 1 ./ plotVars, ...
  'fitWeights', 1 ./ plotVars, ...
  'note', ['Compatibility field for makeAcrossOffsetPlots. Variances are ' ...
           'bootstrap variances of offset-level median gProbe/gPrefParent.']);

fitSummary.meta = struct();
fitSummary.meta.analysisDate = datetime('now');
fitSummary.meta.scaleMetric = 'WeibullGain_gProbeOverGPref';
fitSummary.meta.ciType = 'session_bootstrap_percentile';
fitSummary.meta.nBoot = NaN;
fitSummary.meta.bootstrapType = 'probe_session_bootstrap';
fitSummary.meta.angleUnits = 'deg';
fitSummary.meta.normalization = 'gProbeOverGPrefParent_with_0deg_anchor';
fitSummary.meta.offsetKeysDeg = plotOffsets;
fitSummary.meta.notes = 'Plot-compatible summary generated by makeIDRGainAcrossOffsetReadoutFits.';

%% It is at this point that I would like to use makeAcrossOffsetPlots to plot the fit
opts.PlotDir = fullfile(domainFolder(mfilename('fullpath')), 'Plots', 'AcrossProbes', 'ReadoutFits', 'Beta');
opts.NBoot = max(offsetSummary.nBootstrap, [], 'omitnan');
makeAcrossOffsetPlots(fitSummary, opts);

outputFile = fullfile(domainFolder(mfilename('fullpath')), 'Data', 'AcrossOffsetSummaries', 'GainAnalysis', ...
  sprintf('IDRGainAcrossOffsetReadoutFits_%s.mat', opts.Animal));
save(outputFile, 'fitSummary');

end