# MT max-pooling simulations

## Split-signal experimental-design module (August 31, 2026)

`runSplitSignalDesignSweep` evaluates the proposed two-direction split-signal
experiment, including 90- versus 120-degree separation, prestep pedestal size,
an optional vector-balancing prestep stream, step-only independent noise,
adaptation sensitivity, losing-stream gains, model likelihoods, and session
requirements. See `README_SPLIT_SIGNAL.md` and run
`validateSplitSignalSimulation` before a full sweep.

`runSplitSignalModelRecovery` performs the revised, fair comparison at
90-degree separation with full prestep pedestals and no balancing stream. Max,
population-vector magnitude, and population squared-energy models share the
same MT front end and correlated pre-pooling noise; thresholds are recalibrated,
decision scale and max softness are refitted, BIC corrects for the max family's
extra parameter, and every model is used in turn to generate data.
See `README_SPLIT_SIGNAL_RECOVERY.md` and run
`validateSplitSignalModelRecovery` first.

This MATLAB package implements educational stages of the unified IDR/IDQ/IDS
max-pooling analysis: an evenly spaced aligned `n`-sweep, a dense-bank
paired-probe change-detection sweep, a correlation-length sweep, and a
response-dependent pooled-Poisson variance stage. It also includes temporal-max
and drift-template stages retained as historical comparisons, and a current
drift-cancellation plus firing-dependent adaptation stage.

Run:

```matlab
cd MaxPoolingSimulation
results = runMaxPoolingNSweep;
pairedResults = runPairedProbeOffsetSweep;
correlationResults = runPairedProbeCorrelationSweep;
poissonResults = runPairedProbePoissonSweep;
temporalMaxResults = runPairedProbeTemporalMaxSweep;
driftTemplateResults = runPairedProbeDriftTemplateSweep;
adaptationResults = runPairedProbeAdaptationSweep;
```

The headers in the primary run files are the authoritative model
specifications. They describe the equations, noise placement, calibration,
gain recovery, diagnostics, and later planned stages.

## Sparse aligned n-sweep

- `n = 2:20` evenly spaced direction mechanisms by default.
- One zero-mean scalar noise integral is centered exactly on each mechanism.
- One mechanism/stream is aligned with the drift direction.
- Every stream affects every mechanism through a centered Gaussian MT tuning
  function with default sigma 37.5 degrees.
- Independent Gaussian mechanism noise is added after step integration and
  before hard-max pooling.
- No downstream noise or lapses. A common pedestal is available as an
  invariance control; it adds the same response to every evenly spaced
  mechanism and cannot affect choices.
- Fixed-step and separately 75%-matched simulations.
- Simultaneous probit recovery of all stream gains.
- Optional close-to-max counterfactual diagnostic.

## Dense paired-probe sweep

- Dense 1-degree direction bank with smooth circular internal noise.
- Three imposed streams: drift and probes at `+/-offset`.
- Every mechanism receives its tuning-weighted mean step coherence minus its
  tuning-weighted mean prestep coherence before max pooling.
- Constant base drift and probe pedestals cancel exactly in this first-pass
  linear change detector.
- Gaussian stream variables represent scalar changes in epoch-mean coherence;
  there is no frame-resolved temporal simulation.
- Independent and yoked probe conditions.
- A single internal-noise SD is calibrated at a reference condition so a fixed
  20% drift increment gives 75% correct, then held fixed everywhere.
- All offsets are analyzed at the same 20% increment; performance and the
  independently measured 75% threshold are both reported.
- Fixed-winner versus winner-switch decompositions.
- Winner categories for drift, probe, and other/intermediate mechanisms.

`BaseDriftCoherence`, `DriftPedestal`, and `ProbePedestal` remain explicit
configuration fields as invariance controls. Changing them cannot affect this
stage because the same constants appear in the two epoch means and subtract
away before the direction-mechanism maximum.

## Internal-noise correlation-length sweep

`runPairedProbeCorrelationSweep` repeats the delta-coherence paired-probe stage
over nominal circular correlation lengths of 10, 15, 20, 30, 45, 60, and 90
degrees. Circular-mean removal remains enabled so the sweep changes only the
smoothness/dimensionality parameter within the current noise model.

At each length, the internal-noise SD is recalibrated once at the common 20%
increment/75% reference condition and then held fixed across all offsets and
yoking modes. The output reports both the nominal correlation length and the
effective covariance rank after circular-mean removal. Default effective ranks
range from approximately 19.5 at 10 degrees to 2.4 at 90 degrees.

The sweep uses a reduced default offset set and 12,000 analysis trials to keep
runtime manageable. Modify `pairedProbeCorrelationSweepDefaults.m` for a full
20,000-trial replication.

## Pooled-Poisson variance sweep

`runPairedProbePoissonSweep` replaces the homoscedastic internal-noise SD with
the variance expected from prestep and step spike counts. Absolute firing rates
are set by a random-motion baseline plus coherence-dependent MT response, but
the max pool still receives the step-minus-prestep response difference.

One effective independent-neuron pool size is calibrated at the 7.071% base,
90-degree independent-probe reference condition. It is then held fixed across
7.071% and 40% base drift, all offsets, and both yoking conditions. This makes
the cross-base comparison a model prediction. The output includes performance
at a fixed 20% increment, separately measured 75% thresholds, recovered gains,
winner locations, absolute-rate diagnostics, and clamp fractions.

The dense bank does not receive independent noise at every degree. Its
standardized count fields retain the 30-degree circular correlation structure;
Poisson marginal SD is then applied locally according to each mechanism's
absolute rate and epoch duration.

## Separate prestep and step maxima

This stage is retained for comparison but is no longer the active model.  Its
winner-direction diagnostic showed that a 40% base drift effectively eliminates
120-degree probe winners, predicting stronger negative probe gain in the
high-base condition - opposite to IDR.

`runPairedProbeTemporalMaxSweep` applies the direction maximum separately to
the absolute prestep and step response landscapes and then subtracts the two
scalar maxima within each spatial patch. This differs from taking the maximum
after subtracting each mechanism's prestep response. Consequently, constant
baseline activity cancels when the same mechanism wins both epochs, while the
drift and probe pedestals continue to determine competitive gaps and temporal
winner changes.

The first pass uses one integrated sample per epoch, homoscedastic smooth
Gaussian mechanism noise, and duration scaling from the 250-ms step to the
750-ms prestep. One step-noise SD is calibrated at the low-base 90-degree
independent condition and held fixed across both base drifts, every offset, and
both yoking modes. Outputs include full 3-by-3 prestep-to-step winner-transition
matrices and a decomposition of stream influence by temporal winner identity.
The temporal-max stage also saves the complete binned distribution of step
winners across direction, together with the signed drift and per-stream probe
sensitivity in each bin and their winner-probability-weighted contributions.
This makes it possible to see which angular regions drive a recovered gain
positive or negative instead of relying only on drift/probe/other categories.

## Perceptual drift-template cancellation (intermediate comparison)

`runPairedProbeDriftTemplateSweep` removes an
explicit fraction of the perceived sustained drift profile before direction
max pooling:

```text
E_m = absolute step response_m - alpha * drift baseline * drift profile_m
```

There is no prestep maximum and no difference of scalar epoch maxima.  Probe
pedestals remain in the step landscape, as do the 20% drift increment and the
scalar drift/probe coherence fluctuations.  The cancellation is additive, not
divisive, and therefore does not impose Weber scaling.

The default focused sweep uses base drifts 7.071% and 40%, offsets 120, 170,
and 180 degrees, independent and yoked probes, and cancellation fractions from
0 to 1.  One internal-noise SD is calibrated at full cancellation in the
low-base 90-degree reference condition and then held fixed.  Outputs show
recovered gains, performance, separately measured thresholds, winner-category
fractions, and complete winner-direction distributions and signed
contributions.  The main question is whether one cancellation range can retain
negative low-base 120-degree gain while restoring high-base far-probe winners
and preserving the 170-180-degree rebound.

## Drift cancellation plus firing-dependent adaptation

`runPairedProbeAdaptationSweep` is the current stage. It cancels the perceived
sustained drift exactly, retains both probe pedestals, and then applies a
nonnegative suppressive field centered on the prestep drift direction. The
field is motivated by the preferred/null MT response displacement after a
strong prestep. Its strength is parameterized as a coherence-equivalent
aftereffect at a 60% prestep and scales linearly with the actual base drift:

```text
A(base) = AftereffectAt60 * (BaseDriftCoherence + DriftPedestal) / 60
```

The default adaptation width equals the 37.5-degree MT tuning sigma. Far
mechanisms are spared rather than explicitly enhanced. The max can nevertheless
give them a relative advantage because drift-sensitive mechanisms are
selectively suppressed. The default sweep compares aftereffects from 0 through
8% in 1% increments, base drifts 7.071 and 40%, the five current IDS offsets
10, 20, 45, 120, and 170 degrees, and independent versus yoked probes.

One internal-noise SD is calibrated so the reliable 17-session IDS mean
threshold is reproduced at the low-base, 120-degree, independent, zero-
adaptation reference: 75% correct for a 15.3% step. That SD is held fixed
across every offset, base, yoking mode, and adaptation strength. The 170-degree
IDS mean (11.2%) is displayed as a three-session preliminary estimate rather
than used for calibration. The empirical-threshold figure also shows the IDS
means [13.0 11.9 16.2 15.3 11.2] and the provisional 21-22% IDR range.

Outputs include fixed-step performance, separately measured 75% thresholds,
recovered gains, fixed-winner/winner-switch decompositions, winner categories,
complete angular winner distributions with signed contributions, and an
explicit matched-condition table of 40%-minus-7.071% threshold changes. A
roughly 5-coherence-point threshold increase is the provisional cross-
experiment target. This directly tests whether measured adaptation can restore
high-base 120-degree probe wins while producing a plausible threshold cost.

## Useful modifications

```matlab
% Faster smoke run
cfg = maxPoolingDefaults;
cfg.NValues = 2:8;
cfg.NTrials = 10000;
cfg.NCalibrationTrials = 10000;
results = runMaxPoolingNSweep(cfg);

% Explore more internal mechanism noise
results = runMaxPoolingNSweep('InternalNoiseSD', 2);

% Save elsewhere and suppress visible figures
results = runMaxPoolingNSweep('OutputDir', '/path/to/output', ...
    'Visible', 'off');

% Selected paired-probe offsets
pcfg = pairedProbeDefaults;
pcfg.OffsetsDeg = [10 25 45 120 170 179];
pairedResults = runPairedProbeOffsetSweep(pcfg);

% Noise-correlation sweep
ccfg = pairedProbeCorrelationSweepDefaults;
ccfg.CorrelationLengthsDeg = [15 30 60];
correlationResults = runPairedProbeCorrelationSweep(ccfg);

% Pooled-Poisson stage with one calibrated pool size held across base drift
poissonCfg = pairedProbePoissonDefaults;
poissonCfg.BaseDriftCoherences = [7.071 40];
poissonResults = runPairedProbePoissonSweep(poissonCfg);

% Max each epoch first, then form the temporal difference
temporalCfg = temporalMaxDefaults;
temporalResults = runPairedProbeTemporalMaxSweep(temporalCfg);

% Intermediate comparison: remove a fraction of the drift template
templateCfg = driftTemplateDefaults;
templateResults = runPairedProbeDriftTemplateSweep(templateCfg);

% Current model: exact drift cancellation plus firing-dependent adaptation
adaptationCfg = adaptationSweepDefaults;
adaptationResults = runPairedProbeAdaptationSweep(adaptationCfg);

% Use a supplied internal-noise SD without the one-time calibration
pcfg.CalibrateInternalNoise = false;
pcfg.InternalNoiseSD = 8;
pairedResults = runPairedProbeOffsetSweep(pcfg);
```

The output directories contain tabular summaries, operating points, complete
MAT results, and stage-specific PDF figures.

## Files

- `runMaxPoolingNSweep.m`: primary entry point and complete model description.
- `maxPoolingDefaults.m`: editable default configuration.
- `centeredMTSensitivity.m`: circular centered MT sensitivity.
- `makeSparseTrialBank.m`: common-random-number trial generation.
- `calibrateMaxPoolingStep.m`: 75% step calibration.
- `simulateSparseMaxPooling.m`: hard-max choices and gain recovery.
- `fitProbitIRLS.m`: toolbox-independent probit regression.
- `computeCloseToMaxDiagnostic.m`: paired max-gap perturbation analysis.
- `plotSparseNSweepResults.m`: figures.
- `saveSimulationFigures.m`: PDF export helper.
- `validateSparseMaxPooling.m`: short invariant/calibration self-check.
- `runPairedProbeOffsetSweep.m`: paired-probe entry point and model description.
- `pairedProbeDefaults.m`: paired-probe configuration and offset grid.
- `makeDensePairedTrialBank.m`: common random inputs and dense bank.
- `generateCircularGPNoise.m`: smooth circular internal-noise fields.
- `circularNoiseEffectiveRank.m`: implemented covariance participation ratio.
- `buildPairedProbeResponses.m`: delta-coherence response landscapes.
- `calibratePairedProbeNoise.m`: one-time internal-noise calibration at the
  fixed reference increment.
- `calibratePairedProbeStep.m`: fixed-noise measurement of each 75% threshold.
- `simulatePairedProbeCondition.m`: choices, gains, winners, decomposition.
- `decomposeStreamInfluence.m`: fixed-winner/winner-switch effects.
- `classifyDenseWinners.m`: winner-location diagnostics.
- `plotPairedProbeSweep.m`, `savePairedProbeFigures.m`: paired-probe figures.
- `validatePairedProbeSimulation.m`: paired-probe self-check.
- `runPairedProbeCorrelationSweep.m`: correlation-length sweep entry point and
  complete stage description.
- `pairedProbeCorrelationSweepDefaults.m`: correlation-sweep configuration.
- `plotPairedProbeCorrelationSweep.m`: combined sweep figures.
- `validatePairedProbeCorrelationSweep.m`: short sweep self-check.
- `runPairedProbePoissonSweep.m`: pooled-Poisson entry point and complete
  model specification.
- `pairedProbePoissonDefaults.m`: rate, duration, base-drift, and output
  configuration.
- `makePoissonPairedTrialBank.m`: common stream draws and correlated prestep
  and step count fields.
- `preparePoissonPairedProbeCondition.m`: absolute rates, exact mean changes,
  and response-dependent count variance.
- `calibratePoissonPoolSize.m`: one-time effective-pool calibration.
- `calibratePoissonPairedProbeStep.m`: fixed-pool 75% threshold measurement.
- `simulatePoissonPairedProbeCondition.m`: choices, gains, and winners.
- `plotPairedProbePoissonSweep.m`, `savePairedProbePoissonFigures.m`: figures.
- `validatePairedProbePoissonSweep.m`: short pooled-Poisson self-check.
- `runPairedProbeTemporalMaxSweep.m`: separate-epoch-max entry point and
  complete model specification.
- `temporalMaxDefaults.m`: temporal duration, base drift, and output settings.
- `makeTemporalMaxTrialBank.m`, `buildTemporalMaxResponses.m`: separate epoch
  noise fields and absolute response landscapes.
- `calibrateTemporalMaxNoise.m`, `calibrateTemporalMaxStep.m`: one-time noise
  calibration and fixed-noise threshold measurement.
- `simulateTemporalMaxCondition.m`: difference-of-max choices, gains, winner
  identities, and transition matrices.
- `decomposeTemporalMaxInfluence.m`: same/different temporal winner and
  perturbation-switch contributions.
- `plotPairedProbeTemporalMaxSweep.m`, `savePairedProbeTemporalMaxFigures.m`:
  stage figures.
- `validatePairedProbeTemporalMaxSweep.m`: short temporal-max self-check.
- `runPairedProbeDriftTemplateSweep.m`: intermediate drift-template entry point and
  complete model specification.
- `driftTemplateDefaults.m`: base drifts, cancellation fractions, focused
  offsets, and output configuration.
- `buildDriftTemplateResponses.m`: residual absolute step landscapes after
  selective drift-template cancellation.
- `calibrateDriftTemplateNoise.m`, `calibrateDriftTemplateStep.m`: one-time
  noise calibration and fixed-noise threshold measurement.
- `simulateDriftTemplateCondition.m`: choices, recovered gains, decompositions,
  and winner diagnostics.
- `computeWinnerDirectionDiagnostic.m`: shared binned winner distributions and
  signed winner-weighted contributions.
- `plotPairedProbeDriftTemplateSweep.m`,
  `savePairedProbeDriftTemplateFigures.m`: intermediate-stage figures.
- `validatePairedProbeDriftTemplateSweep.m`: cancellation and base-invariance
  checks for the intermediate drift-template stage.
- `runPairedProbeAdaptationSweep.m`: current entry point and complete model
  description for exact cancellation plus prestep firing-dependent adaptation.
- `adaptationSweepDefaults.m`: empirical IDS/IDR threshold constraints, base
  drifts, fine aftereffect sweep, adaptation width, and output settings.
- `buildAdaptationResponses.m`: retained probe-pedestal landscape after exact
  drift cancellation and the drift-centered suppressive adaptation field.
- `calibrateAdaptationNoise.m`, `calibrateAdaptationStep.m`: common noise
  calibration and condition-specific threshold measurement.
- `simulateAdaptationCondition.m`: gains, winner categories, decompositions,
  and angular signed contributions.
- `plotPairedProbeAdaptationSweep.m`, `savePairedProbeAdaptationFigures.m`:
  current-stage figures and PDF export.
- `validatePairedProbeAdaptationSweep.m`: short structural, base-invariance,
  and adaptation-direction self-check.

Before a long run, the optional self-check is:

```matlab
validateSparseMaxPooling
validatePairedProbeSimulation
validatePairedProbeCorrelationSweep
validatePairedProbePoissonSweep
validatePairedProbeTemporalMaxSweep
validatePairedProbeDriftTemplateSweep
validatePairedProbeAdaptationSweep
```
