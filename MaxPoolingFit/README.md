# MaxPoolingFit

MATLAB implementation of the first-fit dense-direction hard-max model for the
IDR, IDQ, and IDS choice data.

## Frozen model

For RDK `j` and direction mechanism `m`,

```text
H(theta-d) = G(theta-d) - mean_direction(G)
Rstep(j,m) = sum_s Cstep(j,s) H(theta_m - d(j,s))
Ppre(j,m)  = sum_s Cpre(j,s) H(theta_m - d(j,s))
E(j,m)     = Rstep(j,m) - a Ppre(j,m) + sigmaDomain * eta(j,m)
M(j)       = max_m E(j,m)
```

`G` is a circular Gaussian with fixed sigma 37.5 degrees. `H` is its
fixed-density random-dot redistribution response: a coherent dot replaces a
uniformly directed random dot rather than being added to an unchanged random
pool. The omitted isotropic baseline is identical across mechanisms and RDKs
and cancels from the decision. The same centering is applied to step, prestep,
and adaptation responses. `eta` is a
unit-marginal-variance periodic Gaussian field with fixed correlation length
30 degrees. The predicted choice is the RDK with the larger maximum. The primary
model uses one shared adaptation value, `a = 0.56`, selected by a 256-draw
profile likelihood and then held fixed. Adaptation subtracts a fraction of the
complete directional prestep population response. Drift and probe pedestals
therefore adapt nearby mechanisms symmetrically, and a strongly visible prestep
drift produces more absolute adaptation without a domain-specific rule.
Internal-noise SD is domain-specific (IDR, IDQ, or IDS) and shared across
animals within each domain.

The historical additive response remains available explicitly as
`CoherenceResponseMode="additive"`, but new fits default to
`"random-dot-redistribution"`. Because this correction changes deterministic
evidence, adaptation, noise, and correlation length must be re-estimated. It
does not require rebuilding the trial cache.

Reprofile the four design groups in a separate output folder:

```matlab
redistributionFolder = fullfile(outputFolder,"RandomDotRedistribution");

redistributionProfiles = runMaxPoolingDesignGroupProfiles( ...
    cachePath,groupedResult,redistributionFolder, ...
    'CoherenceResponseMode',"random-dot-redistribution", ...
    'AdaptationGrid',0:0.1:1,'NDraws',64);

redistributionResult = runMaxPoolingDesignGroupFit( ...
    cachePath,groupedResult,redistributionProfiles, ...
    redistributionFolder);
```

### Alternative directional-change decision rule

The primary v5 model compares the absolute adapted step maxima shown above. A
diagnostic v6 rule instead constructs direction-specific change evidence before
max pooling:

```text
Echange(j,m) = Rstep(j,m) - Rpre(j,m) - a Ppre(j,m)
Mchange(j)    = max_m Echange(j,m)
```

It then compares the two RDK change maxima. Because this is a different
evidence variable, its shared adaptation and domain-specific internal-noise SDs
must be re-estimated. Generate its initial profile with:

```matlab
changeProfile = runMaxPoolingAdaptationProfile( ...
    cachePath,result,outputFolder, ...
    'DecisionRule',"directional-change-max", ...
    'AdaptationGrid',0:0.1:1,'NDraws',64);
```

The profile outputs use the `MaxPoolingDirectionalChange` prefix and do not
replace the original adaptation-profile files.

The code uses exact hard maxima. Choice probabilities are numerical integrals
over fixed antithetic Gaussian-field draws. No lapse, decision bias, downstream
noise, eligibility field, or soft-max temperature is included.

### Low-confidence detection and guessing diagnostic

An optional diagnostic preserves the absolute adapted maximum in each RDK but
adds a no-detection state. Each noisy RDK maximum is centered on that RDK's
deterministic no-change maximum. If neither centered maximum reaches
`kappa*noiseSD`, the model guesses RDK 1 with probability 0.5; otherwise it
retains the original comparison of the two absolute adapted maxima.
`kappa=-Inf` is the exact original no-guess model.

Profile one shared `kappa` while reoptimizing domain-specific noise SDs:

```matlab
guessProfile = runMaxPoolingGuessCutoffProfile( ...
    cachePath,result,outputFolder, ...
    'KappaGrid',[-Inf 0:0.25:3], ...
    'NDraws',64);
```

The output includes domain-level likelihood, noise, accuracy, and guess
fractions at every cutoff, plus condition-level diagnostics at the best cutoff.

### Changed and unchanged coefficients by exact step

To test whether the RDK asymmetry is concentrated near threshold, estimate the
four noise coefficients independently within every domain-animal-offset-step
condition:

```matlab
stepCoefficients = runMaxPoolingStepNoiseCoefficients( ...
    cachePath,result,outputFolder,'MinimumTrials',100);
```

The same session-fixed-effect regression is applied to observed choices and
primary-model probabilities. Outputs report raw binomial coefficients per
percentage-point coherence, avoiding psychometric-reference normalization.
Plots show the changed coefficient and the sign-reversed unchanged coefficient
on the same axis. Cells with fewer than five sessions are retained descriptively
but marked `inferenceReliable=false`.

For sign-consistent cells (positive changed coefficient and negative unchanged
coefficient), the same run also computes the heuristic gain-implied change-RDK
win fraction

```text
changed / (changed + abs(unchanged))
```

and its delta-method standard error. The
`MaxPoolingGainImpliedWinFraction.pdf` figure plots this quantity against
accuracy for observed choices and primary-model probabilities, separately for
drift and probe noise. The dashed identity line is the prediction of a simple
exclusive winner-gating interpretation; it is a diagnostic rather than an
assumption of the fitted model. The corresponding weighted summaries are saved
in `MaxPoolingGainImpliedWinFractionAgreement.csv`.

Test whether changed- and unchanged-RDK gains have the same temporal shape
during the coherence step:

```matlab
temporalGain = runMaxPoolingTemporalGainDiagnostics( ...
    cachePath,result,outputFolder);
```

The default analysis selects IDR-Meetz 25-degree probe trials. It uses the
existing 30-sample framewise traces in the cache and therefore does not rebuild
the cache or refit the max-pooling model. The 250-ms step is reduced to five
fixed 50-ms bins. Changed-drift, unchanged-drift, changed-probe, and
unchanged-probe predictors from all five bins enter one joint fractional-
binomial regression with session fixed effects. Observed choices and primary-
model probabilities are analyzed identically. The unchanged coefficients are
sign-reversed for display, and a cluster-robust Wald test asks whether their
normalized temporal profile is proportional to the changed-RDK profile.
Session-bootstrap intervals are shown on the normalized profiles.

The run saves raw and normalized profile PDFs, long-form coefficient and shape-
test CSV files, and `MaxPoolingTemporalGainDiagnostics.mat`. Select a different
condition with `Domain`, `Animal`, `ProbeOffsetDeg`, or `StepMagnitudePC`.

Separate noise that determines the hard-max winner from errors introduced after
categorical selection:

```matlab
noisePlacement = runMaxPoolingNoisePlacementSweep( ...
    cachePath,result,outputFolder);
```

The default diagnostic uses the exact IDR-Meetz 25-degree, 23% increment
condition. It sweeps pre-max Gaussian mechanism noise from 0 to 24% coherence.
At each value it estimates changed/unchanged drift and probe coefficients, then
analytically solves the symmetric post-selection flip rate required to restore
observed accuracy. The flipped probabilities are passed back through the same
gain regression, so coefficient-ratio invariance is checked rather than
assumed. Values whose upstream accuracy is below observed accuracy are marked
infeasible because symmetric flips cannot improve performance.

Outputs include the complete sweep and observed targets as CSV and MAT files,
plus overview and absolute-gain PDFs. Ratio RMSE selects upstream noise using
both drift and probe ratios after the accuracy-matching flip; the corresponding
pre-flip ratio error is also retained. Absolute gains, log loss, and the
required flip rate remain independent tests of the resulting account.

Two ordered diagnostics test why the primary model overpredicts unchanged-RDK
gain. First, refit the explicit directional-change rule at its boundary optimum
`a = 0` and apply the identical exact-condition gain regression:

```matlab
decisionRuleGain = runMaxPoolingDecisionRuleGainDiagnostic( ...
    cachePath,result,outputFolder);
```

This reports the likelihood cost and domain noise estimates alongside observed,
absolute-step, and directional-change drift/probe coefficients for the default
IDR-Meetz 25-degree, 23% step condition.

Second, test whether directionally correlated internal variability is shared
between the two RDKs:

```matlab
sharedNoise = runMaxPoolingSharedRDKNoiseSweep( ...
    cachePath,result,outputFolder);
```

The first-pass grid uses cross-RDK correlations `0:0.25:1`, pre-max noise SDs
`[0 5 10 15 20 25]`, and 64 paired draws. At every grid point it solves the
downstream flip rate required to match observed accuracy and reruns the gain
regression. Paired draws are essential: all-cross-pairs integration assumes
independent RDK noise and would erase the imposed correlation. Refine the grid
and draw count only if the coarse sweep reveals a promising region.

## Entry points

```matlab
addpath('/path/to/MaxPoolingFit')
cfg = maxpoolfit.defaultConfig();
cache = maxpoolfit.loadCache('MaxPoolingTrialCache.mat');
fit = maxpoolfit.fitModel(cache, cfg);
prediction = maxpoolfit.predict(cache, fit, cfg);
```

The standard two-pass run saves a development fit, a higher-precision final
fit, fitted-trial predictions, group-resolved calibration diagnostics, and a
PDF comparing the development and final domain-noise estimates:

```matlab
result = runMaxPoolingFit('MaxPoolingTrialCache.mat','FitOutput');
```

To preserve the profile-likelihood provenance inside the saved primary result,
pass the refined profile when rerunning the fit:

```matlab
result = runMaxPoolingFit( ...
    'MaxPoolingTrialCache.mat','FitOutput', ...
    'Adaptation',0.56,'AdaptationProfile',adaptationProfileFine);
```

Profile attenuation of the complete prestep population representation while
reoptimizing the three domain-specific noise SDs at every value of `a`:

```matlab
adaptationProfile = runMaxPoolingAdaptationProfile( ...
    'MaxPoolingTrialCache.mat',result,'FitOutput', ...
    'AdaptationGrid',0:0.1:1,'NDraws',64);
```

This profile uses `Rstep-a*Ppre`: zero retains the complete step response,
whereas one subtracts the entire directional response already driven during
the prestep. `Ppre` is a vector over direction mechanisms, not a scalar.
It saves MAT, CSV, and delta-log-likelihood PDF outputs.

Assess calibration on the primary-fit trials without refitting:

```matlab
calibration = runMaxPoolingCalibration( ...
    'MaxPoolingTrialCache.mat',result,'FitOutput');
```

This saves group, reliability-bin, session, and condition tables. Reliability
bins contain approximately equal trial counts within each domain-animal group.
Condition summaries use probe offset and step magnitude when those cached trial
fields are available. Calibration intercept and slope have ideal values zero
and one, respectively. Accuracy residuals are observed minus predicted.

For the performance-matched diagnostic model, retain the completed primary fit
and calibrate one noise SD for each domain-animal-offset-step-magnitude cell:

```matlab
performanceMatched = runPerformanceMatchedMaxPoolingFit( ...
    'MaxPoolingTrialCache.mat',result,'FitOutput');
```

This does not replace the three-domain primary likelihood model. It minimizes
the difference between observed and predicted accuracy within each experimental
condition so subsequent gain-recovery analyses can compare choices at matched
behavioral set points. The saved condition table reports the original and
performance-matched residuals, fitted noise SD, attainable accuracy range,
boundary status, and likelihood.

After performance matching, diagnose near-opponent competition with:

```matlab
diagnostics = runMaxPoolingMechanismDiagnostics( ...
    'MaxPoolingTrialCache.mat',result,performanceMatched,'FitOutput');
```

The default comparison includes all well-sampled IDS offsets at or beyond 160
degrees, all well-sampled IDR offsets at or beyond 170 degrees, and the nearest
45-degree step-magnitude control for each domain-animal combination. It reports
drift-region versus probe-region winners within each RDK, changed-minus-
unchanged RDK margins, winner directions relative to drift, and attainable
accuracy as noise varies. Diagnostic winner statistics use a fixed random sample
of at most 1,000 trials per condition and 64 internal-noise draws by default.

Compare observed gains with strict primary-model and performance-matched gains:

```matlab
gainComparison = runMaxPoolingGainComparison( ...
    'MaxPoolingTrialCache.mat',result,performanceMatched,'FitOutput');
```

The same binomial choice-noise regression is applied to observed choices and
both sets of model probabilities. Separate changed- and unchanged-RDK drift
and probe coefficients are estimated at every offset with session fixed
effects. Coefficients are divided by the observed no-step-noise psychometric
slope to express gains in equivalent step-coherence units. If the ProbeSession
cache omits no-step-noise trials, the reference slope is instead estimated from
the step-noise trials with probe-offset fixed effects; the reference table
records the method used for every domain-animal group. IDR uses its one
yoked probe-noise trace; IDQ/IDS use the sum of the two independent probe
traces, so the constrained coefficient is the gain per probe stream. The run
saves a long-form CSV, reference-psychometric table, MAT result, and separate
changed- and unchanged-RDK comparison PDFs.

The same run also saves two compact progress plots overlaying IDR, IDQ, and IDS
changed-side gain as a function of probe offset. Filled symbols are observed
values, open symbols are fitted values, and marker shape identifies domain.
Separate PDFs show the primary domain-noise prediction and the
condition-performance-matched prediction. Each includes drift gain, probe gain,
and probe/drift gain by animal:

- `MaxPoolingGainComparison_ChangedAcrossDomains_Primary.pdf`
- `MaxPoolingGainComparison_ChangedAcrossDomains_PerformanceMatched.pdf`

Observed gains supported by fewer than three sessions are marked as
preliminary and plotted without a confidence bar. A one-session
cluster-robust interval is not identifiable; its exported gain SE and CI are
therefore `NaN` rather than a spuriously narrow interval.

To test whether prestep adaptation increases disproportionately with the
directional prestep response, run the nested exponent profile:

```matlab
exponentProfile = runMaxPoolingAdaptationExponentProfile( ...
    cachePath,result,outputFolder, ...
    'ExponentGrid',0:0.1:1,'NDraws',64);
```

The tested response is

```text
Rstep(m) - a*Ppre(m)*(Ppre(m)/P0)^q
```

with the existing `a=0.56` held fixed, `P0=40%` coherence, and one noise SD
reoptimized per domain. `q=0` is exactly the current linear adaptation model.
If the coarse profile favors `q>0`, refine the exponent with 256 draws before
jointly reassessing `a` and `q`.

The power-subtraction profile changes a direction-dependent offset but does
not reduce local step-response gain. Test genuine prestep-dependent response
compression with:

```matlab
gainProfile = runMaxPoolingAdaptationGainProfile( ...
    cachePath,result,outputFolder, ...
    'GainGrid',0:0.1:1,'NDraws',64);
```

This profiles `b` in

```text
(Ppre(m)-a*Apre(m)) + exp(-b*Ppre(m)/P0)*(Rstep(m)-Ppre(m))
```

while holding `a=0.56` and `P0=40%` fixed and reoptimizing one noise SD per
domain. `b=0` exactly recovers the current model. Positive `b` preserves the
adapted prestep representation while reducing the local gain of new step and
external-noise evidence in mechanisms activated strongly during the prestep.

For the final joint check, sparsely profile subtractive adaptation `a` and
increment-gain strength `b`, reoptimizing domain noise at every grid cell:

```matlab
gainProfile2D = runMaxPoolingAdaptationGain2DProfile( ...
    cachePath,result,outputFolder,'NDraws',64);
```

The default grid is `a = [0 0.3 0.56 0.8 1]` by
`b = [0 0.1 0.25 0.5 1]`. The existing `a` and the nested `b=0` model are
always inserted if omitted. MAT, long-form CSV, and likelihood-surface PDF
outputs use the stem `MaxPoolingIncrementAdaptationGain2DProfile`.

## Incremental experimental cache

The standard cache route is incremental:

```matlab
projectRoot = '/Users/Shared/Data/IDReadout';
cachePath = fullfile(projectRoot,'Data','MaxPoolingTrialCache.mat');

% Preview only: no MAT file is opened and nothing is written.
[~,report,manifest] = syncExperimentalMaxPoolingCache(projectRoot,cachePath, ...
    'DryRun',true);
report.actions       % one row per source: extract, replace, missing, or reuse

% Extract new/modified files, recompute kernels from cached framewise values,
% and atomically update the cache.
[cache,report] = syncExperimentalMaxPoolingCache(projectRoot,cachePath);
```

The full-prestep adaptation model introduces cache schema v2. When upgrading
an existing v1 cache, run one explicit restartable rebuild:

```matlab
[cache,report] = syncExperimentalMaxPoolingCache( ...
    projectRoot,cachePath,'Rebuild',true);
```

The previous final cache remains untouched until all v2 fragments have been
extracted, consolidated, validated, and atomically saved.

For a routine update and refit in one command:

```matlab
[result,report] = updateAndRunMaxPoolingAnalysis( ...
    projectRoot,cachePath,'FitOutput');
```

The saved fit result includes the source manifest and cache update time, so
successive monitoring analyses retain the exact session inventory they used.

Each scan compares the absolute source path, byte count, modification time,
and extractor version. Unchanged ProbeSessions are never opened. A new file is
appended; a modified file replaces only its own cached trials. If the adapter
version changes, every affected file is deliberately refreshed.

Every completed ProbeSession extraction is immediately saved as a separate
atomic restart checkpoint in `<cachePath>.staging`. If extraction later stops
because of an error or interruption, rerunning the same command reports those
sources as `resume` and loads their checkpoints instead of reopening the raw
ProbeSessions. A checkpoint is accepted only when source path, file size,
modification time, and extractor version all match the current scan.

After all sources have been processed, the checkpoints are consolidated,
leave-one-session-out kernels are recomputed, and the final cache is saved
atomically. Checkpoints are deleted only after that final save succeeds. To
discard an interrupted build deliberately, use:

```matlab
'ClearCheckpoints',true
```

ProbeSession discovery goes through the project's `selectAnalysisFiles.m` in
its filename-only `LoadMetadata=false` mode. Consequently the existing hard
filename and path exclusions remain the single authoritative list, but routine
scans do not load headers from unchanged MAT files. The updated selector must
be installed at the project root or otherwise available on the MATLAB path.

Missing files stop a real update by default. After inspecting the dry-run
report, choose explicitly if needed:

```matlab
'MissingFilePolicy','retain'  % keep cached trials
'MissingFilePolicy','remove'  % remove their cached trials
```

After first enabling centralized exclusions, a previously cached hard-excluded
session will appear as `missing`. Review `report.actions`, then use `remove`
once to delete those old cached rows deliberately.

The cache is written to a temporary MAT file and moved into place only after
all extraction, kernel recomputation, and validation succeed. Thus a failed
new-file extraction leaves the previous final cache untouched while retaining
the successfully extracted per-session restart checkpoints.

The quick fingerprint intentionally avoids reading or hashing large MAT files.
It will not detect a content change that preserves both byte count and file
modification time; normal project file creation and replacement do not preserve
both, so this is an appropriate speed/safety tradeoff. `Rebuild=true` forces a
complete refresh when desired.

For a complete synthetic check:

```matlab
exampleSyntheticFit
```

Before building the experimental cache, audit one ProbeSession per domain:

```matlab
report = validateRepresentativeProbeSessions([idrFile idqFile idsFile]);
```

This validates field conventions, trial selection, increment/decrement coding,
and the assumption that stored coherence traces are zero-centered noise.

Run the lightweight tests with:

```matlab
results = runtests('tests');
```

## Cache contract

The cache is deliberately independent of the IDR/IDQ/IDS raw-file layouts.
It contains one row per valid completed trial, including no-step-noise trials.
The primary fit selects `valid & hasStepNoise & isIncrement`.

Required fields:

- `cache.schemaVersion`: `"MaxPoolingTrialCache-v2"`.
- `cache.trials`: table with `animal`, `domain`, `sessionID`, `trialID`,
  `choiceRDK`, `changeRDK`, `correct`, `valid`, `hasStepNoise`, and
  `kernelFamily`, plus `stepSignIndex` (1=decrement, 2=increment) and
  `isIncrement`. Experimental gain comparison additionally requires
  `baseCohPC`, `stepEndpointCohPC`, `stepMagnitudePC`, and `probeOffsetDeg`.
- `cache.streamDirDeg`: `nTrials x 2` cell array. Each cell is a vector of
  stream directions for one RDK on one trial.
- `cache.stepCohPC`: `nTrials x 2` cell array. Each cell contains the
  corresponding kernel-weighted *total* step-period coherences in percentage
  points, including probe pedestals.
- `cache.preCohPC`: `nTrials x 2` cell array aligned with `streamDirDeg`.
  Each cell contains the mean total prestep coherence of every stream,
  including drift and probe pedestals. These vectors generate `Ppre`.
- `cache.driftDirDeg`: `nTrials x 2` numeric array.
- `cache.preDriftCohPC`: retained as an audit field containing the nominal
  prestep drift coherence; the model adaptation term uses `preCohPC` instead.

Optional framewise fields used to construct leave-one-session-out kernels are
described in `maxpoolfit.crossfitKernels`.

`maxpoolfit.buildCache` accepts a project-specific reader function and writes
the combined versioned cache. The reader is the only component that should
know the three raw formats. This keeps stimulus extraction auditable and keeps
the fitted model identical across domains.

## Drift-randomization design groups

Before interpreting domain-specific parameters, profile four mutually
exclusive design groups:

```matlab
designProfiles = runMaxPoolingDesignGroupProfiles( ...
    cachePath,result,outputFolder, ...
    'AdaptationGrid',0:0.05:1,'NDraws',64);
```

The groups are IDR, IDS, IDQ-randomized, and IDQ-fixed. For IDQ,
`sessionHeader.randomBaseDir` is authoritative when available; the number of
cached `driftDirectionDeg` values outside the fixed 30/150/270-degree triplet
supplies an audit and a fallback for historical files. Merely observing more
than one drift direction does not imply random base orientation, because the
fixed-base condition selects among those three directions. The original
`domain` field is not changed.

Only groups represented among the primary trials are fitted. If either IDQ
regime is absent, its independent profile is skipped and the within-IDQ model
comparison is retained with `estimable=false`, `NaN` likelihoods, and an
explanatory note; the available across-design comparisons continue.

This diagnostic does not require rebuilding the trial cache. It reads only
the small IDQ session headers and joins their session-level design metadata to
the existing cached trials. Outputs include the session-design audit, four
independent adaptation/noise profiles, and nested likelihood summaries for:

1. pooled versus separate fixed/randomized IDQ fits; and
2. pooled, shared-adaptation/separate-noise, and fully separate fits for
   randomized IDQ versus IDS.

The comparisons are in-sample diagnostics. A useful split should subsequently
be confirmed by held-out-session prediction, calibration, and gain residuals.

The selected grouped model uses the profile-supported compromise: separate
adaptation for IDQ-fixed and IDR, one adaptation shared by IDQ-randomized and
IDS, and four independently fitted internal-noise SDs. The adaptation values
are selected directly from the saved 64-draw profiles; they are not hardcoded
by the runner. Refit at 256 draws and rerun calibration and gain diagnostics:

```matlab
groupedOutputFolder = fullfile(outputFolder,'DesignGroupFit');

groupedResult = runMaxPoolingDesignGroupFit( ...
    cachePath,result,designProfiles,groupedOutputFolder);

groupedCalibration = runMaxPoolingCalibration( ...
    cachePath,groupedResult,groupedOutputFolder);

groupedPerformanceMatched = runPerformanceMatchedMaxPoolingFit( ...
    cachePath,groupedResult,groupedOutputFolder);

groupedGainComparison = runMaxPoolingGainComparison( ...
    cachePath,groupedResult,groupedPerformanceMatched,groupedOutputFolder);
```

For these outputs, `domain` in diagnostic tables and figures is deliberately
relabelled as `IDQ-fixed`, `IDQ-randomized`, `IDS`, or `IDR`. The original
cached domain remains unchanged.

Profile the angular correlation length of the internal Gaussian noise field
without changing MT direction tuning:

```matlab
correlationProfile = runMaxPoolingCorrelationLengthProfile( ...
    cachePath,groupedResult,groupedOutputFolder, ...
    'CorrelationLengthGridDeg',[5 10 15 20 30 45 60 90], ...
    'NDraws',64,'MinimumGainSessions',3);
```

At each correlation length, the three selected adaptation values remain fixed
and all four design-group noise SDs are reoptimized. The profile reports choice
likelihood and changed-side gain error separately. Gain RMSE is unweighted
across retained domain-animal-offset cells; cells with fewer than three sessions
are saved descriptively but excluded from the objective. The initial sweep is
diagnostic. A narrower interval should be rerun with 256 draws only if likelihood
and gain geometry identify a promising common region.

## Numerical strategy

The likelihood is reproducible but not mathematically smooth: fixed draws plus
a hard maximum give a finely stepped simulated likelihood. `fitModel` optimizes
one noise SD independently for each domain, using bounded scalar searches on
log noise SD. The common adaptation value was selected by profile likelihood
and is held fixed at 0.56 while the domain-noise values are optimized.
Development and final draw counts are configured separately.

The two RDK maximum distributions are sampled independently. Their comparison
uses every cross-pair of sampled maxima, rather than only paired draws, to
reduce numerical variance.

The hard-max rule is an increment model: its primary selector is
`valid & hasStepNoise & isIncrement`. Decrement trials remain in the cache for
provenance and later model extensions but do not enter this fit.

## Integration note

Before rebuilding the experimental cache, the updated adapter should be checked
against representative IDR, IDQ, and IDS ProbeSessions. The validation report
now displays the mean prestep drift and probe coherences, making pedestal
handling directly auditable. The schema-v2 rebuild is restartable; subsequent
routine runs use the incremental manifest route above.
## Independent and response-scaled mechanism noise

After completing a design-group fit, first test whether adding independent
homoscedastic mechanism variance changes the choice likelihood and changed-side
gain geometry:

```matlab
noiseProfile = runMaxPoolingIndependentNoiseProfile( ...
    cachePath,result,outputFolder, ...
    'IndependentFractionGrid',[0 .05 .1 .2 .4 .7 1], ...
    'VarianceExponentGrid',0,'NDraws',64);
```

Here `lambda` is the fraction of pointwise variance assigned to independent
mechanism noise. The remaining variance has the configured circular directional
correlation. `lambda = 0` exactly recovers the prior model. After selecting one
or two useful lambda values, test response-dependent variance with, for example,
`'VarianceExponentGrid',[0 .5 1 1.5 2]`. The exponent specifies variance
scaling, so `p = 1` is the Poisson-like case. Overall, group, gain, and
within-RDK winner-fraction CSV files are written to the output folder.

## Gain-targeted directional-width fit

To select effective directional width from the changed-side gain geometry in
IDQ-randomized and IDS while retaining trial-level choice performance as a
guardrail, run:

```matlab
widthFit = runMaxPoolingDirectionalWidthGainFit( ...
    cachePath,result120,outputFolder, ...
    'DirectionalWidthGridDeg',[10 15 20 25 30 37.5 45 60], ...
    'TargetGroups',["IDQ-randomized" "IDS"], ...
    'MinimumGainSessions',3,'NDraws',64);
```

At each width, one noise SD per design group is refitted by choice likelihood.
The primary mechanistic objective then gives equal weight to the observed-versus-
predicted changed drift and probe gains in every retained animal-design-offset
cell. Cells with fewer than three sessions are plotted but excluded from the
gain objective. The output reports gain-optimal and choice-optimal widths
separately and saves observed-versus-fitted gain curves.
### Joint gain-targeted mechanism fit

After a completed design-group fit, jointly optimize effective directional
width, prestep adaptation, and directional-noise correlation length for the
changed-side gains in randomized IDQ and IDS:

```matlab
jointGainFit = runMaxPoolingJointGainFit( ...
    cachePath,result,outputFolder);
```

The runner accepts the grouped result, its `finalFit`, or the selected fit
from the directional-width analysis. If the older primary `result` remains in
the workspace, it automatically uses `MaxPoolingDesignGroupFitResults.mat`
from `outputFolder` when that file is present.

At every candidate point, one upstream noise SD per target design group is
profiled by trial-level choice likelihood. The search uses fixed Monte Carlo
draws, space-filling initialization, and bounded multistart `fminsearch`.
The output retains all evaluations and identifies both the gain optimum and
the choice-likelihood optimum, together with their Pareto frontier.

Decompose the fitted opponent-axis effect into deterministic response gaps,
within-RDK winners, probe-noise-induced winner switches, and final gain:

```matlab
opponentAxis = runMaxPoolingOpponentAxisDecomposition( ...
    cachePath,jointGainFit,outputFolder);
```

The diagnostic compares the actual paired probes with two single-180-degree
controls that retain each actual independent probe trace in turn.
