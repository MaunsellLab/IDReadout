# IDJ session analysis and plotting

Place the fourteen files from `Code` in your existing `IDJ/Code` folder. This
includes an updated `updateIDJ.m`; retain your existing `convertIDJData.m` and
`makeIDJSessionHeader.m`. The shared lab utilities `domainFolder`, `validFolder`,
`selectAnalysisFiles`, and the converter's existing dependencies must remain
on the MATLAB path. There are no dependencies on IDL- or IDQ-specific functions.
The MATLAB features used are the same modern table, string, tiledlayout, and
exportgraphics features used by the supplied IDL pipeline.

## Run

```matlab
updateIDJ('Animal', 'Meetz');
updateIDJ('Animal', 'Neesha');
```

For already converted sessions, skip conversion:

```matlab
makeIDJSessionAnalyses('Animal', 'All');
makeIDJSessionSummaryPlots('Animal', 'All');
```

Rebuild analyses and plots:

```matlab
makeIDJSessionSummaryPlots('Animal', 'All', ...
    'ReplaceAnalysis', true, 'Replace', true);
```

`updateIDJ('ReplaceAnalysis', true)` also rebuilds analyses and plots;
`updateIDJ('Replace', true)` rebuilds only plots. Existing conversion behavior is
unchanged: conversion checks all `.dat` files, and Animal filters later stages.

Inputs remain `IDJ/Data/FullSessions/*.mat`, including the supplied files whose
names start with `IDR6`. The compact outputs are
`IDJ/Data/SessionAnalysis/*_sessionAnalysis.mat`; PDFs are
`IDJ/Plots/FullSessions/*_IDJSummary.pdf`.

## Files and behavior

| File | Purpose |
| --- | --- |
| `makeIDJSessionAnalysis.m` | Extract one compact session |
| `makeIDJSessionAnalyses.m` | Build/cache compact session files |
| `makeIDJSessionSummaryPlots.m` | Build/cache session PDFs |
| `plotIDJSessionSummary.m` | Four allocation panels, threshold comparison, and session text |
| `fitIDJSessionPsychometrics.m` | Allocation-specific thresholds, shared Weibull beta and lapse |
| `idjWeibullP.m` | IDJ psychometric probability function |
| `idjPlotFolder.m` | Direct plots to IDJ/Plots |
| `updateIDJ.m` | Run conversion, extraction, and plotting |

The summaries preserve the IDL design: separate zero-lapse and free shared-lapse
fits, zero-lapse curves in the allocation panels, 75% thresholds with Hessian
intervals, and counts/symbols showing trials and configuration epochs. The lapse
bounds remain [0, 0.15]. Allocation labels now have three components. Missing
allocation conditions in a partial session are permitted; fits use observed
conditions. Sessions with no certified correct/wrong trials fail explicitly.

## Allocation and noise contract

- `trialTable` keeps the recorded `alloc0`, `alloc1`, `alloc2`, all three endpoint
  coherences and directions, and `topUpTrial`.
- `unsignedAllocPC` is the nominal largest allocation: 100, 70, 50, or 100/3.
  `minorAllocPC` is the nominal per-minor-stream allocation: 0, 15, 25, or 100/3.
- The recorded 0.33330/0.33335/0.33335 condition is grouped as exact equal thirds.
  Recorded fractions and coherences are preserved; grouping does not alter the
  stimulus records or the noise arrays. `highStreamIndex` is 0 for equal thirds,
  otherwise a MATLAB stream index 1, 2, or 3.
- Noise arrays use `[frame, analysis trial, stream]`. At 120 Hz with a 750-ms
  prestep and 250-ms step, full arrays have 120 frames and step arrays 30 frames.
  Event times are filled forward exactly as in IDL, separately for all streams.
- `streamTable` has three rows per behavioral trial, ordered by trial and then
  stream. It records `analysisTrialRow`, `sourceTrialIndex`, `streamIndex`,
  `recordedAllocation`, nominal `allocationPC`, and `role` (major/minor/equal).
  Pooling minor rows yields 2N samples for N trials at an unequal allocation;
  pooling equal rows yields 3N samples for N equal-thirds trials. Stream samples
  share their trial's outcome; the table does not make them independent trials.
- Behavioral fits use `trialTable` once per trial, never `streamTable`.
- `stepTable.nAllocationConfigurations` counts configured forms; IDL's
  two-signed-forms assertion has been removed. Its high/low endpoint summaries
  follow nominal major/minor roles; equal-thirds endpoints are averaged to avoid
  elevating one stream because of rounding. Raw endpoints remain in trialTable.
- Top-up trials are always excluded before trial, configuration, stream-table,
  or noise extraction. This also excludes them from psychometrics. The retained
  `topUpTrial` column is always false. Both supplied sessions contain no
  otherwise-analysis-eligible top-up trials.

To pool change-side noise according to streamTable ordering:

```matlab
A = makeIDJSessionAnalysis(sourcePath);
S = A.streamTable;
N = reshape(permute(A.changeNoiseByFrameTrial, [1 3 2]), ...
    numel(A.stepTimeMS), []);
minorNoise = N(:, S.role == "minor");
equalThirdsNoise = N(:, S.role == "equal");
```

## Verification

The supplied Meetz and Neesha sessions contain 1,336 and 1,095 certified
correct/wrong trials, respectively. Independent Python checks confirmed the
recorded allocations, 120-degree separation, step/block agreement, frame-noise
coverage for every trial and all six patch/stream channels, expected 2N/3N
pooling counts, and finite fits under the same Weibull parameterization.
These checks validate the input data and numerical conventions; they do not
execute the MATLAB implementation.

MATLAB is unavailable in the build environment, so MATLAB execution and PDF
rendering remain unverified. A MATLAB integration check is included:

```matlab
addpath('/your/path/to/IDJ_pipeline/Tests');
validateIDJSessionPipeline('/your/path/to/IDR6_Meetz_20261020.mat');
validateIDJSessionPipeline('/your/path/to/IDR6_Neesha_20261020.mat');
```

It runs extraction, checks source noise ordering and sample counts, fits both
models, and creates/closes a hidden summary figure without exporting a PDF.

## Across-session kernels, gains, and model predictions

This update adds:

- `plotIDJAcrossSessionBehavior.m`: IDJ adaptation of the supplied IDL function.
- `fitIDJModelRecoveryAcrossSessions.m`: three-stream adaptation of the IDL
  model-recovery dependency, preserving the temporal-weight estimator.
- `idjAllocationLabel.m` and `idjWeibullAlphaForThreshold.m`: shared helpers.
- `Tests/validateIDJAcrossSessionPipeline.m`: integration check without saving.

Install all files from `Code`, including the revised `updateIDJ.m`. The new
`AcrossSessions` option defaults to true, so the complete pipeline now runs as:

```matlab
updateIDJ('Animal', 'All');
```

For existing compact session analyses, run only the new stages:

```matlab
result = fitIDJModelRecoveryAcrossSessions('Animal', 'All', 'MakePlots', false);
[summary, figures] = plotIDJAcrossSessionBehavior( ...
    'Animal', 'All', 'ModelResult', result);
```

The fitter saves `IDJ/Data/AcrossSessions/ModelRecovery.mat`. The plotter can
load that file automatically when `ModelResult` is omitted:

```matlab
[summary, figures] = plotIDJAcrossSessionBehavior('Animal', 'All');
```

Run the fitter again whenever the selected sessions or their data change.
Run both stages with matching selection options. No IDL ModelRecovery file
can be substituted; the IDJ plotter checks the three-stream model contract.
For session-only updates, use `updateIDJ('AcrossSessions', false)`.

Two figures are made for each animal, plus a Combined group when both animals
are selected. PDFs are written to `IDJ/Plots/AcrossSessions`:

- `IDJFunctions_<group>.pdf`: raw and session-aligned psychometrics, allocation
  thresholds, and within-session threshold ratios.
- `IDJKernels_<group>.pdf`: seven fraction kernels and the relative-gain panel.

As in IDL, alignment divides steps by the geometric mean of the four session
75% thresholds. The allocation order is 100/0/0, 70/15/15, 50/25/25, equal thirds.
The across-session stage currently requires all four allocations in each
session, matching the supplied IDL implementation.

The kernel fractions are 0, 15, 25, 100/3, 50, 70, and 100 percent. Each kernel
shows correct-minus-error mean coherence noise across the entire prestep and
step, with session-by-step noise means removed. The major stream contributes
one sample per trial; minor streams are pooled as two samples per trial; equal
thirds pool all three streams. Panel `n` is the number of stream samples **per
patch**. `kernel.nTrials` retains the behavioral count and
`kernel.nStreamSamples` records these 1N/2N/3N sample counts.

For the scalar gain regression, each trial contributes one outcome and one
integrated predictor per patch/stream class. The minor predictor is the SUM of
the two minor noises. At equal thirds the predictor is the SUM of all three.
This constrains exchangeable streams to share a **per-stream** gain, retaining
the same units as the major-stream gain. Using an averaged predictor would
instead multiply the reported coefficient by two or three. Descriptive kernels
are averaged to retain single-stream coherence units; gain coefficients receive
no additional division. Both observed choices and per-trial model probabilities
are projected through the same regression after session-by-step centering.

Temporal weights use the existing IDL procedure: noisy 100/0/0 Change-major
trials from OTHER sessions of the same animal, normalized to unit L1 mass and
shared across all streams and allocations. If no informative training session
exists, the weights are uniform over the 250-ms step. With the supplied sample
files there is one session per animal, so both use this fallback; the Combined
plot does not train one animal's kernel from the other animal.

Model predictions preserve L1, L2, Max, and Square pooling, now over three
streams. L2 is the norm of the three channel responses, rather than the length
of a two-dimensional motion vector. Integration includes recorded noise and
pedestals and precedes pooling. A Gaussian decision sigma is fitted separately
for each session/model using all allocations by default (`SigmaFitMode='all'`).
The alternative `SigmaFitMode='reference'` calibrates sigma at 100/0/0 and scores
divided allocations. The model overlays are conditional expected regression
coefficients from fitted trial probabilities, as in IDL.

Relative gains, including model predictions, divide by that group's observed
100%-Change coefficient. If the reference is not finite and positive, relative
gains are NaN and the figure says they are unavailable. Noise-free sessions
remain excluded from the primary across-session analysis by default; use
`IncludeNoiseFreeSessions=true` consistently in both stages for secondary
behavioral displays. Trials without any recorded noise do not enter the noise
kernel/gain statistics.

The fitter also preserves `KernelMethod='probit'` and `'logit'`,
`ReferenceAllocationPC`, and session-level likelihood summaries. Its new
`MakePlots=false` suppresses its own model-recovery figures while retaining
numerical fitting and saving. `ExportPDF=false` suppresses its PDF export;
`SaveResult=false` suppresses its MAT write. The plotter retains `ExportPDF`,
`Docked`, `Visible`, `ModelResult`, `ModelRecoveryPath`, and selection options.

Verification for this update: all 14 MATLAB files pass MISS_HIT syntax/lint
checks. Independent calculations on both supplied sessions and their combined
group produced finite kernels and gains and verified the sample multiplicities,
frame-integration ordering, and exchangeable-stream normalization. Synthetic
checks verified tied per-stream coefficient scaling and LOO isolation. MATLAB
execution and PDF rendering are still unverified in this environment.

After installation, this runs the new MATLAB integration check without writing
results or PDFs:

```matlab
addpath('/your/path/to/IDJ_pipeline/Tests');
summary = validateIDJAcrossSessionPipeline('Animal', 'All');
```

## P-norm profiles and the sum-of-squares comparison

`pilotIDJPNormProfile.m` ports the supplied `pilotIDLPNormProfile.m` to three
streams. It retains the q=1/p grid, per-session sigma profiling, same-animal
LOO temporal kernels, whole-session bootstrap, lapse options, optional
Session100 scheme, allocation-fit table, and optional parametric recovery.
The minimum is **three distinct eligible sessions of one animal**, even with
`NBootstrap=0`. Unlike the descriptive across-session analysis, this profile
retains the original strict nonzero LOO-kernel requirement; it does not silently
substitute a uniform kernel when training fails.

```matlab
R = pilotIDJPNormProfile('Animal','Meetz','NBootstrap',1000);
R = pilotIDJPNormProfile('Animal','Neesha','NBootstrap',1000);
```

For a quick run without bootstrap:

```matlab
R = pilotIDJPNormProfile('Animal','Meetz','NBootstrap',0);
```

Normal `updateIDJ` behavior is retained. To add profiles for all selected
animals to the update, explicitly opt in:

```matlab
updateIDJ('Animal','All','PNormProfile',true,'NBootstrap',1000);
```

`PNormProfile` defaults to false because the bootstrap can take substantial
time. This option also requires three eligible sessions for every selected
animal; call the profile directly for just an animal with enough sessions.

The upper-left LL panel includes a horizontal dashed line at

```
Square LL - best p-norm LL
```

Its value is fitted using exactly the same selected sessions, per-session
LOO kernels, and lapse settings as the p-norm profile. Square is defined as
`sum(ChangeStream.^2,2) - sum(NoChangeStream.^2,2)` AFTER temporal integration,
retaining all three streams, recorded coherence pedestals, and noise. It has
its own fitted sigma per session, in squared coherence units. No L2 square
root or extra baseline subtraction is applied.

For the primary JointAll scheme, both models fit sigma/lapse on all trials and
score all trials. For optional Session100, both fit sigma/lapse on 100/0/0
trials and score divided allocations. A separate Square line is plotted in
each scheme's LL panel. Its y-axis limits include the line even when Square
outperforms the best p-norm fit. Comparisons are made within a scheme; the
two schemes' absolute LLs have different scoring sets.

The summary table adds:

| Column | Definition |
| --- | --- |
| `bestPNormLL` | Maximum pooled p-norm LL on the supplied q grid |
| `squareLL` | Pooled Square LL under matching fitting/scoring conditions |
| `relativeLLatSquare` | `squareLL - bestPNormLL`, the dashed-line height |
| `deltaLLBestPNormVsSquare` | `bestPNormLL - squareLL`; positive favors p-norm |

Per-session Square sigmas, lapse fractions and LLs are in `R.squareModel`.
`R.bootstrapSquareLL` and `R.bootstrapDeltaLLBestPNormVsSquare` retain the
comparison on every accepted bootstrap draw. Every draw refits the kernels
and both models' session parameters; all copies of the evaluated source
session remain excluded from its kernel training. No extra random draws are
introduced by fitting Square, preserving the p-norm bootstrap/recovery RNG
separation. The p-norm family fits q in addition to its session parameters;
these raw LL differences are descriptive and are not complexity penalties
or significance tests.

Result/PDF names remain animal-specific beneath IDJ:
`Data/AcrossSessions/PNormPilot_<animal>.mat` and
`Plots/AcrossSessions/PNormPilot_<animal>.pdf`. Lapse-run suffixes remain
`_LapseSession` and `_LapseFixedXXpc`. These are IDJ-folder outputs and do not
replace IDL-folder results.

Examples of retained options:

```matlab
R = pilotIDJPNormProfile('Animal','Meetz','IncludeSession100',true);
R = pilotIDJPNormProfile('Animal','Meetz','LapseMode','session');
R = pilotIDJPNormProfile('Animal','Meetz','LapseMode','fixed','LapseRate',.10);
```

The lapse parameter is the fraction of random-guess trials. Thus .10 means
10% guessing and a 95% performance ceiling, as in the source analysis. Entirely
step-noise-free sessions are excluded by default, consistent with the other
IDJ stages; `IncludeNoiseFreeSessions=true` restores them subject to the strict
LOO requirements. Excluded files are recorded in
`R.excludedNoiseFreeSourceFiles`. The global top-up-trial exclusion also applies.

New plotting controls are `MakePlots`, `Docked`, and `Visible`. Use
`'MakePlots',false,'SaveResult',false` for numerical-only runs without writes.
`ExportPDF=false` suppresses export while leaving the figure enabled.

Verification: all 16 delivered MATLAB files pass syntax/lint checks. Independent
numerical checks verified three-stream p-norm pooling (including near-max and
zero rows), the distinction between L2 and Square, matching fit/score masks for
both schemes, fixed-lapse comparisons, q-grid independence of the Square score,
and exclusion of all held-source copies from kernel training. To exercise the
comparison with the supplied single-session fixtures, these checks split each
session into disjoint trial fragments solely for software testing. Those
fragments and their fitted values are not independent experimental sessions
and must not be used for scientific inference. MATLAB execution, the optional
recovery run, and actual MATLAB PDF rendering remain unverified here.

A MATLAB smoke test is included. It writes and removes temporary test-only
trial fragments and runs a small bootstrap to check returned fields and RNG
restoration:

```matlab
validateIDJPNormProfile('/your/path/to/IDR6_Meetz_20261020.mat');
```

## Mandatory top-up exclusion (compact analysis version 2)

Every analysis now excludes trials whose recorded `trial.data.topUpTrial`
flag is true. Exclusion occurs alongside certification/EOT selection, BEFORE
reading step/configuration quantities or extracting noise. This applies to
psychometrics, descriptive kernels, temporal-kernel training, gain regression,
model recovery, p-norm profiles, session bootstrap, and parametric recovery.
There is no option to include flagged trials.

`trialTable`, `streamTable`, all four noise arrays, and configuration/step
counts are constructed from the same retained source-trial indices. Source
indices and trial numbers are preserved. A top-up trial need not match the
regular step grid because it is removed before that validation. The top-up
flag is required and must be scalar 0/1 on otherwise-eligible trials.

Compact outputs now use `sessionAnalysis.version=2`, include
`topUpTrialsExcluded=true`, and report `nExcludedTopUpTrials`: the number of
certified correct/wrong trials removed specifically by the top-up filter.
`nExcludedTrials` still counts all excluded trials, including other EOT and
certification exclusions; the session header's total trial count remains the
original recorded total. Session plots display the excluded top-up count.

The new `idjHasTopUpExclusion.m` helper checks the version, policy marker, and
retained flags. Both batch builders automatically regenerate older compact
caches. Direct downstream entry points reject older/flagged analyses and ask
for a rebuild. The standalone psychometric fitter also rejects input tables
containing flagged trials. Old ModelRecovery results cannot be used for gain
prediction overlays: the result must certify the new top-up exclusion, because
its kernels could have been trained on flagged trials in another session.

After copying all updated Code files, rebuild existing outputs once:

```matlab
updateIDJ('Animal','All','ReplaceAnalysis',true);
```

This regenerates compact analyses, session plots, model recovery, and
across-session functions/kernels/gains. Also rerun p-norm profiles for each
animal with enough sessions:

```matlab
R = pilotIDJPNormProfile('Animal','Meetz','NBootstrap',1000);
R = pilotIDJPNormProfile('Animal','Neesha','NBootstrap',1000);
```

Or, when every selected animal has at least three eligible sessions:

```matlab
updateIDJ('Animal','All','ReplaceAnalysis',true,'PNormProfile',true);
```

No FullSession reconversion is needed: the flag is already present in the
converted trial structures. The supplied sample sessions had no eligible
top-ups, so their ordinary trial counts are unchanged.

Verification for this change: all 18 delivered MATLAB files pass syntax/lint
checks. Independent checks injected top-ups into each supplied session,
including all four allocations and error trials, and deliberately gave those
trials invalid steps. Exclusion occurred before step checks and retained exact
source-index alignment for all six step-noise channels. The MATLAB regression
test additionally checks full-noise alignment, stream tables, counters,
configuration/step counts, stale-cache rejection, and the standalone fitter
safeguard. MATLAB execution remains unverified in this environment.

```matlab
validateIDJTopUpExclusion('/your/path/to/IDR6_Meetz_20261020.mat');
```
