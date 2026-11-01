# Run locally and share compact results

The recordings stay on your hardware. There is no need to upload the roughly
0.5-GB ProbeSessions or the full trial cache.

## First run

Put this package on the MATLAB path alongside the existing IDS analysis code
and the current **selectAnalysisFiles.m**. Keep using that project function;
this package does not install a replacement selector.

```matlab
runMTReadoutPNormLocal
```

Choose the **IDS** folder or **IDS/Data** once in the folder chooser. The
launcher finds its `Probe*/ProbeSessions` directories, calls your selector,
and remembers the root for subsequent runs. Selecting an individual
ProbeSessions folder also works. A directly selected converted-session
folder works if its headers carry the selector's SourceDomain metadata.

Selection explicitly uses Animal=Meetz, SourceDomain=IDS, IncludeIDQ=false,
LoadMetadata=true, and Bin180Into179=false. Your current hard exclusions
remain in effect. The measured 179/180 offsets are retained separately.
Missing metadata errors come from the existing selector; selection never
falls back to guessing the task from filenames.

Outputs go into `MTReadoutPNormOutput` under the selected root. After a fit,
send **MTReadoutPNormShare.zip** from that folder. It contains aggregate
coverage, fit parameters, likelihood/performance by offset, and recovered
gains/plots when available. It also records any missing-summary or gain
recovery error. It excludes raw recordings, frame traces, individual-trial
inputs and predictions, and the local input cache.

## A short check before the full run

```matlab
validateMTReadoutPNormLocal
runMTReadoutPNormLocal('MaxFiles',2)
```

The validation uses small synthetic fixtures in a temporary folder and your
real selector. The two-file check verifies local I/O and fitting. Its compact
results explicitly identify it as a subset run; gain recovery using a full
observed summary is skipped. File selection is deterministic, so the first
two selected files may not cover both 45 and 120 degrees.

An extraction/coverage-only check is also available:

```matlab
runMTReadoutPNormLocal('ExtractOnly',true)
```

For the full fit, run the unqualified command again. This is optional; the
launcher does not require a staged approval before the full run.

## Large inputs and repeated runs

The first run extracts the two 250-ms windows one source file at a time.
It loads only the relevant top-level variables, although loading a
ProbeSession's probeTrials cell can still require substantial memory for
that one source file. Compact per-file fragments are saved immediately.

Unchanged fragments are reused on subsequent runs. Source path, byte count,
modification timestamp, extraction version, animal, and window settings
guard reuse. Changes to p, tuning width, noise, or integration precision do
not require re-reading the raw files. The selector still reads the small
headers to apply your current exclusions. To deliberately re-extract:

```matlab
runMTReadoutPNormLocal('RebuildCache',true)
```

The fitter constructs directional population responses in trial batches
(128 by default), rather than allocating a whole-dataset response tensor.
Runtime still scales with trials, directions, draws, and candidate p values.
Each likelihood evaluation reports elapsed seconds. A fit checkpoint is
saved after each completed p; interrupted fits are not automatically
resumed, but input extraction is reusable.

To change the root, use a folder chooser again:

```matlab
runMTReadoutPNormLocal('ChooseFolder',true)
```

## Gain recovery

The launcher looks for existing `acrossSummary` structs under each selected
offset's `AcrossSessionSummaries` directory. It uses only saved summaries
explicitly marked as Meetz and containing native IDS sessions. The gain
adapter uses their exact existing predictors and observed psychometric
calibration. It never substitutes a different regression.

If suitable summaries have not been saved, create them through the existing
IDS pipeline:

```matlab
S = plotIDSSummary('Animal','Meetz','IncludeIDQ',false, ...
    'NoisePredictorType','looStep');
finishMTReadoutPNormGains('Summaries',S)
```

Choose the completed result MAT file in the chooser. This adds gains and
regenerates the compact ZIP **without repeating the forward fit or reading
the large source recordings**. With saved matching summaries, simply call
`finishMTReadoutPNormGains` with no arguments.

A saved fit remains available if gain recovery fails. The share ZIP records
the failure instead of labeling the gain comparison complete. All summary
trials must have matching model predictions. Selected source files and the
observed summaries must therefore use the same exclusions and trial set.

## Model and validation scope

The agreed model is unchanged: full directional responses; 37.5-degree
tuning; 1-degree bank; equal 250-ms pre/post windows; independent internal
noise fields in each window and RDK; the same p everywhere within a
condition; and signed differences of the two pooled responses. See README
for the physical-rate floor and the interpretation of the isotropic dot
response. No new temporal-correlation or task-specific adaptation parameter
was introduced for the local workflow.

MATLAB/Octave is not available in the development environment. The earlier
Python mathematical checks and two-session pilot were executed; the new
MATLAB local-workflow integration checks are supplied for your MATLAB
environment and have not been represented as executed here.
