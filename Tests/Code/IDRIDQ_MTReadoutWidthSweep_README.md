# Gaussian readout-width sweep

Package revision 1.1 fixes MATLAB dynamic table-variable construction in the
summary plots. Completed member simulations and gain recoveries from revision
1.0 remain valid and are resumed without recomputation.

This controller runs the same threshold calibration, actual-trial replicated
choice simulation, and matched gain recovery at Gaussian readout widths of
1, 5, 10, and 20 degrees. Five choice replicates are generated per width by
default.

Each width is a separate manifest-backed synthetic run. A rerun loads complete
stages and resumes gain-recovery checkpoints. Once all four runs finish, a
fifth summary run declares the member choice and recovery MAT files as hashed
sources and writes combined tables and a three-page PDF. The summary also
reports the approximate effective tuning width and the forward-model
preferred:null sensitivity ratio at every readout width.

No experimental or derived arrays are copied between run folders. The source
summary files are read through the existing manifest validation.

Install all package files together in `Tests/Code`, then run:

```matlab
clear runIDRIDQMTReadoutWidthSweep formatReadoutWidthRunID
clear recoverIDRIDQMTUpstreamNoiseGains matchSyntheticTrialRows analyzeIDRMatchedGains
rehash
results=runtests("testMTReadoutWidthSweepCore.m");
assertSuccess(results)
```

The package includes the corrected gain-recovery wrapper and matched IDR
fitter. It expects the previously installed upstream calibration/choice core
and the repository's `fitIDQNoiseGain.m`. A useful preflight is:

```matlab
which calibrateIDRIDQMTUpstreamNoise -all
which simulateIDRIDQUpstreamNoiseChoices -all
which recoverIDRIDQMTUpstreamNoiseGains -all
which fitIDQNoiseGain -all
```

Start the sweep with a new, unique prefix:

```matlab
sweep=runIDRIDQMTReadoutWidthSweep( ...
  'RunIDPrefix',"MTReadoutWidth_20260722_001", ...
  'ReadoutWidthsDeg',[1 5 10 20], ...
  'NReplicates',5);
```

If MATLAB stops, issue the identical call again. Existing manifests are
verified; completed calibration and choice stages are loaded; gain recovery
continues from the last committed replicate.

The 200 MB quota applies independently to each member run. The summary run has
a 100 MB quota. Full MT and candidate latents are not saved.

The same seed is intentionally used for every width. This reduces avoidable
Monte Carlo differences between conditions, although different covariance
ranks mean the member simulations are not perfectly paired random draws.
