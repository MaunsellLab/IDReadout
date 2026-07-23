# IDQ–IDR deterministic MT candidate diagnostic

This stage evaluates the actual step-rectangular coherence distributions
without internal variability or simulated choices. It is deliberately a
forward-model audit rather than a behavioral fit.

## Model contract

- One signed, mean-subtracted MT population response is computed for each
  patch and trial.
- Signed IDQ noise predictors are deviations about a positive coherence
  pedestal. The specified pedestal is placed in all three physical directions
  on both patches, and the step is added to the changed-patch drift pedestal.
- Gaussian readouts are equally normalized to unit peak.
- Candidate directions include every physical component direction and its
  opponent. Coincident candidates are collapsed and retain all role labels.
- Opponent signs and relative magnitudes arise from MT direction tuning; no
  fixed preferred:null divisor or opponent gain is imposed.
- Signed and candidate-rectified activations are evaluated separately.
- A hard maximum is taken within each patch. No internal noise, patch-choice
  rule, or stochastic response is added at this stage.
- IDR two-probe yoking is preserved by applying the per-candidate probe
  coherence at both `+theta` and `-theta`. The 180-degree condition has one
  physical probe stream. The sole 1-degree condition remains excluded.

## Isolated run

Create a new synthetic context and declare both experimental summaries as
manifest sources. Use `SaveLatents=true` if the next choice-generation stage
will consume the candidate arrays.

```matlab
idqSource = fullfile(domainFolder(mfilename('fullpath'),'IDQ'), ...
  'Data','AcrossSessionSummaries','IDQ_AcrossSideSummary.mat');
idrSource = fullfile(domainFolder(mfilename('fullpath'),'IDR'), ...
  'Data','AcrossSessionSummaries','IDR_SideGainAnalysis.mat');

context = makeAnalysisContext("synthetic", ...
  'RunID',"IDRIDQ_MTTrialCandidates_001", ...
  'Seed',1729,'SaveLatents',true);
initializeSyntheticRun(context, ...
  'ModelName',"IDRIDQDeterministicMTTrialCandidates", ...
  'ModelParameters',struct('sigmaMTDeg',37.5, ...
    'sigmaReadoutDeg',[15 20 25], ...
    'idqPedestalPC',7), ...
  'SourceFiles',[string(idqSource);string(idrSource)]);

summary = analyzeIDRIDQMTTrialCandidates(context, ...
  'SigmaMTDeg',37.5,'SigmaReadoutDeg',[15 20 25], ...
  'IDQPedestalPC',7, ...
  'LatentSigmaReadoutDeg',20);
```

For the existing IDQ data, `sessionHeader.noiseAmplitudePC` is 7% in every
session. This is both the binary-noise half amplitude and the positive
pedestal: each physical stream therefore alternates between 0% and 14%
coherence, and the three streams can consume at most 42% of the dots.
`IDQPedestalPC` remains mandatory because the across-session summary retains
the signed deviations but not this generating header field. A vector indexed
by `trialTable.sessionIndex` is also accepted for future datasets in which the
amplitude differs between sessions.

The analysis refuses undeclared or subsequently changed source files. Summary
MAT/CSV, exact candidate-role winner frequencies, signed-noise association
tables, and the PDF are written only inside this run. The signed-noise table is
descriptive conditioning on the observed multistream trials, not a causal
single-stream perturbation. Trialwise latent arrays are written only when the
context was created with `SaveLatents=true`.

## Tests

```matlab
results = runtests(["testMTGaussianForwardCore.m", ...
  "testIDRIDQMTTrialCandidates.m"]);
assertSuccess(results)
```
