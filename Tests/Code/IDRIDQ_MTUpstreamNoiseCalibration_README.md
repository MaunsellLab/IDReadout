# Upstream directional-noise calibration

This is the primary internal-noise model. It replaces independent additive
MT-response noise; it is not layered on top of it.

For each patch and trial, a signed Gaussian perturbation is drawn over a dense
bank of latent motion input directions. The existing MT tuning transformation
acts on physical coherence and internal directional noise identically. Gaussian
readouts then pool the resulting MT responses. No separate MT-response or
readout noise is added.

IDQ temporal contract:

- Every trial has three 7% directional pedestals throughout prestep.
- On `hasStepNoise` trials, the pedestals and signed noise continue during the
  step and the physical step is added to the changed-patch drift stream.
- On `~hasStepNoise` psychometric trials, the directional noise and pedestals
  terminate at step onset. Thus the decision-epoch no-change baseline is zero.
- The present rectangular calibration uses the step decision epoch. The
  prestep distribution is retained for later temporal-weighting models; actual
  prestep noise fluctuations are not reconstructed.

Run tests:

```matlab
results=runtests("testMTUpstreamNoiseCore.m");
assertSuccess(results)
```

Initialize a fresh manifested run and calibrate:

```matlab
contextUpstream=makeAnalysisContext("synthetic", ...
  'RunID',"IDRIDQ_MTUpstreamNoise_001", ...
  'Seed',1729,'SaveLatents',false);
initializeSyntheticRun(contextUpstream, ...
  'ModelName',"DenseUpstreamDirectionalNoiseCalibration", ...
  'ModelParameters',struct('sigmaMTDeg',37.5, ...
    'sigmaReadoutDeg',20,'candidateSpacingDeg',2, ...
    'inputSpacingDeg',1,'idqPreStepPedestalPC',7, ...
    'nMonteCarlo',2000), ...
  'SourceFiles',[idqSource;idrSource]);

calibration=calibrateIDRIDQMTUpstreamNoise(contextUpstream, ...
  'SigmaMTDeg',37.5,'SigmaReadoutDeg',20, ...
  'CandidateSpacingDeg',2,'InputSpacingDeg',1, ...
  'IDQPreStepPedestalPC',7,'NMonteCarlo',2000);
```
