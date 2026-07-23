# Dense-bank MT internal-noise calibration

This stage calibrates internal variability without yet generating synthetic
behavioral datasets.

- Every trial uses the same dense direction bank; the default is 180 readouts
  at 2-degree spacing.
- Independent homoscedastic Gaussian variability is defined over the MT
  preferred-direction population on each patch.
- Projection through overlapping Gaussian readouts creates the candidate
  covariance. Nearby candidates are therefore strongly correlated.
- A separate noise SD is calibrated for each experimental session so the
  dense-bank hard-max comparison reaches 75% correct at that session's fitted
  no-noise coherence threshold.
- IDQ calibration has zero coherent prestep signal because the 7% pedestals
  occur only on coherence-noise trials. IDR uses the session's median physical
  prestep coherence.
- This homoscedastic model does not require an absolute firing-rate baseline.
  A later response-dependent or Poisson-like variant would require one.

Create a new manifested run using the same two experimental summary sources:

```matlab
context = makeAnalysisContext("synthetic", ...
  'RunID',"IDRIDQ_MTInternalNoise_001", ...
  'Seed',1729,'SaveLatents',false);
initializeSyntheticRun(context, ...
  'ModelName',"DenseGaussianMTInternalNoiseCalibration", ...
  'ModelParameters',struct('sigmaMTDeg',37.5, ...
    'sigmaReadoutDeg',20,'candidateSpacingDeg',2, ...
    'nMonteCarlo',2000), ...
  'SourceFiles',[idqSource;idrSource]);

calibration=calibrateIDRIDQMTInternalNoise(context, ...
  'SigmaMTDeg',37.5,'SigmaReadoutDeg',20, ...
  'CandidateSpacingDeg',2,'NMonteCarlo',2000);
```

Tests:

```matlab
results=runtests("testMTInternalNoiseCalibrationCore.m");
assertSuccess(results)
```
