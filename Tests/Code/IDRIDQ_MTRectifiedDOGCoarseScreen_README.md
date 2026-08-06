# Rectified-MT DOG coarse screen

This package adds a discrete MT increment-template mode to the existing DOG
screen. It does **not** impose a zero absolute firing-rate floor. In rectified
mode, the unit positive-coherence template for each motion direction is

```text
T_rectified(direction) = max(T_signed(direction),0)
```

before its signed physical coherence coefficient is applied. Thus negative
physical noise remains negative. The same rectified directional templates are
used for experimental stimulus components and upstream internal noise.

The model sequence is:

```text
rectified MT direction increments
  -> unit-center DOG readout
  -> candidate zero rectification
  -> within-patch hard max
  -> patch-max comparison
```

The signed mode remains the default and is unchanged. `SurroundMagnitude=0`
is the rectified-MT Gaussian baseline. The default DOG grid, session-specific
performance calibration, actual trial inventories, 20 choice draws per trial,
fractional expected-choice gain fit, and gCQ-only selection criterion are the
same as in the completed signed-template screen.

Install every package file in `Tests/Code`, replacing included older files,
then run:

```matlab
clear mtPopulationTemplate ...
  makeGaussianReadoutBank ...
  evaluateDenseMTStepCandidates ...
  factorUpstreamDirectionalNoise ...
  calibrateIDRIDQMTUpstreamNoise ...
  simulateUpstreamChoiceCondition ...
  simulateIDRIDQUpstreamNoiseChoices ...
  fitIDRExpectedChoiceDOGGains ...
  runIDRIDQMTDOGCoarseScreen
rehash

results = runtests([ ...
  "testMTGaussianForwardCore.m", ...
  "testMTNegativeAsymptoteInvariance.m", ...
  "testMTDOGCore.m"]);
assertSuccess(results)
```

Start the directly matched rectified-template screen with a new prefix:

```matlab
screenRectified = runIDRIDQMTDOGCoarseScreen( ...
  'RunIDPrefix',"MTRectifiedDOGCoarse_20260723_001", ...
  'MTTemplateMode',"rectified", ...
  'SigmaCenterDeg',5, ...
  'SigmaSurroundDeg',[20 40 80 160], ...
  'SurroundMagnitude',[.05 .15 .30 .60], ...
  'NChoiceDraws',20);
```

Repeat the identical call after interruption. Each of the 17 model members is
an independently manifested, resumable synthetic run. No completed signed-mode
run is modified or reused as rectified data.

Rectified-mode summary outputs have distinct names:

- `IDRIDQ_MTRectifiedDOGCoarseScreen_Ranking.csv`
- `IDRIDQ_MTRectifiedDOGCoarseScreen_IDRGains.csv`
- `IDRIDQ_MTRectifiedDOGCoarseScreen_Performance.csv`
- `IDRIDQ_MTRectifiedDOGCoarseScreen_NoiseMetrics.csv`
- `IDRIDQ_MTRectifiedDOGCoarseScreen.pdf`

Only IDR gCQ enters the ranking. gCP, gNP, gNQ, and performance remain held
out. IDQ gain recovery remains deferred until the rectified and signed screens
identify models worth full replicated recovery.
