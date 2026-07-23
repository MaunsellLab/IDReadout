# IDR/IDQ MT DOG coarse screen

This package adds the nested unit-center DOG readout

```text
W(delta) = (1+a) G(delta;sigmaCenter) - a G(delta;sigmaSurround)
```

with `sigmaSurround > sigmaCenter`. The `a=0` member is exactly the existing
nonnegative Gaussian. Candidate rectification, upstream directional noise,
hard-max pooling, actual trial inventories, IDQ pedestals, and session-specific
performance calibration are unchanged.

The default screen fixes `sigmaCenter=5 deg`, evaluates surround sigmas
`[20 40 80 160] deg` and magnitudes `[.05 .15 .30 .60]`, and includes one
shared Gaussian baseline: 17 independently manifested member runs.

Each member generates 20 internal-noise choices per actual trial. The IDR gain
model is fitted once to each trial's mean correctness. This is the mean
Bernoulli log likelihood across the 20 binary draws because that likelihood is
affine in the observed outcome. It avoids 20 separate optimizer runs per grid
cell.

Only the observed-SE-weighted discrepancy in `gCQ(theta)` ranks models:

```text
sum_theta ((synthetic_gCQ - observed_gCQ)/observed_SE_gCQ)^2
```

`gCP`, `gNP`, `gNQ`, and performance are reported but held out. IDQ gains are
deferred to ordinary replicated gain recovery for a small number of selected
finalists.

Install every package file in `Tests/Code`, replacing included older files.
Then run:

```matlab
clear makeGaussianReadoutBank ...
  calibrateIDRIDQMTUpstreamNoise ...
  simulateIDRIDQUpstreamNoiseChoices ...
  factorUpstreamDirectionalNoise ...
  fitIDRExpectedChoiceDOGGains ...
  formatDOGRunID ...
  runIDRIDQMTDOGCoarseScreen
rehash

results = runtests([ ...
  "testMTGaussianForwardCore.m", ...
  "testMTNegativeAsymptoteInvariance.m", ...
  "testMTDOGCore.m"]);
assertSuccess(results)
```

Start with a new unique prefix:

```matlab
screen = runIDRIDQMTDOGCoarseScreen( ...
  'RunIDPrefix',"MTDOGCoarse_20260722_001", ...
  'SigmaCenterDeg',5, ...
  'SigmaSurroundDeg',[20 40 80 160], ...
  'SurroundMagnitude',[.05 .15 .30 .60], ...
  'NChoiceDraws',20);
```

Repeat the identical call after a MATLAB interruption. Completed calibration,
choice, and expected-gain stages are loaded; the controller advances to the
first incomplete member. Every member has an independent 200 MB quota and the
summary run has a 100 MB quota.

The same seeds are used to reduce avoidable Monte Carlo variation. Candidate
covariance ranks and bases can differ across DOG parameters, so the draws are
not guaranteed to be exactly paired physical upstream perturbations.

Final summary outputs are:

- `IDRIDQ_MTDOGCoarseScreen_Ranking.csv`
- `IDRIDQ_MTDOGCoarseScreen_IDRGains.csv`
- `IDRIDQ_MTDOGCoarseScreen_Performance.csv`
- `IDRIDQ_MTDOGCoarseScreen_NoiseMetrics.csv`
- `IDRIDQ_MTDOGCoarseScreen.pdf`

The three-page PDF shows the selection criterion, the best DOG versus the
Gaussian baseline for all four IDR gains, and the complete `gCQ` curve family.
