# Constant negative-asymptote test

This revision tests the nested readout

```text
W(delta) = (1+b) Gaussian(delta;sigmaR) - b
```

which has unit weight at the candidate's preferred direction and asymptote
`-b`. The default `b=0` reproduces the previously installed nonnegative
Gaussian exactly.

Because every `mtPopulationTemplate` is mean-subtracted across the MT bank,
the constant term is orthogonal to every signal and upstream-noise template.
Candidate evidence and covariance are only globally scaled. Positive global
scaling commutes with candidate rectification and hard-max comparison, and is
absorbed by performance calibration. Consequently `b` cannot alter choices or
recovered gains in the current model.

Install the files in `Tests/Code`, replacing the included existing files, then
run:

```matlab
clear makeGaussianReadoutBank ...
  calibrateIDRIDQMTUpstreamNoise ...
  simulateIDRIDQUpstreamNoiseChoices ...
  factorUpstreamDirectionalNoise ...
  analyzeMTNegativeAsymptoteInvariance
rehash

results = runtests([ ...
  "testMTGaussianForwardCore.m", ...
  "testMTNegativeAsymptoteInvariance.m"]);
assertSuccess(results)

summary = analyzeMTNegativeAsymptoteInvariance;
```

The diagnostic should report numerical residuals near floating-point zero and
`choiceAgreement=1` for every tested asymptote. No experimental files are read
and no simulation run is created, because this is an algebraic property of the
forward model rather than an empirical fit.

The next non-degenerate extension is a finite-width negative surround, such as
a DOG. Unlike a constant tail, the surround is not proportional to the all-ones
vector and therefore survives the mean-subtracted MT transformation.
