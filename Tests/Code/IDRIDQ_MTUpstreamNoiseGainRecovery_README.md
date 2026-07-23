# Replicated gain recovery

This stage replaces only trial correctness with each replicated upstream-noise
hard-max choice vector. It retains the experimental rectangular-step
predictors, physical coherence signs, session thresholds, and fitted pooled
psychometric shapes.

The first recovery stage fits:

- IDQ change-side drift and pooled non-drift gains;
- IDQ no-change-side drift and pooled non-drift gains;
- IDR `gCP`, `gCQ`, `gNP`, and `gNQ` independently at every retained offset.

The one 1-degree IDR session remains excluded because only trials represented
in the synthetic choice inventory enter recovery.

The function saves a compact checkpoint after every completed replicate. A
rerun resumes automatically. Final tables and the two-page PDF are written
only after every replicate is complete, and final outputs are never
overwritten.

Install all files together in `Tests/Code`. Run the focused tests:

```matlab
clear recoverIDRIDQMTUpstreamNoiseGains matchSyntheticTrialRows
rehash
results=runtests("testMTGainRecoveryCore.m");
assertSuccess(results)
```

Optional one-replicate smoke test (it checkpoints, then returns):

```matlab
recovery=recoverIDRIDQMTUpstreamNoiseGains(contextUpstream, ...
  'MaxNewReplicates',1);
```

Resume and finish all remaining replicates:

```matlab
recovery=recoverIDRIDQMTUpstreamNoiseGains(contextUpstream);
```

The existing project function `fitIDQNoiseGain.m` and the previously supplied
`analyzeIDRMatchedGains.m` must be on the MATLAB path. The recovery deliberately
uses those same fitters rather than introducing a second likelihood
implementation.
