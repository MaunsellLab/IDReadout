# Split-signal model-recovery simulation

This is the inferential stage of the split-signal design analysis. It fixes the
stimulus at 90-degree carrier separation, full prestep pedestals, no opposing
balancing stream, and 50:50, 62.5:37.5, and 75:25 allocations.

Max, population-vector magnitude, and population squared-energy observers all
receive the same centered Gaussian MT response landscape and the same smooth
direction-correlated pre-pooling noise. Each model's 75%-correct threshold is
calibrated separately. Recovery trials span 0.8, 1.0, and 1.2 times threshold.

For each generating model, an independent simulation with correlated
pre-pooling noise calibrates a smooth choice-probability link. Recovery data are
then generated on a separate external-noise bank. Every simulated experiment
refits the decision scale for all candidates and refits the log-sum-exp
temperature for the noisy-max family. Model selection uses BIC by default, so
the extra max-family softness parameter is penalized rather than giving max an
unfair flexibility advantage. This avoids giving any model a unique per-trial
Monte Carlo fingerprint. The output is a full true-model by fitted-model
recovery curve, rather than a test confined to max-generated data.

The correlated pre-pooling noise is integrated directly when calibrating each
observer's threshold and its independent smooth probability link. The link's
decision scale is then re-estimated in every pseudoexperiment. The default has
no lapse process; a lapse term should be introduced only if empirical lapse
rates make it necessary.

The module also retains the five individual 50-ms external-noise frames and
reports frame-resolved psychophysical gains with empirical standard errors.

Run:

```matlab
validateSplitSignalModelRecovery
results = runSplitSignalModelRecovery;
```

The validation run is intentionally tiny and checks structure, not numerical
accuracy. The full run is computationally heavier than the original design
sweep because it integrates pre-pooling noise for every candidate model.
