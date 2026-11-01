# Split-signal max-pooling design simulation

This module compares candidate stimulus designs for a coherence-increment task in
which every signal trial increases coherence in two direction streams. The total
mean increment is held fixed while its allocation changes. Independent binary
noise is present only during the five 50-ms step frames and is integrated before
the direction readout.

The primary simulation includes a dense 1-degree MT bank, centered Gaussian
coherence slopes (sigma 37.5 degrees), smooth direction-correlated internal noise
(30-degree correlation length), and hard max pooling between direction
mechanisms. The competing summaries are scalar linear sum, vector magnitude,
and summed squared stream responses. Internal noise makes a deterministic max
score inadequate, so the max family uses a shared log-sum-exp temperature chosen
from a prespecified grid. Each model is given a fitted intercept and slope before
likelihoods are compared.

Pedestal variables are explicit. The sweep compares full, half, and absent
prestep pedestals; 5% and 7.071% step-noise half-amplitudes; optional opposing
prestep motion that zeros the physical vector sum; 90- and 120-degree carrier
separations; and no versus modest firing-dependent adaptation. The balancing
stream remains constant during the step and therefore affects the simulated max
readout only through its prestep adaptation landscape.

When the prestep pedestal is smaller than the step pedestal, the pedestal rise
is applied to both RDKs. The reported `TotalStepCoherences` is the total mean
changed-side increment, including the two pedestal rises. The coherence left for
unequal allocation is therefore `totalStep - 2*commonPedestalRise`. This makes
the reduced-prestep design tradeoff explicit rather than treating the extra
pedestal as free signal.

Run:

```matlab
validateSplitSignalSimulation
summary = runSplitSignalDesignSweep;
```

For an initial faster run:

```matlab
summary = runSplitSignalDesignSweep( ...
    'NDesignTrials', 5000, ...
    'NPowerReplicates', 50, ...
    'MaxSessions', 10);
```

Outputs are written to `SplitSignalDesign_Output`. The CSV contains every
condition; the MAT file preserves trial-summary structures and power curves.
Session estimates are planning estimates, not a stopping rule. They assume
2,000 usable trials per session and should be rerun with the empirical session
yield and lapse behavior before data collection.

Important interpretation: the simulator uses exact template cancellation of the
mean prestep response. Persistent prestep pedestals affect the primary readout
through adaptation and physical stimulus constraints. Any increase needed to
reach the step pedestal is modeled explicitly as a common step on both RDKs.
This assumption is exposed by `TemplateCancellationFraction` and should be
revisited when the exact stimulus implementation is fixed.
