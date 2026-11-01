# Preliminary split-signal simulation results

These values come from an independent numerical mirror used to check the MATLAB
implementation and should be treated as design guidance, not final power
calculations. The primary case used 90-degree carriers, 7.071% prestep
pedestals, five independent 50-ms binary-noise frames, sigma=37.5-degree
centered MT tuning, 30-degree-correlated internal noise with marginal SD 8, a
4/60 adaptation fraction, and no opposing balancing pedestal.

## Feasibility

The coarse threshold grid placed 75%-correct total coherence at approximately:

- 50:50 allocation: 38-40%
- 62.5:37.5 allocation: about 36%
- 75:25 allocation: approximately 30-32%

Thus the equal split appears measurable, but its threshold may be roughly twice
the current single-direction IDS/IDQ threshold. Coherence-budget feasibility
must be checked against the exact dot-allocation code before committing to the
design.

At the nearest-to-threshold grid points, the simulated losing/winning noise-gain
ratio was about 0.09 for 62.5:37.5 and near zero for 75:25. Adding a constant
opposing prestep pedestal changed performance by less than one percentage point
in this particular adaptation model. This small effect is assumption-dependent:
the balancing stream adds a third adapted sector even though it zeros the
physical vector sum.

## Initial model-separation estimate

The noisy-max family was represented by a log-sum-exp summary with its
temperature selected from a fixed grid. Joint fits used all three allocations
and nine total-coherence levels, with one intercept and slope per model. With
2,000 usable trials per session and DeltaLL >= 5 as the criterion, a 500-fold
trial bootstrap reached 90% power at approximately:

- linear sum: 1 session
- vector magnitude: 2 sessions
- summed squared responses: 4 sessions

These estimates are optimistic. They omit lapses, across-session variation,
uncertainty in the sensory front end, and model misspecification. Their useful
message is the ordering: linear pooling should be easy to reject, vector
magnitude harder, and squared energy likely to determine the amount of data
needed.

## Interpretation

The gain measurement carries most of the max-versus-energy leverage. At
90-degree separation, vector magnitude and squared energy require positive gain
from a positive losing component; noisy max drives that gain toward zero. At
120 degrees, vector magnitude itself can generate zero or negative losing gain,
making that geometry less diagnostic despite its practical similarity to IDQ.
