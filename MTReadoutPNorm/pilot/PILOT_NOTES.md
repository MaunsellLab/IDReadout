# Initial Meetz IDS pilot — October 23, 2026

The Python numerical reference was executed on two available converted raw
recordings, with their existing IDS coherence corrections. The earlier
August 5 file is v7.3 and was not used in this Python pilot; MATLAB's reader
supports that format. Raw recordings are not included in this package.

| Recording | Offset | Valid completed trials | Noisy increment trials |
|---|---:|---:|---:|
| IDR4_Meetz_20260806.mat | 120° | 348 | 173 |
| IDR4_Meetz_20260807.mat | 45° | 1,246 | 779 |
| Total | | 1,594 | 952 |

Both recordings have the imposed coherence histories needed for the
500–750-ms prestep estimate and the 750–1000-ms step estimate. Both use
randomized drift and a corrected 7.071% drift base/probe pedestal.

One internal-noise SD was fitted across these sessions for each pooling
rule. The pilot used 64 field draws per RDK at fitting and 256 draws with a
different seed for reevaluation. Every draw contains independent prestep
and step fields; all RDK pairs are integrated after calculating the signed
pooled changes. Scores are preliminary numerical evaluations, not final
exponent estimates or a model comparison with sampling uncertainty.

| Pooling | Fitted internal SD | Reevaluation log likelihood | Predicted noisy-trial accuracy |
|---|---:|---:|---:|
| p=5 | 5.251 | −511.686 | 75.89% |
| Hard max | 6.620 | −507.498 | 75.32% |

SD is in normalized MT response units, per unit per 250-ms estimate. It is
not a coherence threshold. Observed pooled noisy-trial accuracy is about
77.84%. The SD maximizes choice likelihood; this does not force equality of
the mean predicted and observed accuracies.

| Offset | Observed noisy accuracy | p=5 prediction | Hard-max prediction |
|---|---:|---:|---:|
| 45° | 78.18% | 77.46% | 76.41% |
| 120° | 76.30% | 68.82% | 70.46% |

The common-noise calibration describes the 45° session more closely than
the short 120° session. This is a discrepancy to revisit with additional
sessions, not evidence that a particular extra parameter is needed. No
offset-specific tuning, internal noise, adaptation, or eligibility was added.

Absolute-rate clipping affected 0.0239% of evaluated rate samples for p=5
and 0.1305% for hard max. Neither unit changes nor pooled changes were
rectified.

## Checks actually executed

- Known p-norm values and a zero-response pool.
- Signed negative differences: the smaller decrease wins instead of tying.
- Exact equality of the two RDKs gives 0.5 choice probability.
- Noiseless p=1 fixed-budget population-sum invariance, within quadrature precision.
- A common rotation of responses and noise preserves decisions.
- Pedestals affect the finite-p pooled change after prestep subtraction.
- Direction-grid check with a matched continuous noise realization: mean
  absolute choice-probability change from 1° to 2° was 3.45e−7 for p=5 and
  7.25e−5 for hard max. Maximum changes were 6.10e−5 and 3.66e−4, respectively.
  This check used 128 draws with a separate seed and did not refit noise.

## What remains before the primary gain comparison

1. Add Meetz IDS recordings across the offset range, especially near zero
   and 170–180°, and more 120° sessions. Prefer the complete corrected raw
   IDS session folder or equivalent ProbeSessions containing both windows.
2. Supply the matching current IDS across-session summary MAT files, or
   run plotIDSSummary in the existing project to make them. These freeze
   the exact predictors and psychometric calibration; the new adapter then
   runs the existing gain fitter on observed and predicted outcomes.
3. Execute the MATLAB implementation and its analytic checks in MATLAB.
   MATLAB/Octave was not available in this execution environment. The
   Python reference tests the model mathematics independently; it is not
   a claim that the MATLAB package has already passed a runtime test.
4. Increase integration precision and check refitted SD stability before
   interpreting differences between exponents or reporting confidence limits.

No measured probe gains or recovered gain-curve matches are claimed from
this two-session pilot. The exact current summary structs were not among
the available inputs, so gain recovery was not approximated with a different
regression. An IDR/yoking adapter is a subsequent step.
