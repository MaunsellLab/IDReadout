# MT readout: initial population p-norm model

This is a separate implementation of the model agreed on October 23, 2026.
It does not change the earlier MaxPoolingFit package or the IDS analyses.

The October 24 local-workflow update is documented in
**LOCAL_QUICKSTART.md**. Run on your hardware and send the compact
`MTReadoutPNormShare.zip`; the large recordings remain local.

## Model

1. Gaussian direction tuning, sigma 37.5 degrees, sampled every 1 degree.
2. Construct each unit's **full positive response** to the recorded coherent
   streams plus the expected isotropic distribution of remaining dots.
   A 100%-coherent preferred stimulus contributes 100 response units. An
   isotropic stimulus contributes about 26.1 response units. Additional
   spontaneous/untuned firing is initially zero and explicitly configurable.
   Actual random-dot angle realizations are not recovered from coherence
   traces; their remaining direction distribution is modeled as isotropic.
3. Average responses separately in 500–750 ms and 750–1000 ms. Both windows
   are 250 ms; integration precedes pooling. No temporal kernel weights the
   forward response in this initial version.
4. Add four independent Gaussian directional noise fields: one per RDK and
   window. All have the same marginal SD and circular correlation length
   (30 degrees initially). No cross-window or cross-RDK correlation is fitted.
5. Floor negative **absolute rates** at zero, then calculate the ordinary
   p-norm across directions. This physical rate bound is not rectification
   of firing-rate changes or pooled differences. The clipped fraction is
   reported so its influence remains visible.
6. For each RDK, subtract its prestep p-norm from its step p-norm. Keep this
   difference signed and choose the RDK with the larger difference. Exact
   ties receive probability 0.5. No additional adaptation, eligibility,
   downstream noise, bias, or lapse parameter is included.

The same p applies to every direction, both RDKs, and both epochs. Comparing
p values is an exploratory model check, not a claim about volitional control.
At p=1 the noiseless uniformly sampled population sum is approximately
constant for a fixed dot budget; this limit is not a linear sum of abstract
stimulus streams. Noise clipping can change that limiting behavior.

## MATLAB usage

Add this directory, not just +mtpnorm, and the existing IDS analysis code
(including your current selectAnalysisFiles.m) to the MATLAB path:

```matlab
validateMTReadoutPNormLocal;
result = runMTReadoutPNormLocal;
```

Choose the IDS or IDS/Data folder once; the launcher remembers it and calls
selectAnalysisFiles on the ProbeSessions folders. The explicit-file
runMTReadoutPNorm interface remains available for programmatic use. See the
quickstart for two-file checks, incremental caching, folder changes, and
gain completion without repeating the forward fit.

The reader accepts converted raw IDS sessions (`trials`) or ProbeSessions
(`probeTrials`). Use either representation of a given trial; overlapping
sources are rejected. It retains all valid completed trials for prediction,
but fits noise only on noisy increments. Instruction, uncertified, incomplete,
and recognized top-up trials are excluded. Trial descriptors must already
contain the established IDS step-coherence correction. The reader never
applies that correction a second time. A correction-metadata presence flag
is included in the input audit; it does not validate an unknown converter.

The reader supports MATLAB v5 and v7.3 files through MATLAB's load. It keeps
actual imposed pre/post coherence fluctuations, positive probe pedestals,
physical RDK side coding, session IDs, and parent trial indices. It requires
the agreed timing and randomized drift, and rejects unphysical coherence
budgets instead of silently repairing them. The current reader is for IDS;
an IDR adapter preserving yoking is still needed for the subsequent task test.

## Noise calibration and numerical precision

One noise SD is fitted across the supplied Meetz IDS sessions for each p,
using their observed binary choices and the imposed stimulus fluctuations.
Noise units are the normalized MT response units above, per unit per
250-ms estimate. They are not behavioral threshold units.

The runner evaluates a coarse SD grid and refines the best neighborhood.
Common integration draws are used across parameter settings. Within each
RDK, pre and post draws are independent. After forming each signed change,
all independent pairs between the two RDK draw sets are compared. The
nominal number of pairs is draws squared; these pairs are not independent
Monte Carlo replicates and should not be counted as such for precision.

Default fit/final draws are 128/512. The final score reevaluates the fitted
SD at higher integration precision; it does not claim the SD has already
been reoptimized at that precision. Increase fit draws and repeat before
interpreting likelihood differences. Compare seeds and the direction grid
before interpreting an estimated p. An upper-bound SD fit, substantial
rate-floor clipping, or a structured per-offset sensitivity mismatch should
be diagnosed before adding model flexibility.

Outputs: trial cache, input audit CSV, fitted SD/likelihood/accuracy CSV,
per-trial choice probabilities, and configuration, in the result MAT file.
The new cache preserves the final prestep window explicitly; the previous
cache's mean over the entire prestep cannot substitute for it.

## Exact recovery through the existing IDS gain analysis

The new package deliberately does not invent a replacement gain regression.
Use the existing Meetz IDS summaries, which contain the frozen predictors
and observed psychometric calibration:

```matlab
% Run in the current IDS project environment:
summaries = plotIDSSummary('Animal','Meetz','IncludeIDQ',false, ...
    'NoisePredictorType','looStep');
recovered = recoverMTReadoutPNormGains(result,summaries);
save('/path/to/PNormOutput/MTReadoutPNormGains.mat','recovered','-v7.3');
```

You can alternatively pass these structs as 'Summaries' to the main runner.
All trials in a summary must have matching forward predictions; matching
uses sourceSessionID/parentTrialIdx, and checks step size, changed side, and
observed correctness. Supply the same sessions and increment selection to
both analyses. The adapter calls your current fitIDSNoiseGain with the exact
existing predictors, session thresholds, Weibull shape, and lapse. Its only
change is replacing correct/error outcomes by model p(correct). This gives
the expected Bernoulli log likelihood without extra choice-sampling noise.

Observed gains are refitted with the same call. Probe coefficients use the
existing separate and pooled-probe conventions and normalization. Do not
interpret the predicted fit's Hessian CIs as uncertainty in the forward
model parameters; session/bootstrap uncertainty is a later analysis.
The current IDS project functions and their toolboxes must be on the path.

## Validation and pilot

`validateMTReadoutPNorm.m` contains MATLAB analytic checks. MATLAB/Octave was
not available in the execution environment, so these MATLAB files have not
been executed there. `pilot/reference_pilot.py` is an independently written
Python numerical reference of the same model and was executed on the
converted Meetz August 6 (120 degrees) and August 7 (45 degrees) recordings.
Its random-number generator differs from MATLAB; individual draws will not
match. See `pilot/PILOT_NOTES.md` for actual results and limitations.

The raw recordings and existing IDS analysis functions are not bundled.
Additional Meetz sessions across near-zero and opponent offsets, plus the
matching current IDS summary structs, are used locally for the primary gain
test; these large inputs do not need to be transferred. The new local
selector/cache/batching/export MATLAB integration checks are supplied in
validateMTReadoutPNormLocal.m and require execution in MATLAB.
