# Validation performed during development

The plain-numeric MAT-v7/ZIP format and reader were exercised with an
independent Python extraction of the available corrected Meetz IDS
recordings from August 6 and August 7, 2026:

- 2 sources, 1,594 valid trials, 952 noisy increment trials.
- All 20 full-trial 50-ms noise bins, physical side coding, and stream
  directions were retained.
- Both direct 250-ms total-coherence averages were reconstructed from the
  binned noise and deterministic baselines within single-precision tolerance.
- A MAT save/load round trip with only one retained trial verified handling
  of squeezed dimensions.
- A deliberately corrupted yoking assignment was rejected by the validator.

The new MATLAB exporter was reviewed but not executed: MATLAB/Octave is not
available in this development environment. The included
`validateMTReadoutExploration` creates small four-task MATLAB fixtures and
tests the real selector, randomized-IDQ filtering and stream reordering,
IDR 179/180-degree yoking conventions, IDJ top-up exclusion/equal-thirds,
plain serialization, temporal reconstruction, and reuse of cached exports.
It writes to a temporary folder and restores the remembered export root.

This distinguishes numerical-format validation on real inputs from MATLAB
runtime validation. No export of the user's complete local dataset, and no
new population-model search, is claimed before the input bundle is returned.
