# Project adapter boundary

The included unified reader has signature:

```matlab
fragment = readMaxPoolingProbeSession(sourcePath)
```

The returned fragment must satisfy the cache contract in the top-level
README. It can then be used as:

```matlab
addpath('/path/to/MaxPoolingFit/adapters')
cache = maxpoolfit.buildCache(files,@readMaxPoolingProbeSession,outputPath, ...
    'Replace',true);
```

Before running the full extraction, validate one representative file from
each domain by confirming:

1. Trial selection and choice/change-side fields match existing summaries.
2. Step magnitude is endpoint minus prestep coherence.
3. Framewise stream values are correctly identified as total coherence or
   zero-centered noise.
4. Probe pedestals are included exactly once in `stepCohPC`.
5. `preCohPC` contains one prestep coherence per displayed stream, aligned
   with `streamDirDeg`; drift and probe pedestals are included exactly once.
6. IDR paired probes remain separate direction streams even though their
   temporal fluctuations are yoked.
7. IDQ/IDS independent streams remain independent.
8. The corrected IDS trial descriptors are treated as authoritative.
9. Kernel-weighted rectangular predictors reproduce current analyses before
   fitting the max-pooling model.

The current IDR `makeProbeSessions` products contain noise trials only. They
are sufficient for the primary fit, but an IDR FullSession reader will be
needed if later reconstruction of IDR no-step-noise psychometrics must rely on
the new cache alone. IDQ/IDS ProbeSessions retain their no-step-noise trials.

The adapter follows the established source-code conventions and should first
be run on the three representative ProbeSessions. Inspect each returned
`extractionAudit`: noisy traces should span negative and positive values and
have a mean near zero. Then compare its rectangular stream predictors with the
current session analyses before extracting the full dataset.

Historical ProbeSessions use two `trialIdx` contracts. The adapter detects
whether `trialIdx` indexes every `probeTrials` entry or only the valid analysis
trials, resolves the parent indices accordingly, and records the detected
contract in `extractionAudit.trialIdxContract`.

`readMaxPoolingProbeSession` declares an explicit extractor version. The
incremental cache manifest stores this version for every source file. Changing
the extractor version automatically marks older cached files for re-extraction,
preventing a mixed-schema cache after extraction logic changes.
