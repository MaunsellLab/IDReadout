# IDS analysis

This package implements the first offset-specific IDS analysis. It preserves
the IDQ psychometric normalization, leave-one-session-out temporal weighting,
and shared across-session gain fitting.

## Required existing functions

The IDS folder should have access to the existing Common Code and IDQ fitting
helpers, including `domainFolder`, `validFolder`, `convertIDSData`,
`fitIDQInitialSessionThresholds`, `fitIDQAcrossAlignedWeibull`,
`fitIDQSessionWeibullFixedBeta`, `idqWeibullAlphaForThreshold`, and
`idqWeibullP`.

## Pipeline

1. `convertIDSData` creates `Data/FullSessions` files.
2. `makeIDSProbeSessions` creates one raw probe session per unsigned offset.
3. `makeIDSSessionAnalyses` excludes instruction/invalid trials, reconstructs
   the three independent streams, and creates compact session analyses.
4. `plotIDSSummary` creates one across-session MAT file and PDF per offset.

The fixed stream order is drift, drift+offset (Probe0), drift-offset (Probe1).
Psychometric coherence is `abs(stepCohPC-baseCohPC)`, not the endpoint.

With only one session at an offset, `looStep` uses that session's own kernel,
matching the current IDQ fallback. Once two or more sessions exist, each
session's predictor uses a kernel computed from the other sessions only.
