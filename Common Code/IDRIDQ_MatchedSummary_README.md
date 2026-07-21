# Matched IDQ-IDR summary report

Place these functions on the IDReadout MATLAB path:

- `makeIDRIDQMatchedSummary.m`
- `plotIDRIDQMatchedSummary.m`

The report reads already-fitted files from the domain-specific summary
folders:

- `IDQ/Data/AcrossSessionSummaries`: all `IDQ_*.mat` sources below
- `IDR/Data/AcrossSessionSummaries`: `IDR_SideGainAnalysis.mat`

- `IDQ_AcrossSideSummary.mat`
- `IDQ_SignedNoiseSummary.mat`
- `IDQ_NoChangeSignedNoiseSummary.mat`
- `IDQ_OpponentPoolingSummary.mat`
- `IDQ_OpponentPNormSummary.mat`
- `IDQ_InteractionSummary.mat`
- `IDR_SideGainAnalysis.mat`

Run:

```matlab
plotIDRIDQMatchedSummary
```

The function writes to `Common Code`:

- `Common Code/Data/AcrossSessionSummaries/IDRIDQ_MatchedSummary.mat`
- `Common Code/Plots/AcrossSessionSummaries/IDRIDQ_MatchedSummary.pdf`

The default folders are resolved with the requested-domain form of
`domainFolder`, so no path arguments are needed when the IDReadout MATLAB
Project is open.

The compact MAT file contains fitted and binned report quantities, not the
large duplicated trialwise noise matrices. No psychometric, gain, signed,
interaction, or pooling model is refitted. The one new descriptive
calculation reconstructs IDQ physical candidate margins from the saved
trialwise rectangular-step predictors:

- Change side: `max(nPlus,nMinus) - (stepCoh + nD)`
- No-change side: `max(nPlus,nMinus) - nD`

These are physical-coherence overlap calculations. They are not estimates
of neural candidate overlap, which additionally depends on MT tuning,
divisive normalization, internal variability, and within-patch selection.

Custom paths can be supplied when needed:

```matlab
plotIDRIDQMatchedSummary( ...
  'IDQSummaryFolder', '/path/to/IDQ/Data/AcrossSessionSummaries', ...
  'IDRSummaryFolder', '/path/to/IDR/Data/AcrossSessionSummaries', ...
  'CommonSummaryFolder', '/path/to/Common Code/Data/AcrossSessionSummaries', ...
  'OutputPath', '/path/to/IDRIDQ_MatchedSummary.pdf', ...
  'Visible', 'on');
```
