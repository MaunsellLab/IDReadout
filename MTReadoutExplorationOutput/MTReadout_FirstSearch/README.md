# MT readout — completed first exploration batch

Start with `MTReadout_FirstSearch.md` and `MTReadout_FirstSearch.png`.
The PDF is the vector version of the three-panel figure.

All 2,006 coarse settings, eight all-trial IDS/IDQ confirmations, task-transfer
checks, source audits, and code are included. No joint solution is claimed.

To rerun, supply the original MTReadoutExplorationInputs.zip via the
MT_INPUT_BUNDLE environment variable and follow the report commands.
Required packages: numpy, scipy, matplotlib. No MATLAB/Octave is needed.
The derived caches are omitted and rebuild automatically. The original
MATLAB analysis snapshots are included for provenance, not execution.

Numerical logs use Python JSON conventions: Infinity denotes a hard max;
NaN denotes an undefined diagnostic ratio. Confirmed gain estimates are finite.

The IDR transfer gains are exploratory window-noise/probit diagnostics,
not exact replications of the unavailable historical IDR gain analysis.
