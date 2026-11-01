# Compact inputs for MT readout model exploration

Add this folder and the existing project `selectAnalysisFiles.m` to MATLAB's
path, then run:

```matlab
exportMTReadoutExploration
```

Choose the common **IDReadout project root**, containing IDS, IDQ, and IDJ.
Historical IDR data may be in that root's own Data directory, or in an IDR
subdirectory. The root is remembered. The exporter defaults to **Meetz** and
**50-ms bins over the entire prestep and step**. It imposes no p and performs
no population simulation, noise calibration, or likelihood optimization.

Send **MTReadoutExplorationOutput/MTReadoutExplorationInputs.zip**. The command
prints its location. Check the printed coverage table: randomized IDQ at
120 degrees should now contribute additional trials. `Coverage.csv` and
`ExportNotes.csv` identify missing tasks, excluded fixed-IDQ records, or
files that could not be exported. Partial coverage is explicit.

## Existing analyses and large files

The exporter calls your current selector on real source directories, using
its filename/path exclusions. It then filters animal and source task by
their headers. `LoadMetadata=false` is deliberate: legacy compact products
have different header layouts, so these few metadata fields are read in
the export adapter. No alternative hard-exclusion list is installed.

For each selected ProbeSession/FullSession it prefers the corresponding
existing compact SessionAnalysis when that product is at least as recent
as the source and contains the required full noise histories. Otherwise it
reads the relevant source variables one file at a time. Existing two-window
simulation caches cannot supply the full temporal histories; they are not
silently substituted for those histories. The first export can therefore
take time for sources lacking suitable compact analyses.

Each completed source shard is cached immediately. Repeat exports reuse
unchanged shards. Source/compact paths, sizes, modification times, bin size,
task, and adapter version guard reuse. The ZIP includes only the current
selection, not every old cached shard. Nothing in the original data or
existing analyses is changed, and missing analyses are not regenerated.

## Task contracts

- **IDS:** drift, +probe, -probe streams; corrected trial-descriptor base
  and endpoint coherence are used. Pedestals are added exactly once.
- **IDQ:** only `randomBaseDir=true` records. Native absolute streams are
  reordered according to the actual trial dirIndex into drift, +120, -120.
  Already-normalized IDS-compatible IDQ records are also supported, with
  their original IDQ source label retained. Fixed-base IDQ is excluded.
  Equivalent canonical/legacy whole-session representations are deduplicated.
- **IDJ:** all three recorded directions, allocations, endpoints, and noise
  streams are retained. Top-ups are excluded. Compact IDJ inputs must have
  the current top-up exclusion contract; otherwise the original trial flags
  are used. Raw 0.33330/0.33335/0.33335 fractions remain unchanged, with a
  separate nominal equal-thirds label.
- **IDR:** physical patch assignment, preferred direction, base/endpoint
  levels, and the actual yoked probe trace are retained. Paired probes are
  two separately directed streams sharing the recorded trace. A 180-degree
  single probe is not duplicated. ProbeSessions can contain only noisy
  trials; noise-free IDR controls absent from them are not fabricated.

Animal, source task, parent-session identity, source trial index, physical
choice/change side, observed correctness, step-noise flag, directions,
pedestals, and allocations travel with every source shard. No choices are
resimulated. Unknown metadata or ambiguous trial indices produce a recorded
export error rather than a guessed interpretation. Recorded signed coherence
values are preserved without clipping or fixing a dot budget.

## Numerical format

`sources/*.mat` are ordinary compressed **MAT v7** files containing a plain
`source` struct (plus a local cache fingerprint). No MATLAB table, string,
datetime, opaque object, eye trace, or spike record is needed to read the
stimulus inputs in Python. Key arrays use these axis orders:

| Array | Shape / meaning |
|---|---|
| noise_pc | trial × physical RDK × stream × time bin; zero-centered recorded noise |
| deterministic_pre_pc | trial × RDK × stream; actual prestep levels including pedestals |
| deterministic_step_pc | trial × RDK × stream; actual step levels including pedestals |
| pre_window_total_pc | trial × RDK × stream; direct final-250-ms total coherence average |
| step_window_total_pc | trial × RDK × stream; direct 250-ms step total coherence average |
| directions_deg | trial × RDK × stream; actual stimulus directions |
| stream_mask | active streams; padded inactive entries are NaN |
| yoking_group | shared trace groups within each RDK; not a cross-patch correlation assumption |
| allocation_fraction | original three stream fractions |
| nominal_allocation_fraction | equal-thirds labeling separately from recorded values |

Side indices are **1 or 2**. Stream slots are padded to three. Full-trial
noise plus the appropriate deterministic epoch level reconstructs total
coherence. Bin edges and native frame counts are supplied; 50-ms means are
an approximation for within-bin nonlinear temporal models. The direct
two-window averages retain the current model's integration convention.

`analyses/*.mat` preserve existing saved predictors, psychometric fits,
gain estimates, and model kernels using plain columns/structs. Saved animal
and task selection is retained as-is: an All-animal summary is not relabeled
as Meetz. With trial inputs available, targets can be re-estimated here when
needed. Relevant analysis `.m` files found on the MATLAB path are copied to
`analysis_code` for interpreting those definitions.

## Options and checks

```matlab
% Optional small integration check using synthetic sources only:
validateMTReadoutExploration

% Export all available selected sources (default):
exportMTReadoutExploration

% Use a folder chooser again:
exportMTReadoutExploration('ChooseFolder',true)

% Preserve every native video frame instead of 50-ms bins:
exportMTReadoutExploration('BinMS',0)
```

`Animal='All'`, `Domains`, `ProjectRoot`, `OutputFolder`, `Rebuild=true`,
`IncludeAnalyses=false`, and `MaxFiles` are also supported. The defaults
export all four requested task families for Meetz. A limited MaxFiles run
is explicitly marked as a subset in Schema.json.

The Python loader/validator is in `read_exploration.py`. MATLAB/Octave is not
available in the development environment; the supplied MATLAB integration
check must be executed in your MATLAB environment. The numerical bundle
format and two-window reconstruction are independently checked in Python
using available real converted IDS recordings.
