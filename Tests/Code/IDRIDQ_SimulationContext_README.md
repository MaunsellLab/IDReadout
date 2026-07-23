# Experimental and synthetic analysis contexts

The context layer makes every data access explicitly experimental or part
of one named synthetic run. Existing calls to `domainFolder` remain
experimental and backward-compatible.

## Experimental analysis

```matlab
context = makeAnalysisContext("experimental");
idrRoot = domainFolder(mfilename('fullpath'),"IDR",context);
```

The historical forms still work:

```matlab
callerRoot = domainFolder(mfilename('fullpath'));
idrRoot = domainFolder(mfilename('fullpath'),"IDR");
```

## Synthetic run

```matlab
context = makeAnalysisContext("synthetic", ...
    'RunID',"MTMax_20260721_001", ...
    'Seed',1729, ...
    'MaxRunBytes',2e9, ...
    'SaveLatents',false);

manifest = initializeSyntheticRun(context, ...
    'ModelName',"MTMax", ...
    'ModelParameters',parameters, ...
    'SourceFiles',[idrSourceFile; idqSourceFile], ...
    'CodeVersion',gitCommit);
```

The resulting layout is:

```text
Tests/Simulation Runs/MTMax_20260721_001/
  manifest.mat
  IDQ/Data/
  IDQ/Plots/
  IDR/Data/
  IDR/Plots/
  Common Code/Data/
  Common Code/Plots/
```

Resolve all paths through the context:

```matlab
outputPath = analysisPath(context,"IDR", ...
    "Data","AcrossSessionSummaries","IDR_SideGainAnalysis.mat");
```

With a synthetic context, requests for `IDR`, `IDQ`, or `Common Code` are
routed beneath that run. A request for `Tests` returns the named run root,
not the unrestricted top-level Tests folder. `Tests` cannot be requested
without a validated synthetic context.

`validateAnalysisPath` rejects cross-origin access. The only intentional
experimental read in a synthetic workflow occurs during simulation import;
`initializeSyntheticRun` records source paths, sizes, modification times,
and SHA-256 hashes in the manifest.

Use `assertSimulationRunCapacity` before a material write and after it with
zero additional bytes. The default quota is 2 GB per run. Large MT latent
responses should not be saved unless `SaveLatents` was explicitly enabled.

Run the boundary tests with:

```matlab
results = runtests("testAnalysisContextIsolation.m");
assertSuccess(results)
```
