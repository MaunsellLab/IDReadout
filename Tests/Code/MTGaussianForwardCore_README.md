# Gaussian MT forward-model baseline

This first diagnostic deliberately stops before max pooling and behavioral
choice generation. It documents what follows from:

- mean-subtracted signed MT population templates;
- fixed MT tuning sigma, default 37.5 degrees;
- a nonnegative Gaussian readout with width `sigmaReadoutDeg`;
- no surround, additive baseline, candidate rectification, or pooling.

Create a manifest-backed synthetic run and execute the width sweep:

```matlab
context = makeAnalysisContext("synthetic", ...
    'RunID',"MTGaussianWidth_001", ...
    'Seed',1729, ...
    'MaxRunBytes',2e9);

initializeSyntheticRun(context, ...
    'ModelName',"GaussianMTReadoutWidthSweep", ...
    'ModelParameters',struct( ...
        'sigmaMTDeg',37.5, ...
        'sigmaReadoutDeg',5:0.5:45));

summary = analyzeMTGaussianReadoutWidths(context);
```

The run receives a MAT summary, CSV width table, and PDF diagnostic beneath
its `Common Code/Data/Simulation` and `Common Code/Plots/Simulation`
folders. The principal table reports the preferred response, null response,
and their absolute ratio for every readout width. The offset curves show the
negative far-offset behavior of a fixed Gaussian readout.

Run the pure forward-core tests with:

```matlab
results = runtests("testMTGaussianForwardCore.m");
assertSuccess(results)
```

The next model stage will retain this baseline and add candidate banks,
first without and then with candidate-level rectification, followed by sum,
finite-p, and hard-max pooling.
