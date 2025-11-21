function kernelAverage()

  % fileNames = {"IDRTest00", "IDRTest01", "IDRTest02", "IDRTest03"};
  fileNames = {"IDRTest04", "IDRTest05"};
  rootPath = "/Users/maunsell/Desktop/IDR/Data/Kernels/";
 
  [firstPreStepMS] = integralWindowMS();
  for f = 1:length(fileNames)
    fprintf("\nProcessing %s (%d of %d)\n", fileNames{f}, f, length(fileNames));
    load(rootPath + fileNames{f}, "header", "kernels", "kVars", "trialOutcomes");

    if f == 1
      firstStepMS = header.stepMS.data(1);
      firstVFrames = size(kernels, 3);
      avgKernels = zeros(2, 2, firstVFrames);
      avgKInts = zeros(2, 2);
      sumWeights = zeros(2, 2);             % summed variance weighting over sessions
      nHits = zeros(1, 2);
      nTrials = zeros(1, 2);
      sumRWeights = zeros(1, 2);            % we weight R by its inverse variance
      sumWeightedR = zeros(1, 2);                         
      nSessions = 0;
    end
    frameRateHz = header.frameRateHz.data(1);
    preStepMS = header.preStepMS.data(1);       % sometimes header values are replicated vectors
    msPerVFrame = 1000.0 / frameRateHz;
    stepMS = header.stepMS.data(1);
    vFrames = size(kernels, 3);
    if preStepMS ~= firstPreStepMS || stepMS ~= firstStepMS || vFrames ~= firstVFrames
      error('Files based on different trial or kernel lengths.');
    end
    nFileHits = nan(2, 2);
    nFileTrials = nan(2, 2);
    for s = 1:2
      for p = 1:2
        for v = 1:vFrames
          avgKernels(s, p, v) = avgKernels(s, p, v) + kernels(s, p, v) / kVars(s, p);
        end
        sumWeights(s, p) = sumWeights(s, p) + 1.0 / kVars(s, p);
      end
      nFileTrials(s) = length(trialOutcomes{s});
      nTrials(s) = nTrials(s) + nFileTrials(s);
      nFileHits(s) = nFileTrials(s) - sum(trialOutcomes{s});
      nHits(s) = nHits(s) + nFileHits(s);
    end
    [kIntegrals, R, RVar] = kernelIntegral(kernels, kVars, msPerVFrame);
    for s = 1:2
      weight = 1.0 / RVar(s);
      sumWeightedR(s) = sumWeightedR(s) + R(s) * weight;
      sumRWeights(s) = sumRWeights(s) + weight;
    end
    avgKInts = avgKInts + kIntegrals;
    [~, baseName, ~] = fileparts(header.fileName);
    plotKernels(1, baseName, header, kernels, kVars, kIntegrals, R, RVar, nFileHits, nFileTrials);
    nSessions = nSessions + 1;
  end

  % Compute and display averages
  avgR = sumWeightedR ./ sumRWeights;
  avgRSEM = 1 ./ sumRWeights;
  avgKInts = avgKInts ./ nSessions;
  avgKVars = zeros(2, 2);
  for s = 1:2
    for p = 1:2
      avgKernels(s, p, :) = avgKernels(s, p, :) / sumWeights(s, p);
      avgKVars(s, p) = 1.0 / sumWeights(s, p);           
    end
  end
  plotKernels(2, sprintf("%d Session Average", nSessions), header, avgKernels, avgKVars, avgKInts, avgR, avgRSEM, ...
                nHits, nTrials);
end