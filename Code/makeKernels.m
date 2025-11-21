function makeKernels()

  rootPath = "/Users/maunsell/Desktop/IDR";
  load(rootPath + "/Data/Lablib/IDRTest05.mat", "header", "trials");
  
  % extract values we will want to save
  frameRateHz = header.frameRateHz.data(1);
  % preStepMS = header.preStepMS.data(1);       % sometimes header values are replicated vectors
  % [intStartMS, intDurMS] = integralWindowMS();
  msPerVFrame = 1000.0 / frameRateHz;
  m = round((header.preStepMS.data(1) + header.stepMS.data(1)) / msPerVFrame);  

  prefCohNoisePC = trials{1}.trial.data.prefCohNoisePC;
  probeCohNoisePC = trials{1}.trial.data.probeCohNoisePC;
  kernels = nan(2, 2, m);               % 2 x 2 x m kernel vectors with (inc/dec, pref/probe, vFrame)
  kVars = nan(2, 2);                    % 2 x 2 kernel variance scalar (inc/dec, pref/probe)
  trialOutcomes = cell(1, 2);           % 1 x 2 trial outcome vectors (inc/dec); pref and probe are always the same
  nHits = nan(1, 2);
  nTrials = nan(1, 2);

  for s = 1:2                           % for each coherence step direction (inc/dec)
    [prefMat, probeMat, trialOutcomes{s}] = extractNoiseMatrices(header, trials, s, 0);
    nTrials(s) = length(trialOutcomes{s});
    nHits(s) = nTrials(s) - sum(trialOutcomes{s});     % misses are 1, hits are 0
    [kernels(s, 1, :), kVars(s, 1)] = meanPsychKernel(prefMat, trialOutcomes{s}, prefCohNoisePC);
    [kernels(s, 2, :), kVars(s, 2)] = meanPsychKernel(probeMat, trialOutcomes{s}, probeCohNoisePC);
  end
  [kIntegrals, R, RVar] = kernelIntegral(kernels, kVars, msPerVFrame);
  % display file
  [~, baseName, ~] = fileparts(header.fileName);
  plotKernels(1, baseName, header, kernels, kVars, kIntegrals, R, RVar, nHits, nTrials);
  print(gcf, rootPath + "/Plots/Kernel Plots/" + baseName, '-dpdf');
  save(rootPath + "/Data/Kernels/" + baseName, "header", "kernels", "kVars", "trialOutcomes");
end
