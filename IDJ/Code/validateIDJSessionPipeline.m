function validateIDJSessionPipeline(sourcePath)
% validateIDJSessionPipeline  Check extraction, stream counts, and plotting.
% Run once for each sample FullSession after installing Code on MATLAB's path.
A = makeIDJSessionAnalysis(sourcePath);
assert(idjHasTopUpExclusion(A));
assert(all(~A.trialTable.topUpTrial));
T = A.trialTable;
S = A.streamTable;
n = height(T);
assert(A.nSignalStreams == 3 && height(S) == 3*n);
assert(A.nCorrect + A.nWrong == n);
assert(sum(A.stepTable.nAnalysisTrials) == n);
assert(all(accumarray(S.analysisTrialRow, 1) == 3));
assert(isequal(S.sourceTrialIndex, T.sourceTrialIndex(S.analysisTrialRow)));
assert(isequal(size(A.changeNoiseByFrameTrial), [30 n 3]));
assert(isequal(size(A.fullChangeNoiseByFrameTrial), [120 n 3]));
assert(isequal(size(A.fullNoChangeNoiseByFrameTrial), [120 n 3]));
assert(all(isfinite(A.fullChangeNoiseByFrameTrial), 'all'));
assert(all(isfinite(A.fullNoChangeNoiseByFrameTrial), 'all'));
assert(all(ismember(S.role, ["major" "minor" "equal"])));
for allocation = [100 70 50]
  trialRows = abs(T.unsignedAllocPC - allocation) < 1e-6;
  sampleRows = trialRows(S.analysisTrialRow);
  assert(sum(sampleRows & S.role == "major") == sum(trialRows));
  assert(sum(sampleRows & S.role == "minor") == 2*sum(trialRows));
end
equalTrials = abs(T.unsignedAllocPC - 100/3) < 1e-6;
assert(all(T.highStreamIndex(equalTrials) == 0));
assert(sum(S.role == "equal") == 3*sum(equalTrials));
assert(all(abs(S.allocationPC(S.role == "equal") - 100/3) < 1e-6));

% Verify noise ordering directly against the source, including the third stream.
raw = load(sourcePath, 'trials');
for row = unique([1 n])
  tr = raw.trials{T.sourceTrialIndex(row)};
  for side = {'change', 'noChange'}
    prefix = side{1};
    times = double(tr.([prefix 'TimesMS']).data(:));
    field = ['full' upper(prefix(1)) prefix(2:end) 'NoiseByFrameTrial'];
    noise = A.(field);
    for stream = 1:3
      values = double(tr.(sprintf('%sNoise%dCohsPC', prefix, stream-1)).data(:));
      for frame = 1:numel(A.fullTimeMS)
        absoluteTime = A.fullTimeMS(frame) + double(A.sessionHeader.preStepMS);
        event = find(times <= absoluteTime + 1e-6, 1, 'last');
        assert(abs(double(noise(frame,row,stream)) - values(event)) < 1e-5);
      end
    end
  end
end

[summary, fig] = plotIDJSessionSummary(A, 'ExportPDF', false, 'Visible', 'off');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
assert(summary.freeLapseFit.nTrials == n);
assert(summary.zeroLapseFit.nTrials == n);
assert(sum(summary.freeLapseFit.thresholdTable.nTrials) == n);
assert(all(isfinite(summary.zeroLapseFit.predictedP)));
assert(all(isfinite(summary.freeLapseFit.predictedP)));
assert(summary.freeLapseFit.negLogLikelihood <= ...
  summary.zeroLapseFit.negLogLikelihood + 1e-3);
fprintf('IDJ validation passed: %s (%d analysis trials, %d stream samples).\n', ...
  char(A.sourceSessionID), n, height(S));
end
