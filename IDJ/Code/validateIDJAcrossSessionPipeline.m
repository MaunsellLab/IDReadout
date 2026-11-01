function summary = validateIDJAcrossSessionPipeline(varargin)
% validateIDJAcrossSessionPipeline  Run IDJ fits and plots without saving.
% Inputs are selectAnalysisFiles options, for example 'Animal','All'.
result = fitIDJModelRecoveryAcrossSessions(varargin{:}, ...
  'MakePlots', false, 'SaveResult', false);
assert(strcmp(result.sourceDomain,'IDJ') && result.nSignalStreams == 3);
assert(all(abs(sum(abs(result.kernelByFrameSession),1)-1) < 1e-6));
assert(height(result.predictionTable) == sum(result.sessionTable.nFitTrials));
P = result.predictionTable{:,{'pL1','pL2','pLinf','pSquare'}};
assert(all(isfinite(P),'all') && all(P>=0 & P<=1,'all'));
[summary, figures] = plotIDJAcrossSessionBehavior(varargin{:}, ...
  'ModelResult',result,'ExportPDF',false,'Docked',false,'Visible','off');
cleanup = onCleanup(@() closeFigures(figures)); %#ok<NASGU>
expected = [0;15;25;100/3;50;70;100];
assert(max(abs(summary.signalFractionPC(:)-expected)) < 1e-6);
for group = 1:numel(summary.groups)
  K = summary.groups(group).kernel;
  assert(isequal(size(K.coefficients),[7 120 2]));
  assert(isequal(size(K.predictedRelativeGain),[7 2 4]));
  assert(isequal(K.nStreamSamples,K.nTrials.*[2;2;2;3;1;1;1]));
  assert(isequal(K.nTrials(1:3),K.nTrials([7 6 5])));
  if K.referenceGainAvailable
    assert(abs(K.relativeGain(7,1)-1) < 1e-6);
  end
end
fprintf('IDJ across-session validation passed: %d sessions, %d groups.\n', ...
  height(summary.sessionThresholdTable),numel(summary.groups));
end

function closeFigures(figures)
names = fieldnames(figures);
for k = 1:numel(names)
  H = figures.(names{k});
  if isgraphics(H.functions), close(H.functions); end
  if isgraphics(H.kernels), close(H.kernels); end
end
end
