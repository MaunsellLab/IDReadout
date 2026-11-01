function validateIDJPNormProfile(sourcePath)
% validateIDJPNormProfile  Software smoke check using disjoint trial fragments.
% These fragments are test fixtures, NOT independent experimental sessions.
% Never use their p estimates or bootstrap intervals as scientific results.
A = makeIDJSessionAnalysis(sourcePath);
folder = tempname;
mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
for k = 1:3
  rows = mod((1:height(A.trialTable))'-1,3) == k-1;
  B = A;
  B.trialTable = A.trialTable(rows,:);
  B.changeNoiseByFrameTrial = A.changeNoiseByFrameTrial(:,rows,:);
  B.noChangeNoiseByFrameTrial = A.noChangeNoiseByFrameTrial(:,rows,:);
  B.sourceSessionID = string(A.sourceSessionID) + "_TestFragment" + k;
  reference = abs(B.trialTable.unsignedAllocPC-100) < 1e-6;
  assert(any(reference & B.trialTable.correct) && any(reference & ~B.trialTable.correct));
  sessionAnalysis = B; %#ok<NASGU>
  save(fullfile(folder,char(B.sourceSessionID + "_sessionAnalysis.mat")), ...
    'sessionAnalysis','-v7.3');
end
args = {'SessionFolder',folder,'Animal',A.sessionHeader.animal, ...
  'IncludeSession100',true, ...
  'SaveResult',false,'ExportPDF',false,'MakePlots',false,'Verbose',false};
originalRng = rng;
R = pilotIDJPNormProfile(args{:},'QGrid',[0 .25 .5 .75 1],'NBootstrap',3);
assert(isequal(rng,originalRng),'The caller RNG state was not restored.');
assert(strcmp(R.sourceDomain,'IDJ') && R.nSignalStreams == 3);
assert(isequal(size(R.profileLL),[5 2]));
assert(isequal(size(R.bootstrapSquareLL),[3 2]));
assert(all(isfinite(R.profileLL),'all') && all(isfinite(R.bootstrapSquareLL),'all'));
assert(all(R.squareModel.sigmaBySessionScheme>0,'all'));
peak = max(R.profileLL,[],1)';
square = sum(R.squareModel.sessionLLByScheme,1)';
assert(max(abs(square-R.summaryTable.squareLL)) < 1e-8);
assert(max(abs(peak-square-R.summaryTable.deltaLLBestPNormVsSquare)) < 1e-8);
assert(max(abs(R.summaryTable.relativeLLatSquare+ ...
  R.summaryTable.deltaLLBestPNormVsSquare)) < 1e-8);
% Square is a separate fit; adding q grid points must not change its score.
C = pilotIDJPNormProfile(args{:},'QGrid',[0 .5 1],'NBootstrap',0);
assert(max(abs(C.summaryTable.squareLL-R.summaryTable.squareLL)) < 1e-8);
assert(all(peak>=max(C.profileLL,[],1)'-1e-8));
% All schemes compare precisely the same score sets within a scheme.
assert(isequal(R.nScoredTrialsBySessionScheme(:,1),R.nTrials));
assert(isequal(R.nScoredTrialsBySessionScheme(:,2),R.nDividedTrials));
fprintf('IDJ p-norm smoke check passed on test-only fragments of %s.\n', ...
  char(A.sourceSessionID));
end
