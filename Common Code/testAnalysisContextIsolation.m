function tests = testAnalysisContextIsolation
% testAnalysisContextIsolation  Boundary tests for context-aware paths.
tests = functiontests(localfunctions);
end

function testExperimentalRoots(testCase)
context = makeAnalysisContext("experimental");
projectRoot = string(context.ProjectRoot);
verifyEqual(testCase,string(domainFolder(mfilename('fullpath'),'IDR',context)), ...
  canonicalPath(fullfile(projectRoot,'IDR')));
verifyEqual(testCase,string(domainFolder(mfilename('fullpath'),'IDQ',context)), ...
  canonicalPath(fullfile(projectRoot,'IDQ')));
verifyEqual(testCase,string(domainFolder(mfilename('fullpath'),'Common Code',context)), ...
  canonicalPath(fullfile(projectRoot,'Common Code')));
end

function testTestsRequiresSyntheticContext(testCase)
verifyError(testCase,@() domainFolder(mfilename('fullpath'),'Tests'), ...
  'IDReadout:TestsContextRequired');
context = makeAnalysisContext("experimental");
verifyError(testCase,@() domainFolder(mfilename('fullpath'),'Tests',context), ...
  'IDReadout:TestsContextRequired');
end

function testSyntheticDomainRouting(testCase)
context = makeAnalysisContext("synthetic",'RunID',"isolation_test",'Seed',7);
runRoot = canonicalPath(fullfile(context.ProjectRoot,'Tests', ...
  'Simulation Runs','isolation_test'));
verifyEqual(testCase,string(domainFolder(mfilename('fullpath'),'Tests',context)), ...
  runRoot);
verifyEqual(testCase,string(domainFolder(mfilename('fullpath'),'IDR',context)), ...
  canonicalPath(fullfile(runRoot,'IDR')));
verifyEqual(testCase,string(analysisPath(context,'IDQ','Data','result.mat')), ...
  canonicalPath(fullfile(runRoot,'IDQ','Data','result.mat')));
end

function testSyntheticCannotEscapeRun(testCase)
context = makeAnalysisContext("synthetic",'RunID',"escape_test",'Seed',11);
experimentalPath = fullfile(context.ProjectRoot,'IDR','Data','result.mat');
verifyError(testCase,@() validateAnalysisPath(context,experimentalPath), ...
  'IDReadout:SyntheticPathEscape');
verifyError(testCase,@() analysisPath(context,'IDR','..','escape.mat'), ...
  'IDReadout:UnsafePathComponent');
end

function testExperimentalCannotAccessTests(testCase)
context = makeAnalysisContext("experimental");
testPath = fullfile(context.TestsRoot,'Simulation Runs','any','IDR','Data');
verifyError(testCase,@() validateAnalysisPath(context,testPath), ...
  'IDReadout:ExperimentalTestsAccess');
end

function testSyntheticContextCannotBeForged(testCase)
context = makeAnalysisContext("synthetic",'RunID',"forgery_test",'Seed',13);
context.RunRoot = string(fullfile(context.ProjectRoot,'IDR'));
verifyError(testCase,@() validateAnalysisContext(context), ...
  'IDReadout:InvalidSyntheticRunRoot');
end

function path = canonicalPath(path)
path = string(char(java.io.File(char(path)).getCanonicalPath()));
end
