function results = runSplitSignalModelRecovery(varargin)
% RUNSPLITSIGNALMODELRECOVERY  Fair recovery of max, vector, and energy models.
%
% results=runSplitSignalModelRecovery();

cfg=splitSignalRecoveryDefaults(); cfg=parseOptions(cfg,varargin{:});
if cfg.SaveOutputs && ~exist(cfg.OutputDir,'dir'), mkdir(cfg.OutputDir); end
nModels=numel(cfg.PoolingModels); nAlloc=numel(cfg.AllocationFractions);

calibrationBank=makeSplitSignalTrialBank(cfg.NCalibrationTrials,cfg,cfg.Seed);
thresholds=nan(nModels,nAlloc);
for m=1:nModels
    for a=1:nAlloc
        if cfg.Verbose
            fprintf('Calibrating %s, allocation %.3f\n', ...
                cfg.PoolingModels{m},cfg.AllocationFractions(a));
        end
        thresholds(m,a)=calibrateSplitSignalModelThreshold(calibrationBank, ...
            cfg.PoolingModels{m},cfg.AllocationFractions(a),cfg);
    end
end
if any(~isfinite(thresholds(:)))
    error('At least one model failed to reach target performance within the coherence limit.');
end

recoveryCell=cell(1,nModels); gainsCell=cell(1,nModels);
for trueModel=1:nModels
    if cfg.Verbose
        fprintf('Model recovery for true %s observer\n',cfg.PoolingModels{trueModel});
    end
    design=makeSplitSignalRecoveryDesign(thresholds(trueModel,:),cfg, ...
        cfg.Seed+100*trueModel);
    scores=computeSplitSignalScoreSet(design,cfg);
    link=calibrateSplitSignalTrueLink(thresholds(trueModel,:),trueModel,cfg, ...
        cfg.Seed+1000*trueModel);
    trueProbability=applySplitSignalTrueLink(link,scores,cfg);
    recoveryCell{trueModel}=simulateSplitSignalScoreRecoveryCurve(trueProbability, ...
        scores,trueModel,cfg,cfg.Seed+2000*trueModel);
    gainsCell{trueModel}=estimateSplitSignalRecoveryGains(design,trueProbability, ...
        cfg,cfg.Seed+3000*trueModel);
end
recovery=[recoveryCell{:}]; gains=[gainsCell{:}];

results=struct('config',cfg,'thresholds',thresholds, ...
    'recovery',recovery,'gains',gains);
if cfg.CreateFigures, results.figures=plotSplitSignalModelRecovery(results,cfg); end
if cfg.SaveOutputs
    save(fullfile(cfg.OutputDir,'SplitSignalModelRecovery.mat'),'results','-v7.3');
    writeThresholdCSV(results,fullfile(cfg.OutputDir,'SplitSignalRecoveryThresholds.csv'));
    writeRecoveryCSV(results,fullfile(cfg.OutputDir,'SplitSignalRecoveryProbabilities.csv'));
    writeGainCSV(results,fullfile(cfg.OutputDir,'SplitSignalRecoveryGains.csv'));
    if cfg.CreateFigures
        saveas(results.figures.thresholds,fullfile(cfg.OutputDir,'SplitSignalRecoveryThresholds.pdf'));
        saveas(results.figures.recovery,fullfile(cfg.OutputDir,'SplitSignalModelRecovery.pdf'));
        saveas(results.figures.gains,fullfile(cfg.OutputDir,'SplitSignalRecoveryGains.pdf'));
    end
end
end

function writeThresholdCSV(results,path)
cfg=results.config; fid=fopen(path,'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'generating_model,allocation,threshold\n');
for m=1:numel(cfg.PoolingModels)
    for a=1:numel(cfg.AllocationFractions)
        fprintf(fid,'%s,%.6g,%.6g\n',cfg.PoolingModels{m}, ...
            cfg.AllocationFractions(a),results.thresholds(m,a));
    end
end
end

function writeRecoveryCSV(results,path)
cfg=results.config; fid=fopen(path,'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'sessions,true_model,fitted_model,selection_probability\n');
for t=1:numel(cfg.PoolingModels)
    for s=1:cfg.MaxRecoverySessions
        for f=1:numel(cfg.PoolingModels)
            fprintf(fid,'%d,%s,%s,%.6g\n',s,cfg.PoolingModels{t}, ...
                cfg.PoolingModels{f},results.recovery(t).probability(s,f));
        end
    end
end
end

function writeGainCSV(results,path)
cfg=results.config; fid=fopen(path,'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,['generating_model,allocation,stream1_gain,stream1_se,' ...
    'stream2_gain,stream2_se,equivalent_sessions\n']);
for m=1:numel(cfg.PoolingModels)
    g=results.gains(m);
    for a=1:numel(cfg.AllocationFractions)
        fprintf(fid,'%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n', ...
            cfg.PoolingModels{m},cfg.AllocationFractions(a), ...
            g.stream1Mean(a),g.stream1MeanSE(a), ...
            g.stream2Mean(a),g.stream2MeanSE(a),g.equivalentSessions);
    end
end
end

function cfg=parseOptions(cfg,varargin)
if mod(numel(varargin),2)~=0, error('Options must be name/value pairs.'); end
fields=fieldnames(cfg);
for k=1:2:numel(varargin)
    index=find(strcmpi(varargin{k},fields),1);
    if isempty(index), error('Unknown option: %s',varargin{k}); end
    cfg.(fields{index})=varargin{k+1};
end
end
