function summary = runSplitSignalDesignSweep(varargin)
% RUNSPLITSIGNALDESIGNSWEEP  Compare split-signal stimulus designs and models.
%
% summary=runSplitSignalDesignSweep();
% summary=runSplitSignalDesignSweep('NDesignTrials',5000,'CreateFigures',false);

cfg=splitSignalDefaults(); cfg=parseOptions(cfg,varargin{:});
if cfg.SaveOutputs && ~exist(cfg.OutputDir,'dir'), mkdir(cfg.OutputDir); end

bank=makeSplitSignalTrialBank(cfg.NDesignTrials,cfg,cfg.Seed);
rows=[]; row=0;
for separation=cfg.SeparationsDeg
    for noiseAmplitude=cfg.StepNoiseHalfAmplitudes
        for prestepFraction=cfg.PrestepPedestalFractions
            for balance=cfg.IncludeBalancedPrestep
                if balance && prestepFraction==0, continue; end
                for adaptation=cfg.AdaptationFractions
                    for totalStep=cfg.TotalStepCoherences
                        requiredRise=2*(1-prestepFraction)*noiseAmplitude;
                        if totalStep+1e-9<requiredRise, continue; end
                        for allocation=cfg.AllocationFractions
                            row=row+1;
                            nextResult=simulateSplitSignalDesign(bank,separation, ...
                                allocation,totalStep,noiseAmplitude,prestepFraction, ...
                                balance,adaptation,cfg);
                            if row==1
                                rows=nextResult;
                            else
                                rows(row)=nextResult; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
        end
    end
end

summary=struct(); summary.config=cfg; summary.conditions=rows;
summary.primary=selectPrimary(rows,cfg);
for s=1:numel(cfg.SeparationsDeg)
    separation=cfg.SeparationsDeg(s);
    field=sprintf('sep%d',round(separation));
    summary.joint.(field)=fitSplitSignalJointModels(bank,separation,cfg);
    names={'maxMinusLinear','maxMinusVector','maxMinusEnergy'};
    for k=1:numel(names)
        name=names{k};
        summary.power.(field).(name)=estimateSplitSignalSessions( ...
            summary.joint.(field).pointDeltaLL.(name),cfg,cfg.Seed+1000+10*s+k);
    end
end

if cfg.CreateFigures, summary.figures=plotSplitSignalDesignSweep(summary,cfg); end
if cfg.SaveOutputs
    save(fullfile(cfg.OutputDir,'SplitSignalDesignSummary.mat'),'summary','-v7.3');
    writeSummaryCSV(rows,fullfile(cfg.OutputDir,'SplitSignalDesignSummary.csv'));
    if cfg.CreateFigures
        saveas(summary.figures.performance,fullfile(cfg.OutputDir,'SplitSignalPerformance.pdf'));
        saveas(summary.figures.gains,fullfile(cfg.OutputDir,'SplitSignalGains.pdf'));
        saveas(summary.figures.power,fullfile(cfg.OutputDir,'SplitSignalPower.pdf'));
    end
end
end

function primary=selectPrimary(rows,cfg)
keep=arrayfun(@(r) r.noiseHalfAmplitude==cfg.PrimaryStepNoiseHalfAmplitude && ...
    r.prestepFraction==cfg.PrimaryPrestepPedestalFraction && ...
    r.balancePrestep==cfg.PrimaryBalancePrestep && ...
    abs(r.adaptationFraction-cfg.PrimaryAdaptationFraction)<1e-10,rows);
primary=rows(keep);
end

function writeSummaryCSV(rows,path)
fid=fopen(path,'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,['separation_deg,allocation,total_step,noise_half_amplitude,' ...
    'prestep_fraction,balance,adaptation_fraction,performance,' ...
    'common_pedestal_rise,allocatable_signal,' ...
    'winner_gain,loser_gain,loser_winner_ratio,' ...
    'delta_ll_linear,delta_ll_vector,delta_ll_energy\n']);
for r=rows
    fprintf(fid,'%.6g,%.6g,%.6g,%.6g,%.6g,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n', ...
        r.separationDeg,r.allocation,r.totalStep,r.noiseHalfAmplitude, ...
        r.prestepFraction,r.balancePrestep,r.adaptationFraction,r.performance, ...
        r.commonPedestalRise,r.allocatableSignal, ...
        r.winnerGain,r.loserGain,r.loserWinnerGainRatio, ...
        r.deltaLL.maxMinusLinear,r.deltaLL.maxMinusVector,r.deltaLL.maxMinusEnergy);
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
