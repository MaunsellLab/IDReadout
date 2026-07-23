function summary = analyzeMTGaussianReadoutWidths(context,varargin)
% analyzeMTGaussianReadoutWidths  Signed Gaussian-readout baseline sweep.
%
% This diagnostic intentionally excludes candidate rectification and max
% pooling. It documents the preferred:null ratio and fixed-readout offset
% curve implied by MT tuning plus a nonnegative Gaussian readout.

context=validateAnalysisContext(context);
if context.Mode~="synthetic"
  error('IDReadout:SyntheticContextRequired', ...
    'The MT simulation diagnostic requires a synthetic context.');
end
validateSyntheticManifest(context);

p=inputParser;
addParameter(p,'SigmaMTDeg',37.5,@isPositiveScalar);
addParameter(p,'SigmaReadoutDeg',5:0.5:45, ...
  @(x) isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0));
addParameter(p,'OffsetsDeg',0:180, ...
  @(x) isnumeric(x)&&isvector(x)&&all(isfinite(x)));
addParameter(p,'DisplaySigmaReadoutDeg',[10 20 30], ...
  @(x) isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0));
addParameter(p,'Visible','off',@(x) any(strcmpi(string(x),["on","off"])));
addParameter(p,'SaveOutputs',true,@(x) islogical(x)&&isscalar(x));
parse(p,varargin{:});
opts=p.Results;

mtModel=makeMTReadoutForwardModel('sigmaMTDeg',opts.SigmaMTDeg);
offsets=double(opts.OffsetsDeg(:)');
templates=mtPopulationTemplate(offsets,mtModel);
sigmaR=double(opts.SigmaReadoutDeg(:));
nSigma=numel(sigmaR); nOffsets=numel(offsets);
preferredSensitivity=nan(nSigma,1);
nullSensitivity=nan(nSigma,1);
preferredNullRatio=nan(nSigma,1);
normalizedSensitivity=nan(nSigma,nOffsets);

nullTemplate=mtPopulationTemplate(180,mtModel);
preferredTemplate=mtPopulationTemplate(0,mtModel);
for k=1:nSigma
  bank=makeGaussianReadoutBank(0,sigmaR(k),mtModel);
  a=bank.weightsPhi;
  preferredSensitivity(k)=a*preferredTemplate';
  nullSensitivity(k)=a*nullTemplate';
  preferredNullRatio(k)=preferredSensitivity(k)/abs(nullSensitivity(k));
  normalizedSensitivity(k,:)=(a*templates')/preferredSensitivity(k);
end

widthTable=table(sigmaR,preferredSensitivity,nullSensitivity,preferredNullRatio, ...
  'VariableNames',{'sigmaReadoutDeg','preferredSensitivity', ...
  'nullSensitivity','preferredNullRatio'});
summary=struct();
summary.createdAt=datetime('now');
summary.createdBy=mfilename;
summary.dataOrigin="synthetic";
summary.context=context;
summary.sigmaMTDeg=opts.SigmaMTDeg;
summary.phiDeg=mtModel.phiDeg;
summary.offsetsDeg=offsets;
summary.widthTable=widthTable;
summary.normalizedSensitivity=normalizedSensitivity;
summary.modelDefinition= ...
  "Signed mean-subtracted MT templates; nonnegative Gaussian readout; " + ...
  "fixed preferred-centered readout; no rectification; no pooling";

fig=makePlot(summary,opts.DisplaySigmaReadoutDeg,opts.Visible);
if opts.SaveOutputs
  dataFolder=analysisPath(context,'Common Code','Data','Simulation');
  plotFolder=analysisPath(context,'Common Code','Plots','Simulation');
  if ~isfolder(dataFolder), mkdir(dataFolder); end
  if ~isfolder(plotFolder), mkdir(plotFolder); end
  matPath=fullfile(dataFolder,'MTGaussianReadoutWidthSweep.mat');
  csvPath=fullfile(dataFolder,'MTGaussianReadoutWidthSweep.csv');
  pdfPath=fullfile(plotFolder,'MTGaussianReadoutWidthSweep.pdf');
  summary.outputPaths=string({matPath;csvPath;pdfPath});
  info=whos('summary'); assertSimulationRunCapacity(context,info.bytes);
  save(matPath,'summary','-v7.3');
  writetable(widthTable,csvPath);
  exportgraphics(fig,pdfPath,'ContentType','vector');
  assertSimulationRunCapacity(context,0);
end
if strcmpi(opts.Visible,'off'), close(fig); end
end

function fig=makePlot(summary,displaySigma,visible)
fig=figure('Color','w','Visible',visible,'Units','inches', ...
  'Position',[1 1 10 4.5]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
title(tl,'Gaussian MT Readout: Signed, Unrectified Baseline');

ax=nexttile(tl); hold(ax,'on');
T=summary.widthTable;
plot(ax,T.sigmaReadoutDeg,T.preferredNullRatio,'k-','LineWidth',1.5);
yline(ax,2,'k:'); yline(ax,3,'k:');
xlabel(ax,'Readout \sigma_R (deg)'); ylabel(ax,'Preferred : |null| sensitivity');
grid(ax,'on'); box(ax,'off');

ax=nexttile(tl); hold(ax,'on');
colors=lines(numel(displaySigma)); labels=cell(numel(displaySigma),1);
for k=1:numel(displaySigma)
  [~,row]=min(abs(T.sigmaReadoutDeg-displaySigma(k)));
  plot(ax,summary.offsetsDeg,summary.normalizedSensitivity(row,:), ...
    'LineWidth',1.3,'Color',colors(k,:));
  labels{k}=sprintf('\\sigma_R %.1f deg; P:N %.2f', ...
    T.sigmaReadoutDeg(row),T.preferredNullRatio(row));
end
yline(ax,0,'k:','HandleVisibility','off');
xlabel(ax,'Direction offset (deg)'); ylabel(ax,'Normalized fixed-readout sensitivity');
legend(ax,labels,'Location','best','FontSize',7); grid(ax,'on'); box(ax,'off');
end

function tf=isPositiveScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0;
end
