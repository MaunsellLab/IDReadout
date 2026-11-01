function noiseProfile=runMaxPoolingIndependentNoiseProfile( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGINDEPENDENTNOISEPROFILE Profile structured mechanism noise.
p=inputParser;
p.addParameter('IndependentFractionGrid',[0 .05 .1 .2 .4 .7 1],@isnumeric);
p.addParameter('VarianceExponentGrid',0,@isnumeric);
p.addParameter('NDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('MinimumGainSessions',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end
cache=maxpoolfit.loadCache(cachePath); result=loadInput(resultInput,'result');
assert(isfield(result,'finalFit')&&isfield(result.finalFit,'adaptationByGroup'), ...
  'A completed design-group fit is required.');
noiseProfile=maxpoolfit.profileIndependentNoiseMixture(cache,result.finalFit, ...
  'IndependentFractionGrid',p.Results.IndependentFractionGrid, ...
  'VarianceExponentGrid',p.Results.VarianceExponentGrid, ...
  'NDraws',p.Results.NDraws, ...
  'MinimumGainSessions',p.Results.MinimumGainSessions);
noiseProfile.cachePath=string(cachePath);
save(fullfile(outputFolder,'MaxPoolingIndependentNoiseProfile.mat'), ...
  'noiseProfile','-v7.3');
writetable(noiseProfile.overallProfile,fullfile(outputFolder, ...
  'MaxPoolingIndependentNoiseProfile_Overall.csv'));
writetable(noiseProfile.groupProfile,fullfile(outputFolder, ...
  'MaxPoolingIndependentNoiseProfile_Groups.csv'));
writetable(noiseProfile.gainProfile,fullfile(outputFolder, ...
  'MaxPoolingIndependentNoiseProfile_Gains.csv'));
writetable(noiseProfile.winnerProfile,fullfile(outputFolder, ...
  'MaxPoolingIndependentNoiseProfile_Winners.csv'));
plotProfile(noiseProfile,outputFolder);
end

function value=loadInput(input,name)
if isstruct(input), value=input; return; end
S=load(input,name); assert(isfield(S,name), ...
  'MAT file does not contain variable %s.',name); value=S.(name);
end

function plotProfile(P,outputFolder)
O=P.overallProfile; exponents=unique(O.noiseVarianceExponent,'stable');
colors=lines(numel(exponents));
fig=figure('Color','w','Position',[80 80 1050 760]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
names={'deltaLogLFromBest','absoluteGainRMSE','ratioRMSE','meanNoiseSDPC'};
titles={'Choice likelihood','Absolute changed-side gain error', ...
  'Probe/drift gain-ratio error','Refitted sensory-noise scale'};
ylabels={'Delta log likelihood','Gain RMSE','Ratio RMSE','Noise SD (% coherence)'};
for a=1:4
  ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
  for j=1:numel(exponents)
    q=O.noiseVarianceExponent==exponents(j);
    plot(ax,O.independentNoiseFraction(q),O.(names{a})(q),'o-', ...
      'Color',colors(j,:),'MarkerFaceColor',colors(j,:), ...
      'DisplayName',sprintf('p = %g',exponents(j)),'LineWidth',1.2);
  end
  xlabel(ax,'Independent variance fraction, lambda'); ylabel(ax,ylabels{a});
  title(ax,titles{a});
  if a==1, legend(ax,'Location','best'); end
end
title(tl,sprintf('Independent response-scaled noise; %d draws',P.nNoiseDraws));
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingIndependentNoiseProfile.pdf'),'ContentType','vector'); close(fig);
end
