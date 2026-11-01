function app = exploreIDSDenseMaxGUI(varargin)
% exploreIDSDenseMaxGUI  Interactive dense-direction max-readout display.
%
% This is an educational, unfitted IDS model.  A homogeneous bank of MT-like
% direction mechanisms spans -180:179 degrees.  Each RDK generates a mean
% population response from a drift pedestal, two probe pedestals, the change
% step (change RDK only), and uniformly distributed random dots.  Smooth,
% correlated Gaussian mechanism noise is added before taking the maximum.
%
% For a noise stream at direction theta, the displayed gain is
%
%   E[d max_phi R(phi) / d x(theta)]
%       = sum_phi P(phi wins) h(phi-theta),
%
% where h is the common MT direction-tuning function.  Thus there is no
% separately added "win gain."  The no-change curve is plotted with a minus
% sign because increasing its maximum favors the incorrect RDK.
%
% Population responses are expressed relative to the response to 100%
% preferred-direction coherence.  Consequently, the response to 0% coherence
% has a meaningful fixed level rather than being changed by panel-wise peak
% normalization.  The population response itself is deliberately isolated in
% populationResponse().  That function can be replaced by the response/
% normalization equation from an earlier MT model without changing the GUI or
% max-readout calculation.
%
% Example:
%   app = exploreIDSDenseMaxGUI;
%
% Name-value options:
%   'ProbeOffsetDeg'       default 120
%   'StepCoherencePC'      default 12.3
%   'PedestalCoherencePC'  default 7.071
%   'MTSigmaDeg'           default 37.5
%   'MechanismNoiseSDPC'   default 6.5
%   'NoiseCorrelationDeg'  default 30
%   'NSimulations'         default 4000
%   'RandomSeed'           default 7391

p = inputParser;
p.addParameter('ProbeOffsetDeg',120,@isPositiveOffset);
p.addParameter('StepCoherencePC',12.3,@isNonnegativeScalar);
p.addParameter('PedestalCoherencePC',7.071,@isNonnegativeScalar);
p.addParameter('MTSigmaDeg',37.5,@isPositiveScalar);
p.addParameter('MechanismNoiseSDPC',6.5,@isNonnegativeScalar);
p.addParameter('NoiseCorrelationDeg',30,@isPositiveScalar);
p.addParameter('NSimulations',4000, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 100);
p.addParameter('RandomSeed',7391, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.parse(varargin{:});
opts = p.Results;

dirs = (-180:179)';
plotDirs = (-180:180)';
nDirections = numel(dirs);
nSimulations = round(opts.NSimulations);
pedestalPC = double(opts.PedestalCoherencePC);
sigmaMT = double(opts.MTSigmaDeg);

% H(i,j) is the sensitivity of the mechanism centered at dirs(i) to a
% coherence perturbation centered at dirs(j).
directionDifference = circularDifference(dirs,dirs');
H = exp(-0.5*(directionDifference/sigmaMT).^2);

stream = RandStream('mt19937ar','Seed',round(opts.RandomSeed));
noiseChange = correlatedDirectionNoise(stream,nDirections,nSimulations, ...
  opts.NoiseCorrelationDeg);
noiseNoChange = correlatedDirectionNoise(stream,nDirections,nSimulations, ...
  opts.NoiseCorrelationDeg);
exampleIndex = 1;
accumulation = emptyAccumulation();

fig = uifigure('Name','IDS dense-direction max-readout explorer', ...
  'Position',[40 30 1500 1030],'Color',[0.97 0.97 0.97]);
outer = uigridlayout(fig,[2 1]);
outer.RowHeight = {122,'1x'};
outer.ColumnWidth = {'1x'};
outer.Padding = [10 10 10 10];
outer.RowSpacing = 8;

controlPanel = uipanel(outer,'Title','Controls and fixed assumptions');
controlPanel.Layout.Row = 1;
controlGrid = uigridlayout(controlPanel,[2 10]);
controlGrid.RowHeight = {48,28};
controlGrid.ColumnWidth = {92,'1x',48,105,78,124,78,118,110,'1x'};
controlGrid.Padding = [8 6 8 5];
controlGrid.ColumnSpacing = 7;

uilabel(controlGrid,'Text','Probe offset','HorizontalAlignment','right');
offsetSlider = uislider(controlGrid,'Limits',[0 180], ...
  'MajorTicks',0:30:180,'MinorTicks',[],'Value',opts.ProbeOffsetDeg);
offsetSlider.Layout.Column = 2;
offsetValueLabel = uilabel(controlGrid,'Text',sprintf('%3.0f%c', ...
  opts.ProbeOffsetDeg,char(176)),'HorizontalAlignment','left');
offsetValueLabel.Layout.Column = 3;

uilabel(controlGrid,'Text','Step coherence','HorizontalAlignment','right');
stepField = uieditfield(controlGrid,'numeric','Limits',[0 70], ...
  'Value',opts.StepCoherencePC,'ValueDisplayFormat','%.2f %%');
stepField.Layout.Column = 5;

uilabel(controlGrid,'Text','Mechanism noise','HorizontalAlignment','right');
noiseField = uieditfield(controlGrid,'numeric','Limits',[0 40], ...
  'Value',opts.MechanismNoiseSDPC,'ValueDisplayFormat','%.2f PC');
noiseField.Layout.Column = 7;

uilabel(controlGrid,'Text','Noise corr. width','HorizontalAlignment','right');
corrField = uieditfield(controlGrid,'numeric','Limits',[1 180], ...
  'Value',opts.NoiseCorrelationDeg,'ValueDisplayFormat','%.1f deg');
corrField.Layout.Column = 9;

resampleButton = uibutton(controlGrid,'push','Text','New example trial');
resampleButton.Layout.Column = 10;

assumptionText = uilabel(controlGrid, ...
  'Text',sprintf(['Homogeneous 1%c direction bank; MT SD %.1f%c; pedestal %.3f%%; ' ...
  '%d fixed Monte Carlo noise fields.  Both columns use identical y-axis scales.'], ...
  char(176),sigmaMT,char(176),pedestalPC,nSimulations), ...
  'FontAngle','italic','FontColor',[0.25 0.25 0.25]);
assumptionText.Layout.Row = 2;
assumptionText.Layout.Column = [1 10];

panelGrid = uigridlayout(outer,[1 2]);
panelGrid.Layout.Row = 2;
panelGrid.ColumnWidth = {'1x','1x'};
panelGrid.ColumnSpacing = 10;
panelGrid.Padding = [0 0 0 0];

[changePanel,changeAxes] = makeRDKPanel(panelGrid,'Change RDK',[0.12 0.43 0.72]);
changePanel.Layout.Column = 1;
[noChangePanel,noChangeAxes] = makeRDKPanel(panelGrid,'No-change RDK',[0.76 0.25 0.20]);
noChangePanel.Layout.Column = 2;

offsetSlider.ValueChangingFcn = @(~,event) updateFromSlider(event.Value);
offsetSlider.ValueChangedFcn = @(~,event) updateFromSlider(event.Value);
stepField.ValueChangedFcn = @(~,~) resetAccumulation();
noiseField.ValueChangedFcn = @(~,~) resetAccumulation();
corrField.ValueChangedFcn = @(~,~) regenerateNoise();
resampleButton.ButtonPushedFcn = @(~,~) resampleExample();

updateDisplay(true);

app = struct('Figure',fig,'ProbeOffsetSlider',offsetSlider, ...
  'StepCoherenceField',stepField,'MechanismNoiseField',noiseField, ...
  'NoiseCorrelationField',corrField,'ChangeAxes',changeAxes, ...
  'NoChangeAxes',noChangeAxes);

  function updateFromSlider(value)
    roundedValue = round(value);
    offsetSlider.Value = roundedValue;
    offsetValueLabel.Text = sprintf('%3.0f%c',roundedValue,char(176));
    updateDisplay(true);
  end

  function resetAccumulation()
    accumulation = emptyAccumulation();
    updateDisplay(true);
  end

  function regenerateNoise()
    fig.Pointer = 'watch';
    drawnow;
    noiseChange = correlatedDirectionNoise(stream,nDirections,nSimulations, ...
      corrField.Value);
    noiseNoChange = correlatedDirectionNoise(stream,nDirections,nSimulations, ...
      corrField.Value);
    exampleIndex = 1;
    accumulation = emptyAccumulation();
    fig.Pointer = 'arrow';
    updateDisplay(true);
  end

  function resampleExample()
    exampleIndex = randi(stream,nSimulations);
    updateDisplay(false);
  end

  function updateDisplay(accumulateCurrentPoint)
    if ~isvalid(fig)
      return
    end
    offset = round(offsetSlider.Value);
    stepPC = stepField.Value;
    mechanismNoiseSD = noiseField.Value;

    change = calculateRDK(true,offset,stepPC,mechanismNoiseSD,noiseChange);
    noChange = calculateRDK(false,offset,stepPC,mechanismNoiseSD,noiseNoChange);

    if accumulateCurrentPoint
      accumulation = storeAccumulationPoint(accumulation,offset,change,noChange);
    end

    changeIsOverallWinner = change.exampleMaximum >= noChange.exampleMaximum;

    sharedLimits = commonYLimits(change,noChange);
    plotOneRDK(changeAxes,change,true,offset,sharedLimits, ...
      changeIsOverallWinner,accumulation);
    plotOneRDK(noChangeAxes,noChange,false,offset,sharedLimits, ...
      ~changeIsOverallWinner,accumulation);
    drawnow limitrate;
  end

  function result = calculateRDK(isChange,offset,stepPC,noiseSD,noiseBank)
    if isChange
      thisStep = stepPC;
    else
      thisStep = 0;
    end
    componentDirections = [0; offset; -offset];
    componentCoherences = [pedestalPC+thisStep; pedestalPC; pedestalPC];
    randomPC = max(0,100-sum(componentCoherences));

    meanResponse = populationResponse(componentDirections,componentCoherences, ...
      randomPC,dirs,H);
    noisyResponses = meanResponse + noiseSD*noiseBank;
    [~,winnerIndex] = max(noisyResponses,[],1);
    winnerProbability = accumarray(winnerIndex(:),1,[nDirections 1]) / ...
      nSimulations;

    % Probability-weighted sensitivity of the winning mechanism.  This is
    % the continuous-bank counterpart of a kernel/regression gain.
    sensoryGain = (winnerProbability' * H)';
    if isChange
      decisionGain = sensoryGain;
    else
      decisionGain = -sensoryGain;
    end

    % Express response in units of the peak response produced by 100%
    % preferred-direction coherence.  Unlike panel-wise peak normalization,
    % this preserves the response level for 0% coherence and permits direct
    % comparison of change and no-change responses.
    responseScalePC = 100;
    normalizedMean = meanResponse/responseScalePC;
    exampleResponse = meanResponse + noiseSD*noiseBank(:,exampleIndex);
    exampleResponse = exampleResponse/responseScalePC;
    [exampleMaximum,exampleWinner] = max(exampleResponse);
    zeroCoherenceResponse = mean(H(:,1));

    plusProbeIndex = directionIndex(dirs,offset);
    minusProbeIndex = directionIndex(dirs,-offset);
    plusProbeGain = decisionGain(plusProbeIndex);
    minusProbeGain = decisionGain(minusProbeIndex);
    meanProbeGain = mean([plusProbeGain minusProbeGain]);

    result = struct('componentDirections',componentDirections, ...
      'componentCoherences',componentCoherences,'stepPC',thisStep, ...
      'randomPC',randomPC,'meanResponse',normalizedMean, ...
      'exampleResponse',exampleResponse,'exampleWinner',exampleWinner, ...
      'exampleMaximum',exampleMaximum, ...
      'zeroCoherenceResponse',zeroCoherenceResponse, ...
      'winnerProbability',winnerProbability,'decisionGain',decisionGain, ...
      'plusProbeGain',plusProbeGain,'minusProbeGain',minusProbeGain, ...
      'meanProbeGain',meanProbeGain);
  end

  function limits = commonYLimits(change,noChange)
    coherenceMaximum = max([change.componentCoherences(:); ...
      noChange.componentCoherences(:); change.randomPC/360; ...
      noChange.randomPC/360]);
    responseMaximum = max([change.meanResponse(:);change.exampleResponse(:); ...
      noChange.meanResponse(:);noChange.exampleResponse(:); ...
      change.zeroCoherenceResponse;noChange.zeroCoherenceResponse]);
    winnerMaximum = max([change.winnerProbability(:); ...
      noChange.winnerProbability(:)]);

    limits = struct();
    limits.coherence = [0 max(1,1.08*coherenceMaximum)];
    responseMargin = max(0.02,0.05*responseMaximum);
    limits.response = [0 responseMaximum+responseMargin];
    limits.winner = [0 max(0.005,1.08*winnerMaximum)];
    % Winner sensitivities lie between zero and one before the no-change
    % decision sign is applied.  A fixed symmetric scale makes the two sides
    % immediately comparable as controls are changed.
    limits.gain = [-1.05 1.05];
  end

  function plotOneRDK(ax,R,isChange,offset,limits,isOverallWinner,A)
    panelColor = [0.12 0.43 0.72];
    if ~isChange
      panelColor = [0.76 0.25 0.20];
    end

    % Coherence components.
    cla(ax(1)); hold(ax(1),'on');
    randomPerDegree = R.randomPC/360;
    plot(ax(1),plotDirs,repmat(randomPerDegree,size(plotDirs)),':', ...
      'Color',[0.45 0.45 0.45],'LineWidth',1.2,'DisplayName','Random dots / degree');
    hDrift = stem(ax(1),0,pedestalPC,'filled','Color',[0.15 0.35 0.65], ...
      'LineWidth',1.5,'MarkerSize',5,'DisplayName','Drift pedestal');
    if isChange && R.stepPC > 0
      hStep = stem(ax(1),0,pedestalPC+R.stepPC,'Color',[0.05 0.05 0.05], ...
        'LineWidth',2,'Marker','none','DisplayName','Drift + step');
    else
      hStep = gobjects(0);
    end
    hProbe = stem(ax(1),[offset -offset],[pedestalPC pedestalPC],'filled', ...
      'Color',[0.85 0.48 0.10],'LineWidth',1.5,'MarkerSize',5, ...
      'DisplayName','Probe pedestals');
    hold(ax(1),'off');
    ylabel(ax(1),'Coherence (%)');
    title(ax(1),sprintf('Stimulus components; random dots %.1f%%',R.randomPC));
    ylim(ax(1),limits.coherence);
    if isempty(hStep)
      legend(ax(1),[hDrift hProbe],{'Drift pedestal','Probe pedestals'}, ...
        'Location','northeast','FontSize',8);
    else
      legend(ax(1),[hDrift hStep hProbe], ...
        {'Drift pedestal','Drift + step','Probe pedestals'}, ...
        'Location','northeast','FontSize',8);
    end

    % Population response and one correlated-noise example.
    cla(ax(2)); hold(ax(2),'on');
    plot(ax(2),plotDirs,closeCurve(R.exampleResponse),'Color',[0.72 0.72 0.72], ...
      'LineWidth',1,'DisplayName','Example noisy response');
    plot(ax(2),plotDirs,closeCurve(R.meanResponse),'Color',panelColor, ...
      'LineWidth',2,'DisplayName','Mean response');
    yline(ax(2),R.zeroCoherenceResponse,'--', ...
      '0% coherence / population mean','Color',[0.25 0.25 0.25], ...
      'LineWidth',1.1,'LabelHorizontalAlignment','left', ...
      'LabelVerticalAlignment','bottom','HandleVisibility','off');
    winnerDirection = dirs(R.exampleWinner);
    if isOverallWinner
      winnerFill = [0.20 0.72 0.28];
      winnerName = 'Overall example winner';
    else
      winnerFill = [0 0 0];
      winnerName = 'Example patch winner (loser)';
    end
    plot(ax(2),winnerDirection,R.exampleResponse(R.exampleWinner),'kv', ...
      'MarkerFaceColor',winnerFill,'MarkerEdgeColor','k','MarkerSize',7, ...
      'DisplayName',winnerName);
    hold(ax(2),'off');
    ylabel(ax(2),'Normalized response');
    title(ax(2),'Dense MT population response');
    ylim(ax(2),limits.response);
    legend(ax(2),'Location','northeast','FontSize',8);

    % Winner-direction distribution.
    cla(ax(3));
    plot(ax(3),plotDirs,closeCurve(R.winnerProbability), ...
      'Color',panelColor,'LineWidth',1.8);
    ylabel(ax(3),'P(direction wins)');
    title(ax(3),sprintf('Winner distribution (%d simulations)',nSimulations));
    ylim(ax(3),limits.winner);

    % Probability-weighted sensitivity of the winner.
    cla(ax(4)); hold(ax(4),'on');
    yline(ax(4),0,'Color',[0.5 0.5 0.5]);
    plot(ax(4),plotDirs,closeCurve(R.decisionGain), ...
      'Color',panelColor,'LineWidth',1.5,'DisplayName','All stream directions');
    plot(ax(4),[offset -offset],[R.plusProbeGain R.minusProbeGain],'o', ...
      'Color',[0.85 0.48 0.10],'MarkerFaceColor',[1.0 0.72 0.22], ...
      'MarkerSize',7,'LineWidth',1.2,'DisplayName','Actual probes');
    hold(ax(4),'off');
    ylabel(ax(4),'Decision gain');
    ylim(ax(4),limits.gain);
    title(ax(4),sprintf(['Probe gain: +%g%c %.3f, -%g%c %.3f; ' ...
      'mean %.3f'],offset,char(176),R.plusProbeGain,offset,char(176), ...
      R.minusProbeGain,R.meanProbeGain));
    legend(ax(4),'Location','northeast','FontSize',8);

    % Accumulated probe-gain estimates as the offset slider is explored.
    cla(ax(5)); hold(ax(5),'on');
    yline(ax(5),0,'Color',[0.5 0.5 0.5]);
    [sortedOffsets,order] = sort(A.offset);
    if isChange
      accumulatedGain = A.changeMean(order);
    else
      accumulatedGain = A.noChangeMean(order);
    end
    plot(ax(5),sortedOffsets,accumulatedGain,'o-', ...
      'Color',panelColor,'MarkerFaceColor',panelColor,'MarkerSize',5, ...
      'LineWidth',1.5);
    hold(ax(5),'off');
    ylabel(ax(5),'Mean probe gain');
    xlabel(ax(5),'Probe offset (deg)');
    xlim(ax(5),[0 180]);
    ax(5).XTick = 0:30:180;
    if isChange
      ylim(ax(5),[-0.1 1]);
      ax(5).YTick = [-0.1 0 0.2 0.4 0.6 0.8 1];
    else
      ylim(ax(5),[-1 0.1]);
      ax(5).YTick = [-1 -0.8 -0.6 -0.4 -0.2 0 0.1];
    end
    title(ax(5),'Accumulated probe-gain assessment');

    for a = 1:4
      xlim(ax(a),[-180 180]);
      ax(a).XTick = -180:60:180;
      ax(a).Box = 'on';
      ax(a).FontSize = 9;
      grid(ax(a),'on');
      ax(a).XTickLabel = {};
    end
    ax(5).Box = 'on';
    ax(5).FontSize = 9;
    grid(ax(5),'on');
  end
end

function A = emptyAccumulation()
A = struct('offset',zeros(0,1),'changeMean',zeros(0,1), ...
  'noChangeMean',zeros(0,1));
end

function A = storeAccumulationPoint(A,offset,change,noChange)
index = find(A.offset == offset,1);
if isempty(index)
  index = numel(A.offset)+1;
  A.offset(index,1) = offset;
end
A.changeMean(index,1) = change.meanProbeGain;
A.noChangeMean(index,1) = noChange.meanProbeGain;
end

function index = directionIndex(dirs,direction)
[~,index] = min(abs(circularDifference(dirs,direction)));
end

function [panel,ax] = makeRDKPanel(parent,titleText,titleColor)
panel = uipanel(parent,'Title',titleText,'FontSize',14,'FontWeight','bold', ...
  'ForegroundColor',titleColor);
gridLayout = uigridlayout(panel,[5 1]);
gridLayout.RowHeight = {'1x','1x','1x','1x','1x'};
gridLayout.ColumnWidth = {'1x'};
gridLayout.RowSpacing = 5;
gridLayout.Padding = [7 5 7 7];
ax = gobjects(5,1);
for k = 1:5
  ax(k) = uiaxes(gridLayout);
  ax(k).Layout.Row = k;
end
end

function response = populationResponse(componentDirections,componentCoherences, ...
  randomPC,dirs,H)
% Replace this function when the exact earlier MT normalization is restored.
% For now, coherent components sum through the common MT tuning function.
% Uniformly random dots contribute a direction-independent response floor.
response = zeros(size(dirs));
for k = 1:numel(componentDirections)
  [~,directionIndex] = min(abs(circularDifference(dirs,componentDirections(k))));
  response = response + componentCoherences(k)*H(:,directionIndex);
end
response = response + randomPC*mean(H(:,1));
end

function E = correlatedDirectionNoise(stream,nDirections,nSamples,correlationDeg)
% Generate a stationary circular Gaussian field with approximately the
% requested correlation SD.  Nearby mechanisms therefore do not become an
% infinite set of independent extreme-value competitors as grid spacing falls.
whiteNoise = randn(stream,nDirections,nSamples);
lag = (0:nDirections-1)';
lag = min(lag,nDirections-lag);
filterSD = max(correlationDeg/sqrt(2),0.5);
kernel = exp(-0.5*(lag/filterSD).^2);
kernel = kernel/norm(kernel);
E = real(ifft(fft(whiteNoise,[],1).*fft(kernel),[],1));
E = E/std(E(:));
end

function y = closeCurve(x)
y = [x(:); x(1)];
end

function d = circularDifference(a,b)
d = mod(double(a)-double(b)+180,360)-180;
end

function tf = isPositiveOffset(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 180;
end

function tf = isPositiveScalar(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x) && x > 0;
end

function tf = isNonnegativeScalar(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0;
end
