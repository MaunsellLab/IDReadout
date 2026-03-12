function plotStrategySummary(S, path)
% plotStrategySummary  Make a compact visual summary of strategy-model outputs.
%
% Usage:
%   plotStrategySummary(S)
%   plotStrategySummary(S, path)
%
% Input:
%   S      struct returned by summarizeStrategyModels()
%   path   project root (default folderPath())
%
% Output files:
%   <path>/Summaries/strategy_summary.pdf
%   <path>/Summaries/strategy_summary.fig
%
% Panels:
%   1) Session performance (P(correct), nTrials)
%   2) Session deviance explained for m1..m4
%   3) Session incremental deviance improvement over m1
%   4) Pooled model deviance explained

  % ---- Handle inputs ----
  if nargin < 1 || isempty(S)
    error('plotStrategySummary:MissingInput', ...
        'Input S (from summarizeStrategyModels) is required.');
  end

  if nargin < 2 || isempty(path)
    path = folderPath();
  end
  if isstring(path)
    path = char(path);
  end

  summaryFolder = fullfile(path, 'Summaries');
  if ~exist(summaryFolder, 'dir')
    mkdir(summaryFolder);
  end

  if ~isfield(S, 'sessionTable') || ~istable(S.sessionTable) || isempty(S.sessionTable)
    warning('plotStrategySummary:NoSessionTable', ...
        'No sessionTable found in S. Nothing to plot.');
    return;
  end

  T = S.sessionTable;

  % ---- Sort by session if possible ----
  if ismember('session', T.Properties.VariableNames)
    T = sortrows(T, 'session');
  end

  nSessions = height(T);
  x = 1:nSessions;

  % ---- Labels ----
  if ismember('session', T.Properties.VariableNames)
    sessLabels = string(T.session);
  elseif ismember('file', T.Properties.VariableNames)
    sessLabels = string(T.file);
  else
    sessLabels = "S" + string(x);
  end

  % optionally shorten labels a bit
  % for i = 1:numel(sessLabels)
  %   s = sessLabels(i);
  %   s = erase(s, ".mat");
  %   s = erase(s, "IDReadout_");
  %   sessLabels(i) = s;
  % end
  % ---- Labels ----
  if ismember('session', T.Properties.VariableNames)
      sessLabels = string(T.session);
  elseif ismember('file', T.Properties.VariableNames)
      sessLabels = string(T.file);
  else
      sessLabels = "S" + string(x);
  end
  
  % extract date-like token if present
  for i = 1:numel(sessLabels)
      s = erase(sessLabels(i), ".mat");
      tok = regexp(s, '\d{8}', 'match', 'once');
      if ~isempty(tok)
          % convert YYYYMMDD -> MM/DD
          sessLabels(i) = extractBetween(tok, 5, 6) + "/" + extractBetween(tok, 7, 8);
      else
          sessLabels(i) = erase(s, "IDReadout_");
      end
  end

  % ---- Figure ----
  figure(301); clf;
  set(gcf, 'Color', 'w', 'Name', 'Strategy Summary');
  set(gcf, 'Position', [100 100 1200 850]);

  tl = tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

  % ============================================================
  % Panel 1: session performance
  % ============================================================
  nexttile;
  hold on;

  yyaxis left
  if ismember('pCorrect', T.Properties.VariableNames)
    plot(x, T.pCorrect, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    ylabel('P(correct)');
    ylim([0 1]);
  end
  % if isfield(S, 'pooled') && isfield(S.pooled, 'summary') && isfield(S.pooled.summary, 'pCorrect')
  %     yyaxis left
  %     yline(S.pooled.summary.pCorrect, ':', 'Pooled', 'LabelHorizontalAlignment', 'left');
  % end

  if ismember('failed', T.Properties.VariableNames)
    bad = find(T.failed);
    if ~isempty(bad)
        yyaxis left
        plot(bad, T.pCorrect(bad), 'x', 'MarkerSize', 10, 'LineWidth', 1.5);
    end
  end

  yyaxis right
  if ismember('nTrials', T.Properties.VariableNames)
    plot(x, T.nTrials, '--s', 'LineWidth', 1.2, 'MarkerSize', 4);
    ylabel('nTrials');
    ylim([0 inf]);
  end

  xlim([0.5 nSessions + 0.5]);
  grid on;
  title('Session performance');
  set(gca, 'XTick', x, 'XTickLabel', sessLabels, 'XTickLabelRotation', 45);

  % ============================================================
  % Panel 2: session deviance explained
  % ============================================================
  ax2 = nexttile;
  hold on;

  plottedAny = false;

  if ismember('devExpl_m1', T.Properties.VariableNames)
    plot(x, T.devExpl_m1, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    plottedAny = true;
  end
  if ismember('devExpl_m2', T.Properties.VariableNames)
    plot(x, T.devExpl_m2, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    plottedAny = true;
  end
  if ismember('devExpl_m3', T.Properties.VariableNames)
    plot(x, T.devExpl_m3, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    plottedAny = true;
  end
  if ismember('devExpl_m4', T.Properties.VariableNames)
    plot(x, T.devExpl_m4, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    plottedAny = true;
  end

  if plottedAny
    legendEntries = {};
    if ismember('devExpl_m1', T.Properties.VariableNames), legendEntries{end+1} = 'm1'; end %#ok<AGROW>
    if ismember('devExpl_m2', T.Properties.VariableNames), legendEntries{end+1} = 'm2'; end %#ok<AGROW>
    if ismember('devExpl_m3', T.Properties.VariableNames), legendEntries{end+1} = 'm3'; end %#ok<AGROW>
    if ismember('devExpl_m4', T.Properties.VariableNames), legendEntries{end+1} = 'm4'; end %#ok<AGROW>
    legend(legendEntries, 'Location', 'best');

    % if isfield(S, 'pooled') && isfield(S.pooled, 'summary')
    %   Ps = S.pooled.summary;
    %   if isfield(Ps, 'devExplained_p1')
    %       yline(Ps.devExplained_p1, ':', 'pooled m1', 'LabelHorizontalAlignment', 'left');
    %   end
    %   if isfield(Ps, 'devExplained_p2')
    %       yline(Ps.devExplained_p2, ':', 'pooled m2', 'LabelHorizontalAlignment', 'left');
    %   end
    % end
  end

  xlim([0.5 nSessions + 0.5]);
  yLimits = ylim();
  ylim([0, yLimits(2) * 1.35]);
  ylabel('Deviance explained');
  title('Session model performance');
  grid on;
  set(gca, 'XTick', x, 'XTickLabel', sessLabels, 'XTickLabelRotation', 45);
  
  modelText = {
  'Model definitions:'
  'm1: net evidence (sumDM)'
  'm2: net evidence + run structure'
  'm3: net evidence + early/late weighting'
  'm4: net evidence + cumulative excursions'
  % ''
  % 'Formulas:'
  % 'm1: isCorrect ~ sumDM'
  % 'm2: isCorrect ~ sumDM + longestPosRun + longestNegRun + nSwitch'
  % 'm3: isCorrect ~ sumDM + dEarlyLate'
  % 'm4: isCorrect ~ sumDM + maxCumDM + minCumDM'
  };
  ax = gca;
  pos = ax.Position;
   % [pos(1)+pos(3)*0.55, pos(2)+pos(4)*0.05, pos(3)*0.42, pos(4)*0.3], ...
  text(0.025, 0.975, modelText, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'units', 'normalized', 'fontSize', 9);
    % 'String', modelText, ...
    % 'FitBoxToText', 'on', ...
    % 'BackgroundColor', 'white', ...
    % 'EdgeColor', [0.7 0.7 0.7], ...
    % 'FontSize', 9);
  % annotation(gcf,'textbox', ...
  %  [pos(1), 0, pos(3)*0.42, pos(4)*0.3], ...
  %   'String', modelText, ...
  %   'FitBoxToText', 'on', ...
  %   'BackgroundColor', 'white', ...
  %   'EdgeColor', [0.7 0.7 0.7], ...
  %   'FontSize', 9);
  % 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...


  % ============================================================
  % Panel 3: increment over m1
  % ============================================================
  nexttile;
  hold on;

  plottedAny = false;

  if ismember('incOver_m1_m2', T.Properties.VariableNames)
    plot(x, T.incOver_m1_m2, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    plottedAny = true;
  end
  if ismember('incOver_m1_m3', T.Properties.VariableNames)
    plot(x, T.incOver_m1_m3, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    plottedAny = true;
  end
  if ismember('incOver_m1_m4', T.Properties.VariableNames)
    plot(x, T.incOver_m1_m4, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    plottedAny = true;
  end

  yline(0, '--');

  if plottedAny
    legendEntries = {};
    if ismember('incOver_m1_m2', T.Properties.VariableNames), legendEntries{end+1} = 'm2 - m1'; end %#ok<AGROW>
    if ismember('incOver_m1_m3', T.Properties.VariableNames), legendEntries{end+1} = 'm3 - m1'; end %#ok<AGROW>
    if ismember('incOver_m1_m4', T.Properties.VariableNames), legendEntries{end+1} = 'm4 - m1'; end %#ok<AGROW>
    legend(legendEntries, 'Location', 'best');
  end

  xlim([0.5 nSessions + 0.5]);
  ylabel('\DeltaDeviance over m1');
  title('Increment beyond net evidence');
  grid on;
  set(gca, 'XTick', x, 'XTickLabel', sessLabels, 'XTickLabelRotation', 45);

  % ============================================================
  % Panel 4: pooled summary
  % ============================================================
  nexttile;
  cla;
  axis off;
  title('Pooled summary');

  lines = strings(0,1);

  if isfield(S, 'info')
    if isfield(S.info, 'nSessionsLoaded')
      lines(end+1) = sprintf('Sessions pooled: %d', S.info.nSessionsLoaded); %#ok<AGROW>
    elseif isfield(S.info, 'nFilesLoaded')
      lines(end+1) = sprintf('Files loaded: %d', S.info.nFilesLoaded); %#ok<AGROW>
    end

    if isfield(S.info, 'nTrialsTotal')
      lines(end+1) = sprintf('Trials pooled: %d', S.info.nTrialsTotal); %#ok<AGROW>
    end
  end

  if isfield(S, 'pooled') && isfield(S.pooled, 'summary')
    Ps = S.pooled.summary;

    if isfield(Ps, 'pCorrect')
      lines(end+1) = sprintf('P(correct): %.3f', Ps.pCorrect); %#ok<AGROW>
    end

    lines(end+1) = ""; %#ok<AGROW>
    lines(end+1) = "Pooled deviance explained:"; %#ok<AGROW>

    if isfield(Ps, 'devExplained_p1')
      lines(end+1) = sprintf('  p1: %.6f', Ps.devExplained_p1); %#ok<AGROW>
    end
    if isfield(Ps, 'devExplained_p2')
      lines(end+1) = sprintf('  p2: %.6f', Ps.devExplained_p2); %#ok<AGROW>
    end

    if isfield(S.pooled, 'compare')
      C = S.pooled.compare;
      lines(end+1) = ""; %#ok<AGROW>
      lines(end+1) = "Nested comparisons:"; %#ok<AGROW>

      if isfield(C, 'p1_vs_p0')
        lines(end+1) = sprintf('  p1 vs p0: dDev %.4f, p %.4g', ...
            C.p1_vs_p0.deltaDeviance, C.p1_vs_p0.pValue); %#ok<AGROW>
      end
      if isfield(C, 'p2_vs_p1')
        lines(end+1) = sprintf('  p2 vs p1: dDev %.4f, p %.4g', ...
            C.p2_vs_p1.deltaDeviance, C.p2_vs_p1.pValue); %#ok<AGROW>
      end
    end

    if isfield(Ps, 'sumDM_diffCorrectError')
      lines(end+1) = ""; %#ok<AGROW>
      lines(end+1) = sprintf('Δ mean(sumDM)         : %.4f', Ps.sumDM_diffCorrectError);
    end
    if isfield(Ps, 'longestNegRun_diffCorrectError')
      lines(end+1) = sprintf('Δ mean(longestNegRun) : %.4f', Ps.longestNegRun_diffCorrectError);
    end
    if isfield(Ps, 'nSwitch_diffCorrectError')
      lines(end+1) = sprintf('Δ mean(nSwitch)       : %.4f', Ps.nSwitch_diffCorrectError);
    end
  end

  if isempty(lines)
    lines = "No pooled summary available.";
  end

  text(0.02, 0.98, strjoin(lines, newline), ...
      'Units', 'normalized', ...
      'VerticalAlignment', 'top', ...
      'HorizontalAlignment', 'left', ...
      'FontName', 'Courier', ...
      'FontSize', 11);


  % ---- Save ----
  pdfFile = fullfile(summaryFolder, 'strategy_summary.pdf');
  exportgraphics(gcf, pdfFile, 'ContentType','vector');
  fprintf('  Saved strategy summary: %s\n', pdfFile);
end