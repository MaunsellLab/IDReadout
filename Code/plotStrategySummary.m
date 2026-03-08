function plotStrategySummary(S, basePath)
% plotStrategySummary  Make a compact visual summary of strategy-model outputs.
%
% Usage:
%   plotStrategySummary(S)
%   plotStrategySummary(S, basePath)
%
% Input:
%   S         struct returned by summarizeStrategyModels()
%   basePath  project root (default folderPath())
%
% Output files:
%   <basePath>/Summaries/strategy_summary.pdf
%
% Panels:
%   1) Session performance (P(correct), nTrials)
%   2) Session deviance explained for m1..m4
%   3) Session incremental deviance improvement over m1
%   4) Preferred-kernel trial forecast (50/125/250 ms)
%   5) Beta1 (m1) with SE
%   6) Text summary / model definitions

  % ---- Handle inputs ----
  if nargin < 1 || isempty(S)
    error('plotStrategySummary:MissingInput', ...
        'Input S (from summarizeStrategyModels) is required.');
  end

  if nargin < 2 || isempty(basePath)
    basePath = folderPath();
  end
  if isstring(basePath)
    basePath = char(basePath);
  end

  summaryFolder = fullfile(basePath, 'Summaries');
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

  % Shorten labels to MM/DD when possible
  for i = 1:numel(sessLabels)
    s = erase(sessLabels(i), ".mat");
    tok = regexp(s, '\d{8}', 'match', 'once');
    if ~isempty(tok)
      sessLabels(i) = extractBetween(tok, 5, 6) + "/" + extractBetween(tok, 7, 8);
    else
      s = erase(s, "IDReadout_");
      sessLabels(i) = s;
    end
  end

  % ---- Figure ----
  figure(301); clf;
  set(gcf, 'Color', 'w', 'Name', 'Behavioral Strategy Summary');
  set(gcf, 'Position', [80 60 1500 900]);

  tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

  % ============================================================
  % Panel 1: session performance
  % ============================================================
  ax1 = nexttile;
  hold(ax1, 'on');

  yyaxis(ax1, 'left')
  if ismember('pCorrect', T.Properties.VariableNames)
    plot(ax1, x, T.pCorrect, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    ylabel(ax1, 'P(correct)');
    ylim(ax1, [0 1]);
  end

  if isfield(S, 'pooled') && isfield(S.pooled, 'summary') && isfield(S.pooled.summary, 'pCorrect')
    yline(ax1, S.pooled.summary.pCorrect, ':', 'Pooled', ...
        'LabelHorizontalAlignment', 'left');
  end

  yyaxis(ax1, 'right')
  if ismember('nTrials', T.Properties.VariableNames)
    plot(ax1, x, T.nTrials, '--s', 'LineWidth', 1.2, 'MarkerSize', 4);
    ylabel(ax1, 'nTrials');
    ymax = max(T.nTrials, [], 'omitnan');
    if isfinite(ymax) && ymax > 0
      ylim(ax1, [0 1.05 * ymax]);
    else
      ylim(ax1, [0 1]);
    end
  end

  xlim(ax1, [0.5 nSessions + 0.5]);
  grid(ax1, 'on');
  title(ax1, 'Session performance');
  set(ax1, 'XTick', x, 'XTickLabel', sessLabels, 'XTickLabelRotation', 45);

  % ============================================================
  % Panel 2: session deviance explained
  % ============================================================
  ax2 = nexttile;
  hold(ax2, 'on');

  legendEntries = {};

  if ismember('devExpl_m1', T.Properties.VariableNames)
    plot(ax2, x, T.devExpl_m1, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    legendEntries{end+1} = 'm1';
  end
  if ismember('devExpl_m2', T.Properties.VariableNames)
    plot(ax2, x, T.devExpl_m2, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    legendEntries{end+1} = 'm2';
  end
  if ismember('devExpl_m3', T.Properties.VariableNames)
    plot(ax2, x, T.devExpl_m3, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    legendEntries{end+1} = 'm3';
  end
  if ismember('devExpl_m4', T.Properties.VariableNames)
    plot(ax2, x, T.devExpl_m4, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    legendEntries{end+1} = 'm4';
  end

  if isfield(S, 'pooled') && isfield(S.pooled, 'summary')
    Ps = S.pooled.summary;
    if isfield(Ps, 'devExplained_p1')
      yline(ax2, Ps.devExplained_p1, ':', 'pooled m1', ...
          'LabelHorizontalAlignment', 'left');
    end
    if isfield(Ps, 'devExplained_p2')
      yline(ax2, Ps.devExplained_p2, ':', 'pooled m2', ...
          'LabelHorizontalAlignment', 'left');
    end
  end

  if ~isempty(legendEntries)
    legend(ax2, legendEntries, 'Location', 'best');
  end

  xlim(ax2, [0.5 nSessions + 0.5]);
  ylabel(ax2, 'Deviance explained');
  title(ax2, 'Session model performance');
  grid(ax2, 'on');
  set(ax2, 'XTick', x, 'XTickLabel', sessLabels, 'XTickLabelRotation', 45);

  % ============================================================
  % Panel 3: increment over m1
  % ============================================================
  ax3 = nexttile;
  hold(ax3, 'on');

  legendEntries = {};

  if ismember('incOver_m1_m2', T.Properties.VariableNames)
    plot(ax3, x, T.incOver_m1_m2, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    legendEntries{end+1} = 'm2 - m1';
  end
  if ismember('incOver_m1_m3', T.Properties.VariableNames)
    plot(ax3, x, T.incOver_m1_m3, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    legendEntries{end+1} = 'm3 - m1';
  end
  if ismember('incOver_m1_m4', T.Properties.VariableNames)
    plot(ax3, x, T.incOver_m1_m4, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    legendEntries{end+1} = 'm4 - m1';
  end

  yline(ax3, 0, '--');

  if ~isempty(legendEntries)
    legend(ax3, legendEntries, 'Location', 'best');
  end

  xlim(ax3, [0.5 nSessions + 0.5]);
  ylabel(ax3, 'ΔDeviance over m1');
  title(ax3, 'Increment beyond net evidence');
  grid(ax3, 'on');
  set(ax3, 'XTick', x, 'XTickLabel', sessLabels, 'XTickLabelRotation', 45);

  % ============================================================
  % Panel 4: preferred forecast
  % ============================================================
  ax4 = nexttile;
  hold(ax4, 'on');

  prefLegend = {};
  plottedPref = false;
  displayCeil = 1e9;   % clip only for display

  if ismember('Npred_pref_50', T.Properties.VariableNames)
    y = T.Npred_pref_50;
    y(~isfinite(y) | y <= 0) = NaN;
    y(y > displayCeil) = displayCeil;
    plot(ax4, x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    prefLegend{end+1} = '50 ms';
    plottedPref = true;
  end
  if ismember('Npred_pref_125', T.Properties.VariableNames)
    y = T.Npred_pref_125;
    y(~isfinite(y) | y <= 0) = NaN;
    y(y > displayCeil) = displayCeil;
    plot(ax4, x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    prefLegend{end+1} = '125 ms';
    plottedPref = true;
  end
  if ismember('Npred_pref_250', T.Properties.VariableNames)
    y = T.Npred_pref_250;
    y(~isfinite(y) | y <= 0) = NaN;
    y(y > displayCeil) = displayCeil;
    plot(ax4, x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    prefLegend{end+1} = '250 ms';
    plottedPref = true;
  end

  if plottedPref
    set(ax4, 'YScale', 'log');
    ylim(ax4, [1e3 displayCeil]);
    yticks(ax4, [1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
    yline(ax4, 1e4, ':', '10^4', 'LabelHorizontalAlignment', 'left');
    yline(ax4, 1e5, ':', '10^5', 'LabelHorizontalAlignment', 'left');
    yline(ax4, 1e6, ':', '10^6', 'LabelHorizontalAlignment', 'left');
    yline(ax4, 1e9, ':', '10^9+', 'LabelHorizontalAlignment', 'left');
    legend(ax4, prefLegend, 'Location', 'best');
  end

  xlim(ax4, [0.5 nSessions + 0.5]);
  ylabel(ax4, 'Predicted trials');
  title(ax4, 'Preferred-kernel forecast (+2 SEM)');
  grid(ax4, 'on');
  set(ax4, 'XTick', x, 'XTickLabel', sessLabels, 'XTickLabelRotation', 45);

  % ============================================================
  % Panel 5: beta1 with SE
  % ============================================================
  ax5 = nexttile;
  hold(ax5, 'on');

  if ismember('b_m1_sumDM', T.Properties.VariableNames)
    y = T.b_m1_sumDM;

    if ismember('se_m1_sumDM', T.Properties.VariableNames)
      e = T.se_m1_sumDM;
      errorbar(ax5, x, y, e, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    else
      plot(ax5, x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    end
  end

  yline(ax5, 0, '--');

  xlim(ax5, [0.5 nSessions + 0.5]);
  ylabel(ax5, '\beta_1 (m1)');
  title(ax5, '\beta_1 choice sensitivity');
  grid(ax5, 'on');
  set(ax5, 'XTick', x, 'XTickLabel', sessLabels, 'XTickLabelRotation', 45);

  % ============================================================
  % Panel 6: text summary / model definitions
  % ============================================================
  ax6 = nexttile;
  cla(ax6);
  axis(ax6, 'off');
  title(ax6, 'Summary');

  lines = strings(0,1);

  if isfield(S, 'info')
    if isfield(S.info, 'nSessionsLoaded')
      lines(end+1) = sprintf('Sessions pooled: %d', S.info.nSessionsLoaded);
    elseif isfield(S.info, 'nFilesLoaded')
      lines(end+1) = sprintf('Files loaded: %d', S.info.nFilesLoaded);
    end

    if isfield(S.info, 'nTrialsTotal')
      lines(end+1) = sprintf('Trials pooled: %d', S.info.nTrialsTotal);
    end
  end

  if isfield(S, 'pooled') && isfield(S.pooled, 'summary')
    Ps = S.pooled.summary;

    if isfield(Ps, 'pCorrect')
      lines(end+1) = sprintf('P(correct): %.3f', Ps.pCorrect);
    end

    lines(end+1) = "";

    if isfield(Ps, 'devExplained_p1')
      lines(end+1) = sprintf('Pooled m1 devExpl: %.6f', Ps.devExplained_p1);
    end
    if isfield(Ps, 'devExplained_p2')
      lines(end+1) = sprintf('Pooled m2 devExpl: %.6f', Ps.devExplained_p2);
    end
  end

  % median forecasts
  lines(end+1) = "";
  lines(end+1) = "Median preferred forecast N (+2 SEM):";
  lines(end+1) = sprintf('  50 ms  : %s', localFmtMedian(T, 'Npred_pref_50'));
  lines(end+1) = sprintf('  125 ms : %s', localFmtMedian(T, 'Npred_pref_125'));
  lines(end+1) = sprintf('  250 ms : %s', localFmtMedian(T, 'Npred_pref_250'));

  if ismember('b_m1_sumDM', T.Properties.VariableNames)
    lines(end+1) = "";
    lines(end+1) = sprintf('Median β_1: %s', localFmtMedian(T, 'b_m1_sumDM'));
  end

  lines(end+1) = "";
  lines(end+1) = "Forecast values > 10^9 clipped in plot";

  lines(end+1) = "";
  lines(end+1) = "Model definitions:";
  lines(end+1) = "m1: net evidence (sumDM)";
  lines(end+1) = "m2: net evidence + run structure";
  lines(end+1) = "m3: net evidence + early/late weighting";
  lines(end+1) = "m4: net evidence + cumulative excursions";

  text(ax6, 0.02, 0.98, strjoin(lines, newline), ...
      'Units', 'normalized', ...
      'VerticalAlignment', 'top', ...
      'HorizontalAlignment', 'left', ...
      'FontName', 'Courier', ...
      'FontSize', 8);

  sgtitle(tl, 'Behavioral Strategy Summary');

  % ---- Save ----
  pdfFile = fullfile(summaryFolder, 'strategy_summary.pdf');
  exportgraphics(gcf, pdfFile, 'ContentType', 'vector');
  fprintf('  Saved strategy summary: %s\n', pdfFile);
end


function s = localFmtMedian(T, varName)
% helper for panel 6 text
  s = 'n/a';
  if ~ismember(varName, T.Properties.VariableNames)
    return;
  end

  x = T.(varName);
  x = x(isfinite(x));

  if startsWith(varName, 'Npred_')
    x = x(x > 0);
  end

  if isempty(x)
    return;
  end

  m = median(x, 'omitnan');

  if abs(m) >= 1e5
    s = sprintf('%.2g', m);
  elseif abs(m) >= 1000
    s = sprintf('%d', round(m));
  elseif abs(m) < 1e-3 && m ~= 0
    s = sprintf('%.2e', m);
  else
    s = sprintf('%.4g', m);
  end
end