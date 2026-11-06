function Fig = PlotValueSaccadeBehavior(Behavior, FigName)
%PLOTVALUESACCADEBEHAVIOR Plot high vs low value saccade metrics for one
%session. Rendering only; takes the output of ComputeValueSaccadeBehavior.
%
% Row 1: per-metric trial scatter + mean +/- sem (RT, peak velocity,
%        gaze duration in target window)
% Row 2: cumulative distributions of the same metrics
%
% Inputs:
%   Behavior : struct from ComputeValueSaccadeBehavior
%   FigName  : figure name string
% Output:
%   Fig      : figure handle
%
% Xuefei Yu

% Lab colour convention: object red/blue, location purple/green
if isfield(Behavior, 'TaskCode')
    [ColorGood, ColorBad] = ValueColors(Behavior.TaskCode);
else
    [ColorGood, ColorBad] = ValueColors('Object');
end
Labels  = {'Reaction time (ms)', 'Peak velocity (deg/s)', 'Gaze duration (ms)'};
Metrics = Behavior.MetricNames;

Fig = PrepareFigure(111, 'w', [50, 100, 1200, 700], 'Name', FigName);

for m = 1:numel(Metrics)
    Val = Behavior.(Metrics{m});
    G = Val(Behavior.IsGood);  G = G(~isnan(G));
    B = Val(~Behavior.IsGood); B = B(~isnan(B));
    S = Behavior.Summary.(Metrics{m});

    % --- Trial scatter + mean/sem ---
    ax = subplot(2, numel(Metrics), m); hold(ax, 'on');
    Rs = RandStream('mt19937ar', 'Seed', 0); % deterministic jitter, no global RNG change
    scatter(ax, 1 + 0.15*randn(Rs, size(G)), G, 10, ColorGood, 'filled', 'MarkerFaceAlpha', 0.3);
    scatter(ax, 2 + 0.15*randn(Rs, size(B)), B, 10, ColorBad,  'filled', 'MarkerFaceAlpha', 0.3);
    errorbar(ax, 1, S.GoodMean, S.GoodSem, 's', 'Color', 'k', 'MarkerFaceColor', ColorGood, 'MarkerSize', 9, 'LineWidth', 1.5);
    errorbar(ax, 2, S.BadMean,  S.BadSem,  's', 'Color', 'k', 'MarkerFaceColor', ColorBad,  'MarkerSize', 9, 'LineWidth', 1.5);
    set(ax, 'XTick', [1 2], 'XTickLabel', {sprintf('High (n=%d)', S.GoodN), sprintf('Low (n=%d)', S.BadN)}, ...
        'XLim', [0.4 2.6], 'TickDir', 'out');
    ylabel(ax, Labels{m});
    title(ax, sprintf('AUC=%.2f, p=%.2g', S.AUC, S.P));
    box(ax, 'off');

    % --- Cumulative distribution ---
    ax2 = subplot(2, numel(Metrics), numel(Metrics) + m); hold(ax2, 'on');
    if ~isempty(G)
        [f, x] = ecdf(G); stairs(ax2, x, f, 'Color', ColorGood, 'LineWidth', 1.5);
    end
    if ~isempty(B)
        [f, x] = ecdf(B); stairs(ax2, x, f, 'Color', ColorBad, 'LineWidth', 1.5);
    end
    xlabel(ax2, Labels{m}); ylabel(ax2, 'Cumulative fraction');
    set(ax2, 'TickDir', 'out'); box(ax2, 'off');
    if m == 1
        legend(ax2, {'High value', 'Low value'}, 'Location', 'southeast', 'Box', 'off');
    end
end

annotation(Fig, 'textbox', [0.45 0.955 0.54 0.04], 'String', ...
    sprintf('Gaze window r=%g deg | excluded: artifact %d, landed outside %d, never left %d', ...
    Behavior.Params.TargetWindowRadius, sum(Behavior.Artifact), ...
    sum(Behavior.LandedOutside), sum(Behavior.GazeCensored)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'FontSize', 9);

end
