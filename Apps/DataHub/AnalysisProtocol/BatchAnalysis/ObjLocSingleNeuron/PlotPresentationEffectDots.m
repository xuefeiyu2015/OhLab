function fig = PlotPresentationEffectDots(E, figNum, names)
% PlotPresentationEffectDots  Slide figure: point estimates with 95% CI, one
% row per measure, FEF/SC side by side, labeled reference lines that say what
% each value would mean (e.g. "independent", "identical tuning").
%
% Inputs
%   E — struct with
%       .Est, .Lo, .Hi   R-by-2 (rows = measures, cols = areas)
%       .RowLabels       cellstr{R}
%       .Refs            struct array (.Value, .Label)
%       .XLabel, .XLim, .LogX (logical), .Headline, .Subtitle
%       .pArea           R-by-1 Holm p for area difference (NaN = none)
%       .LegendPos       [x y w h] normalized position of the FEF/SC key
%   figNum, names
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-05

st  = PresentationStyle();
fig = SlideFigure(figNum, E.Headline, st);
ax  = axes(fig, 'Position', [0.24 0.15 0.50 0.58]); hold(ax, 'on');
R   = size(E.Est, 1);
off = [0.12 -0.12];             % FEF above SC, matching the legend order

for k = 1:numel(E.Refs)
    xline(ax, E.Refs(k).Value, '--', 'Color', st.Ref, 'LineWidth', 1.2);
    text(ax, E.Refs(k).Value, R + 0.42, E.Refs(k).Label, 'Color', st.Muted, ...
        'FontSize', st.TickSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
for r = 1:R
    y = R - r + 1;
    for j = 1:2
        yy = y + off(j);
        plot(ax, [E.Lo(r,j) E.Hi(r,j)], [yy yy], '-', 'Color', st.Area{j}, 'LineWidth', 2.5);
        plot(ax, E.Est(r,j), yy, 'o', 'MarkerSize', 10, 'MarkerFaceColor', st.Area{j}, ...
            'MarkerEdgeColor', st.Surface, 'LineWidth', 1.5);
        text(ax, E.Hi(r,j), yy, sprintf('  %.2f', E.Est(r,j)), 'Color', st.Ink, ...
            'FontSize', st.LabelSize, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    end
    if ~isnan(E.pArea(r))
        text(ax, E.XLim(2), y, ['  ' SlidePString(E.pArea(r), names)], 'Color', st.Muted, ...
            'FontSize', st.TickSize, 'VerticalAlignment', 'middle');
    end
end
SlideStyleAxes(ax, st);
if E.LogX
    tk = [0.25 0.5 1 2 4 8 16]; tk = tk(tk >= E.XLim(1) & tk <= E.XLim(2));
    set(ax, 'XScale', 'log', 'XTick', tk, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), tk, 'UniformOutput', false));
end
set(ax, 'YTick', 1:R, 'YTickLabel', flip(E.RowLabels), 'YLim', [0.5 R + 0.55], 'XLim', E.XLim, ...
    'TickLabelInterpreter', 'tex');
ax.YAxis.FontSize = st.LabelSize; ax.YAxis.Color = st.Ink;
xlabel(ax, E.XLabel, 'FontSize', st.LabelSize, 'Color', st.Ink);
SlideAreaLegend(fig, st, names, E.LegendPos);
SlideSubtitle(fig, E.Subtitle, st);
end
