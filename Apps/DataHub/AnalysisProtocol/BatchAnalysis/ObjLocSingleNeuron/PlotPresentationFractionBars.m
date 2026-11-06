function fig = PlotPresentationFractionBars(F, figNum, names, yLabel)
% PlotPresentationFractionBars  Slide figure: % neurons per category, grouped
% by area, direct value labels, chance line, area-difference p above groups.
%
% Inputs
%   F       — struct with .Frac (2 cat x 2 area), .k, .n, .pArea (1 x 2 cat),
%             .Categories, .Chance, .Headline, .Subtitle
%   figNum  — figure number; names — {'FEF','SC'}; yLabel — y-axis label
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-05

st  = PresentationStyle();
fig = SlideFigure(figNum, F.Headline, st);
ax  = axes(fig, 'Position', [0.10 0.13 0.62 0.62]); hold(ax, 'on');

nC = size(F.Frac, 1); w = 0.34; gap = 0.02;
for j = 1:2
    x = (1:nC) + (j - 1.5) * (w + gap);
    bar(ax, x, 100*F.Frac(:,j), w, 'FaceColor', st.Area{j}, 'EdgeColor', 'none');
    for i = 1:nC
        text(ax, x(i), 100*F.Frac(i,j) + 1.2, sprintf('%.0f%%', 100*F.Frac(i,j)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', st.LabelSize, 'Color', st.Ink, 'FontWeight', 'bold');
        if 100*F.Frac(i,j) >= 8      % count inside the bar only when it fits
            text(ax, x(i), 0.8, sprintf('%d/%d', F.k(i,j), F.n(i,j)), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', st.TickSize - 2, 'Color', [1 1 1]);
        end
    end
end
yTop = 100 * max(F.Frac(:)) * 1.35 + 5;
for i = 1:nC
    yl = 100 * max(F.Frac(i,:)) + 7;
    plot(ax, i + [-1 1]*(w+gap)/2, [yl yl], '-', 'Color', st.Muted, 'LineWidth', 1);
    text(ax, i, yl + 0.8, SlidePString(F.pArea(i), names), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', st.TickSize, 'Color', st.Muted);
end
yline(ax, 100*F.Chance, '--', 'Color', st.Ref, 'LineWidth', 1.2, 'HandleVisibility', 'off');
text(ax, nC + 0.55, 100*F.Chance, 'chance', 'Color', st.Muted, 'FontSize', st.TickSize, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
SlideStyleAxes(ax, st);
set(ax, 'XTick', 1:nC, 'XTickLabel', F.Categories, 'XLim', [0.4 nC + 0.6], 'YLim', [0 yTop]);
ylabel(ax, yLabel, 'FontSize', st.LabelSize, 'Color', st.Ink);
SlideAreaLegend(fig, st, names, [0.76 0.55 0.2 0.2]);
SlideSubtitle(fig, F.Subtitle, st);
end
