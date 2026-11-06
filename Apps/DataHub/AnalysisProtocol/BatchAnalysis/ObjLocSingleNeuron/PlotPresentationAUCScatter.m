function fig = PlotPresentationAUCScatter(F, figNum, names)
% PlotPresentationAUCScatter  Slide figure: object vs location AUC, one panel
% per area. Dark = significant for both, mid = one, light = neither. Panel
% text gives observed and noise-corrected r.
%
% Inputs
%   F — struct with .AUCObj{2}, .AUCLoc{2}, .Both{2}, .Either{2}, .rObs, .rDis,
%       .Headline, .Subtitle
%   figNum, names
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-05

st  = PresentationStyle();
fig = SlideFigure(figNum, F.Headline, st);
pos = {[0.08 0.12 0.36 0.62], [0.56 0.12 0.36 0.62]};
for j = 1:2
    ax = axes(fig, 'Position', pos{j}); hold(ax, 'on');
    c  = st.Area{j};
    x = F.AUCObj{j}; y = F.AUCLoc{j};
    none = ~F.Both{j} & ~F.Either{j};
    scatter(ax, x(none), y(none), 28, c, 'filled', 'MarkerFaceAlpha', 0.18, 'MarkerEdgeColor', 'none');
    scatter(ax, x(F.Either{j}), y(F.Either{j}), 34, c, 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none');
    scatter(ax, x(F.Both{j}), y(F.Both{j}), 44, c, 'filled', 'MarkerFaceAlpha', 1, ...
        'MarkerEdgeColor', st.Surface, 'LineWidth', 1);
    xline(ax, 0.5, '-', 'Color', st.Ref); yline(ax, 0.5, '-', 'Color', st.Ref);
    plot(ax, [0 1], [0 1], ':', 'Color', st.Ref, 'LineWidth', 1);
    SlideStyleAxes(ax, st);
    axis(ax, 'square'); set(ax, 'XLim', [0.1 0.9], 'YLim', [0.1 0.9], 'XTick', 0.1:0.2:0.9, 'YTick', 0.1:0.2:0.9);
    xlabel(ax, 'Object value AUC', 'FontSize', st.LabelSize, 'Color', st.Ink);
    if j == 1, ylabel(ax, 'Location value AUC', 'FontSize', st.LabelSize, 'Color', st.Ink); end
    title(ax, sprintf('%s  (n = %d)', names{j}, numel(x)), 'FontSize', st.LabelSize + 2, 'Color', st.Ink);
    text(ax, 0.13, 0.87, sprintf('r = %.2f\nnoise-corrected r = %.2f', F.rObs(j), F.rDis(j)), ...
        'FontSize', st.LabelSize, 'Color', st.Ink, 'VerticalAlignment', 'top', ...
        'BackgroundColor', st.Surface, 'Margin', 2);
    text(ax, 0.88, 0.12, sprintf('dark: both sig (%d)\nmid: one sig (%d)', sum(F.Both{j}), sum(F.Either{j})), ...
        'FontSize', st.TickSize - 1, 'Color', st.Muted, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end
SlideSubtitle(fig, F.Subtitle, st);
end
