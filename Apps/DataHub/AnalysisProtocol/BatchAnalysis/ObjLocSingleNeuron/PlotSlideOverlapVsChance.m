function fig = PlotSlideOverlapVsChance(Os, panelTitles, layout, figNum, headline, subtitle)
% PlotSlideOverlapVsChance  Slide: the Venn as bars (object only / both /
% location only) with the independence expectation drawn on each bar
% (black line = expected, gray band = 95% chance range).
%
% Inputs
%   Os          — cell of ComputeOverlapVsChance outputs, one per panel
%   panelTitles — cellstr, one per panel
%   layout      — [nRows nCols] subplot grid
%   figNum, headline, subtitle
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-05

st = PresentationStyle();
hex = @(h) sscanf(h(2:end), '%2x%2x%2x', [1 3]) / 255;
cols = {hex('#d6402b'), hex('#c98500'), hex('#8a5cc2')};   % object only, both, location only
if layout(1) > 1, st.FigPos = [100 60 1000 760]; end
fig = SlideFigure(figNum, headline, st);
SlideSubtitle(fig, subtitle, st);

nP = numel(Os); w = 0.6;
top = 0.74; bottom = 0.08; hgap = 0.10; vgap = 0.12;
pw = (0.86 - hgap*(layout(2)-1)) / layout(2);
ph = (top - bottom - vgap*(layout(1)-1)) / layout(1);
yMax = max(cellfun(@(O) max([O.Obs(1:3), O.NullHi(1:3)]), Os)) * 1.30;
for p = 1:nP
    O = Os{p};
    r = ceil(p / layout(2)); c = p - (r-1)*layout(2);
    ax = axes(fig, 'Position', [0.08 + (c-1)*(pw + hgap), top - r*ph - (r-1)*vgap, pw, ph]); hold(ax, 'on');
    for k = 1:3
        bar(ax, k, O.Obs(k), w, 'FaceColor', cols{k}, 'EdgeColor', 'none');
        patch(ax, k + [-1 1 1 -1]*w/2*1.1, [O.NullLo(k) O.NullLo(k) O.NullHi(k) O.NullHi(k)], ...
            [0.6 0.6 0.6], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
        plot(ax, k + [-1 1]*w/2*1.1, [O.Exp(k) O.Exp(k)], '-', 'Color', st.Ink, 'LineWidth', 2.5);
        text(ax, k, max(O.Obs(k), O.NullHi(k)) + 0.02*yMax, sprintf('%d', O.Obs(k)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', st.LabelSize, 'FontWeight', 'bold', 'Color', st.Ink);
    end
    if O.pFisher < 0.001, ps = 'p < .001'; else, ps = sprintf('p = %.3f', O.pFisher); end
    text(ax, 2, yMax*0.97, sprintf('both: %d observed vs %.0f chance\nOR = %.2f, %s', O.Obs(2), O.Exp(2), O.OR, ps), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', st.TickSize, 'Color', st.Muted);
    SlideStyleAxes(ax, st);
    set(ax, 'XTick', 1:3, 'XTickLabel', O.Labels(1:3), 'XLim', [0.4 3.6], 'YLim', [0 yMax]);
    if c == 1, ylabel(ax, 'Neurons', 'FontSize', st.LabelSize, 'Color', st.Ink); end
    title(ax, sprintf('%s  (n = %d)', panelTitles{p}, O.N), 'FontSize', st.LabelSize, 'Color', st.Ink, 'FontWeight', 'bold');
end
end
