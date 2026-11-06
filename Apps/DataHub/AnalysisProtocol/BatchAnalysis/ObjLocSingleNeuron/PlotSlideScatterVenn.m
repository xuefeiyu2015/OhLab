function fig = PlotSlideScatterVenn(aucObj, pObj, aucLoc, pLoc, O, alpha, figNum, headline)
% PlotSlideScatterVenn  Slide: object vs location AUC scatter with marginal
% histograms (significant high/low-preferring neurons filled) and a Venn of
% object / both / location coding counts, with the chance expectation.
%
% Inputs
%   aucObj, pObj, aucLoc, pLoc — n-vectors (one area)
%   O        — ComputeOverlapVsChance output for the same neurons
%   alpha    — significance level; figNum; headline (char)
% Output: fig — figure handle
%
% Colors: object red / location purple / both orange; histogram low-preferring
% object blue, low-preferring location green (each pair CVD-validated).
%
% Xuefei Yu, 2026-11-05

st = PresentationStyle(); st.FigPos = [100 60 900 800];
hex = @(h) sscanf(h(2:end), '%2x%2x%2x', [1 3]) / 255;
cObj = hex('#d6402b'); cLoc = hex('#8a5cc2'); cBoth = hex('#c98500');
cObjLow = hex('#2a78d6'); cLocLow = hex('#199e70'); cGray = st.Ref;

fig = SlideFigure(figNum, headline, st);
so = pObj < alpha; sl = pLoc < alpha;
edges = 0:0.05:1;

% --- main scatter
ax = axes(fig, 'Position', [0.09 0.08 0.50 0.52]); hold(ax, 'on');
nn = ~so & ~sl;
scatter(ax, aucObj(nn), aucLoc(nn), 26, cGray, 'LineWidth', 0.8);
scatter(ax, aucObj(so & ~sl), aucLoc(so & ~sl), 34, cObj,  'filled', 'MarkerEdgeColor', st.Surface);
scatter(ax, aucObj(~so & sl), aucLoc(~so & sl), 34, cLoc,  'filled', 'MarkerEdgeColor', st.Surface);
scatter(ax, aucObj(so & sl),  aucLoc(so & sl),  40, cBoth, 'filled', 'MarkerEdgeColor', st.Surface);
xline(ax, 0.5, ':', 'Color', st.Ink, 'LineWidth', 1.5); yline(ax, 0.5, ':', 'Color', st.Ink, 'LineWidth', 1.5);
SlideStyleAxes(ax, st);
set(ax, 'XLim', [0 1], 'YLim', [0 1], 'XTick', 0:0.5:1, 'YTick', 0:0.5:1);
xlabel(ax, 'AUC (Object value)', 'FontSize', st.LabelSize, 'Color', st.Ink);
ylabel(ax, 'AUC (Location value)', 'FontSize', st.LabelSize, 'Color', st.Ink);
text(ax, 0.04, 0.96, sprintf('N = %d', O.N), 'FontSize', st.LabelSize, 'Color', st.Ink, 'VerticalAlignment', 'top');
text(ax, 0.96, 0.04, sprintf('r = %.2f', O.Rho), 'FontSize', st.LabelSize, 'Color', st.Ink, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'BackgroundColor', st.Surface, 'Margin', 1);

% --- top histogram: object AUC
axT = axes(fig, 'Position', [0.09 0.62 0.50 0.14]); hold(axT, 'on');
histogram(axT, aucObj, edges, 'FaceColor', 'none', 'EdgeColor', cGray);
histogram(axT, aucObj(so & aucObj > 0.5), edges, 'FaceColor', cObj,    'EdgeColor', 'none', 'FaceAlpha', 1);
histogram(axT, aucObj(so & aucObj < 0.5), edges, 'FaceColor', cObjLow, 'EdgeColor', 'none', 'FaceAlpha', 1);
SlideStyleAxes(axT, st); set(axT, 'XLim', [0 1], 'XTick', [], 'XColor', 'none');
ylabel(axT, 'Neurons', 'FontSize', st.TickSize, 'Color', st.Muted);

% --- right histogram: location AUC
axR = axes(fig, 'Position', [0.61 0.08 0.14 0.52]); hold(axR, 'on');
histogram(axR, aucLoc, edges, 'Orientation', 'horizontal', 'FaceColor', 'none', 'EdgeColor', cGray);
histogram(axR, aucLoc(sl & aucLoc > 0.5), edges, 'Orientation', 'horizontal', 'FaceColor', cLoc,    'EdgeColor', 'none', 'FaceAlpha', 1);
histogram(axR, aucLoc(sl & aucLoc < 0.5), edges, 'Orientation', 'horizontal', 'FaceColor', cLocLow, 'EdgeColor', 'none', 'FaceAlpha', 1);
SlideStyleAxes(axR, st); set(axR, 'YLim', [0 1], 'YTick', [], 'YColor', 'none');
xlabel(axR, 'Neurons', 'FontSize', st.TickSize, 'Color', st.Muted);

% --- Venn
axV = axes(fig, 'Position', [0.68 0.60 0.28 0.22]); hold(axV, 'on');
th = linspace(0, 2*pi, 201); th(end) = [];   % no duplicate closing vertex
pa = polyshape(-0.55 + cos(th), sin(th)); pb = polyshape(0.55 + cos(th), sin(th));
plot(axV, pa, 'FaceColor', cObj, 'FaceAlpha', 0.85, 'EdgeColor', 'none');
plot(axV, pb, 'FaceColor', cLoc, 'FaceAlpha', 0.85, 'EdgeColor', 'none');
plot(axV, intersect(pa, pb), 'FaceColor', cBoth, 'FaceAlpha', 1, 'EdgeColor', 'none');
vt = {'FontSize', st.LabelSize + 2, 'FontWeight', 'bold', 'Color', [1 1 1], 'HorizontalAlignment', 'center'};
text(axV, -0.95, 0, sprintf('%d', O.Obs(1)), vt{:});
text(axV,  0,    0, sprintf('%d', O.Obs(2)), vt{:});
text(axV,  0.95, 0, sprintf('%d', O.Obs(3)), vt{:});
text(axV, -0.55, 1.25, 'Object', 'FontSize', st.TickSize, 'Color', st.Ink, 'HorizontalAlignment', 'center');
text(axV,  0.55, 1.25, 'Location', 'FontSize', st.TickSize, 'Color', st.Ink, 'HorizontalAlignment', 'center');
text(axV, 0, -1.35, sprintf('both expected by chance: %.0f', O.Exp(2)), 'FontSize', st.TickSize, ...
    'Color', st.Muted, 'HorizontalAlignment', 'center');
axis(axV, 'equal', 'off'); set(axV, 'XLim', [-1.7 1.7], 'YLim', [-1.6 1.5]);
end
