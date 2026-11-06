function fig = PlotLatencySlide(Cols, names, figNum, headline, subtitle)
% PlotLatencySlide  Slide: FEF vs SC timing, one column per signal.
% Top row (if Pop given): population time course (mean +/- SEM) per area with
% half-max onset (dotted) and its 95% bootstrap CI (bar). Bottom row:
% cumulative distribution of single-neuron latencies, medians marked.
%
% Inputs
%   Cols — struct array, one per column:
%          .Title, .t (1-by-T), .PopA/.PopB (ComputePopulationOnset or []),
%          .YLabel, .XLimCurve, .LatA/.LatB (vectors, NaN ignored), .LatLabel,
%          .XLimLat, .pPop, .pLat (Holm p, NaN to omit)
%   names — {'FEF','SC'}; figNum; headline; subtitle
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-06

st = PresentationStyle();
nC = numel(Cols);
hasPop = ~isempty(Cols(1).PopA);
st.FigPos = [100 60 520*nC + 120, 560 + 140*hasPop];
fig = SlideFigure(figNum, headline, st);
SlideSubtitle(fig, subtitle, st);
top = 0.74 - 0.06*~hasPop; bottom = 0.10 + 0.02*~hasPop; hgap = 0.08; vgap = 0.13;
nR = 1 + hasPop;
pw = (0.80 - hgap*(nC-1)) / nC;
ph = (top - bottom - vgap*(nR-1)) / nR;

for c = 1:nC
    C = Cols(c); x0 = 0.08 + (c-1)*(pw + hgap);
    if hasPop
        ax = axes(fig, 'Position', [x0, top - ph, pw, ph]); hold(ax, 'on');
        Ps = {C.PopA, C.PopB};
        yl = [min(cellfun(@(P) min(P.Mean - P.SEM), Ps)), max(cellfun(@(P) max(P.Mean + P.SEM), Ps))];
        yl = yl + [-0.05 0.25] * diff(yl);
        for j = 1:2
            P = Ps{j}; m = P.Mean; s = P.SEM; ok = ~isnan(m);
            fill(ax, [C.t(ok) fliplr(C.t(ok))], [m(ok)+s(ok) fliplr(m(ok)-s(ok))], st.Area{j}, ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(ax, C.t, m, '-', 'Color', st.Area{j}, 'LineWidth', 2);
        end
        for j = 1:2
            P = Ps{j}; if isnan(P.Onset), continue; end
            yb = yl(2) - (0.05 + 0.08*(j-1)) * diff(yl);
            xline(ax, P.Onset, ':', 'Color', st.Area{j}, 'LineWidth', 2, 'HandleVisibility', 'off');
            plot(ax, P.CI, [yb yb], '-', 'Color', st.Area{j}, 'LineWidth', 4, 'HandleVisibility', 'off');
            plot(ax, P.Onset, yb, 'o', 'MarkerSize', 8, 'MarkerFaceColor', st.Area{j}, ...
                'MarkerEdgeColor', st.Surface, 'LineWidth', 1.5, 'HandleVisibility', 'off');
            text(ax, C.XLimCurve(2), yb, sprintf('%s %.0f ms [%.0f, %.0f]  n=%d', names{j}, P.Onset, P.CI, P.n), ...
                'FontSize', st.TickSize - 1, 'Color', st.Ink, 'VerticalAlignment', 'middle', ...
                'HorizontalAlignment', 'right');
        end
        xline(ax, 0, '-', 'Color', st.Ref, 'HandleVisibility', 'off');
        yline(ax, 0, '-', 'Color', st.Ref, 'HandleVisibility', 'off');
        SlideStyleAxes(ax, st);
        set(ax, 'XLim', C.XLimCurve, 'YLim', yl);
        ylabel(ax, C.YLabel, 'FontSize', st.LabelSize, 'Color', st.Ink);
        xlabel(ax, 'Time from target onset (ms)', 'FontSize', st.TickSize, 'Color', st.Muted);
        title(ax, {C.Title, PString(C.pPop, 'onset')}, 'FontSize', st.LabelSize, 'Color', st.Ink, 'FontWeight', 'bold');
    end

    ax = axes(fig, 'Position', [x0, bottom, pw, ph]); hold(ax, 'on');
    Ls = {C.LatA, C.LatB};
    for j = 1:2
        v = sort(Ls{j}(~isnan(Ls{j})));
        if isempty(v), continue; end
        y = (1:numel(v))' / numel(v);
        stairs(ax, [C.XLimLat(1); v], [0; y], '-', 'Color', st.Area{j}, 'LineWidth', 2);
        md = median(v);
        plot(ax, [md md], [0 0.5], ':', 'Color', st.Area{j}, 'LineWidth', 1.5);
        text(ax, C.XLimLat(2), 0.30 - 0.12*(j-1), sprintf('%s median %.0f ms  n=%d', names{j}, md, numel(v)), ...
            'FontSize', st.TickSize - 1, 'Color', st.Ink, 'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'right');
    end
    yline(ax, 0.5, '-', 'Color', st.Ref, 'HandleVisibility', 'off');
    SlideStyleAxes(ax, st);
    set(ax, 'XLim', C.XLimLat, 'YLim', [0 1.05], 'YTick', 0:0.5:1);
    xlabel(ax, C.LatLabel, 'FontSize', st.TickSize, 'Color', st.Muted);
    ylabel(ax, 'Cumulative fraction', 'FontSize', st.LabelSize, 'Color', st.Ink);
    if ~hasPop
        title(ax, {C.Title, PString(C.pLat, 'median')}, 'FontSize', st.LabelSize, 'Color', st.Ink, 'FontWeight', 'bold');
    else
        title(ax, PString(C.pLat, 'single neurons'), 'FontSize', st.TickSize, 'Color', st.Muted, 'FontWeight', 'normal');
    end
end
SlideAreaLegend(fig, st, names, [0.89 0.50 0.10 0.16]);
end

function s = PString(p, what)
if isnan(p),      s = '';
elseif p < 0.001, s = sprintf('%s: p < .001', what);
else,             s = sprintf('%s: p = %.3f', what, p);
end
end
