function fig = PlotOverlapStatsComparison(S, titleStr, figNum)
% PlotOverlapStatsComparison  Compare object-location relationship across
% conditions (rows of S): Spearman rho, overlap odds ratio, sign concordance.
%
% Inputs
%   S        — table from ComputeOverlapStats (rows = conditions to compare)
%   titleStr, figNum
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-05

fig = PrepareFigure(figNum, 'w', [100 100 1000 360], 'Name', titleStr);
K = height(S); x = 1:K;
panels = {S.Rho, S.pRho, 0, 'Spearman \rho (AUC_{obj}, AUC_{loc})'; ...
          S.OddsRatio, S.pFisher, 1, 'Overlap odds ratio (Fisher)'; ...
          S.Concord, S.pConcord, 0.5, 'Fraction same preferred value'};
for k = 1:3
    subplot(1, 3, k); hold on;
    v = panels{k,1}; p = panels{k,2};
    bar(x, v, 0.6, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'BaseValue', panels{k,3});
    yline(panels{k,3}, '--k');
    for i = 1:K
        if isnan(v(i)), continue; end
        text(i, v(i), sprintf('%.2f\np=%.2g\nn=%d', v(i), p(i), S.n(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
    yl = ylim; ylim([yl(1), yl(2) + 0.25*diff(yl)]);
    set(gca, 'XTick', x, 'XTickLabel', S.Subset, 'TickLabelInterpreter', 'none', 'XTickLabelRotation', 20);
    title(panels{k,4});
    box off;
end
sgtitle(titleStr, 'Interpreter', 'none');
end
