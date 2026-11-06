function fig = PlotOtherValueAUCDistribution(E, titleStr, figNum)
% PlotOtherValueAUCDistribution  Histograms of the sign-aligned "other"
% value AUC deviation d = (AUC_other-0.5)*sign(AUC_ref-0.5) for ObjOnly and
% LocOnly neurons, with the equivalence band shaded.
%
% Inputs
%   E        — table from ComputeOtherValueEquivalence
%   titleStr — figure title; figNum — figure number
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-05

fig = PrepareFigure(figNum, 'w', [100 100 800 340], 'Name', titleStr);
m = E.Properties.UserData.Margin;
cols = {[0.55 0.10 0.80], [0.85 0.15 0.15]};   % other = Loc (purple), Obj (red)
edges = -0.5:0.025:0.5;
for r = 1:2
    subplot(1, 2, r); hold on;
    yl = [0 1];
    d = E.Properties.UserData.d{r};
    if ~isempty(d)
        histogram(d, edges, 'FaceColor', cols{r}, 'EdgeColor', 'none');
        yl = ylim;
    end
    patch([-m m m -m], [yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');
    xline(0, '--k');
    if ~isnan(E.AlignedMean(r))
        xline(E.AlignedMean(r), '-', 'Color', cols{r}, 'LineWidth', 2);
    end
    ylim(yl); xlim([-0.5 0.5]);
    xlabel('Other AUC - 0.5 (aligned to reference preference)');
    ylabel('Neurons');
    title(sprintf('%s  n=%d\nmean=%.3f [%.3f %.3f], p=%.3g; TOST p=%.3g', ...
        E.Group{r}, E.n(r), E.AlignedMean(r), E.AlignedCI(r,1), E.AlignedCI(r,2), ...
        E.pAligned(r), E.pTOST(r)), 'FontSize', 8, 'Interpreter', 'none');
    box off;
end
sgtitle(sprintf('%s  (band = \\pm%.2f)', titleStr, m));
end
