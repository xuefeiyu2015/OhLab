function fig = PlotRawVsCorrectedAUC(CObj, CLoc, mask, alpha, titleStr, figNum)
% PlotRawVsCorrectedAUC  Effect of baseline subtraction, per task.
% Top row   : raw vs baseline-corrected response AUC (object | location).
% Bottom row: distribution of the baseline-window AUC (filled = p < alpha).
%
% Inputs
%   CObj, CLoc — structs from ComputeTrialWindowAUC
%   mask       — N-by-1 logical neurons to show
%   alpha      — significance level for the filled histogram
%   titleStr, figNum
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-05

fig = PrepareFigure(figNum, 'w', [100 100 820 700], 'Name', titleStr);
Cs = {CObj, CLoc}; lbl = {'Object value (stable)', 'Location value (oneDR)'};
cols = {[0.85 0.15 0.15], [0.55 0.10 0.80]};
edges = 0:0.025:1;
for k = 1:2
    C = Cs{k}; m = mask(:);
    subplot(2, 2, k); hold on;
    scatter(C.AUCRaw(m), C.AUCCorr(m), 14, cols{k}, 'filled', 'MarkerFaceAlpha', 0.5);
    plot([0 1], [0 1], '-', 'Color', [0.5 0.5 0.5]);
    xline(0.5, ':k'); yline(0.5, ':k');
    axis square; xlim([0 1]); ylim([0 1]);
    xlabel('Response AUC, raw'); ylabel('Response AUC, baseline-corrected');
    title(sprintf('%s\nn sig: raw %d -> corrected %d', lbl{k}, ...
        sum(C.pRaw(m) < alpha), sum(C.pCorr(m) < alpha)));
    box off;

    subplot(2, 2, k + 2); hold on;
    histogram(C.AUCBase(m), edges, 'FaceColor', 'none', 'EdgeColor', cols{k});
    histogram(C.AUCBase(m & C.pBase < alpha), edges, 'FaceColor', cols{k}, 'EdgeColor', 'none');
    xline(0.5, '--k');
    xlabel(sprintf('Baseline AUC [%d %d] ms', C.Params.BaseWin)); ylabel('Neurons');
    title(sprintf('Baseline difference: %d / %d sig', sum(m & C.pBase < alpha), sum(m & ~isnan(C.pBase))));
    box off;
end
sgtitle(titleStr, 'Interpreter', 'none');
end
