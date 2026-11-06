function fig = PlotObjLocAUCScatterBySubset(aucObj, pObj, aucLoc, pLoc, masks, names, alpha, titleStr, figNum)
% PlotObjLocAUCScatterBySubset  Object vs location AUC scatter, one panel per
% subset. Colors: red = object only, purple = location only, green = both,
% open gray = neither (p < alpha, ranksum).
%
% Inputs
%   aucObj, pObj, aucLoc, pLoc — N-vectors
%   masks — cell{K} of N-by-1 logical; names — cellstr{K}
%   alpha, titleStr, figNum
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-05

K = numel(masks);
fig = PrepareFigure(figNum, 'w', [100 100 380*K 400], 'Name', titleStr);
so = pObj(:) < alpha; sl = pLoc(:) < alpha;
cats = {so & ~sl, ~so & sl, so & sl, ~so & ~sl};
cols = {[0.85 0.15 0.15], [0.55 0.10 0.80], [0.20 0.60 0.20], [0.6 0.6 0.6]};
lbls = {'Obj only', 'Loc only', 'Both', 'Neither'};
for k = 1:K
    subplot(1, K, k); hold on;
    for c = [4 1 2 3]
        m = masks{k}(:) & cats{c};
        if c == 4
            scatter(aucObj(m), aucLoc(m), 14, cols{c}, 'DisplayName', sprintf('%s (%d)', lbls{c}, sum(m)));
        else
            scatter(aucObj(m), aucLoc(m), 18, cols{c}, 'filled', 'DisplayName', sprintf('%s (%d)', lbls{c}, sum(m)));
        end
    end
    xline(0.5, ':k', 'HandleVisibility', 'off'); yline(0.5, ':k', 'HandleVisibility', 'off');
    axis square; xlim([0 1]); ylim([0 1]);
    xlabel('Object value AUC'); ylabel('Location value AUC');
    title(sprintf('%s (n=%d)', names{k}, sum(masks{k})), 'Interpreter', 'none');
    legend('Location', 'southeast', 'Box', 'off', 'FontSize', 7);
    box off;
end
sgtitle(titleStr, 'Interpreter', 'none');
end
