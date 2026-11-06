function fig = PlotPreTargetAUC(P, titleStr, figNum)
% PlotPreTargetAUC  Grouped bars: pre-target AUC (each task aligned by its own
% 100-200 ms preference) per group, mean +/- SEM, chance line at 0.5.
%
% Inputs
%   P        — table from ComputePreTargetAUC
%   titleStr — figure title; figNum — figure number
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-05

fig = PrepareFigure(figNum, 'w', [100 100 600 380], 'Name', titleStr);
hold on;
G = height(P);
M = [P.ObjMean, P.LocMean];
S = [P.ObjSEM,  P.LocSEM];
pv = [P.ObjP, P.LocP];
b = bar(1:G, M, 0.8, 'EdgeColor', 'none', 'BaseValue', 0.5);
b(1).FaceColor = [0.85 0.15 0.15];
b(2).FaceColor = [0.55 0.10 0.80];
for k = 1:2
    x = b(k).XEndPoints;
    errorbar(x, M(:,k), S(:,k), 'k.', 'LineWidth', 1, 'HandleVisibility', 'off');
    for g = 1:G
        if ~isnan(pv(g,k)) && pv(g,k) < 0.05
            text(x(g), M(g,k) + S(g,k) + 0.005, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
        end
    end
end
yline(0.5, '--k', 'HandleVisibility', 'off');
set(gca, 'XTick', 1:G, 'XTickLabel', P.Group, 'TickLabelInterpreter', 'none');
ylabel(sprintf('AUC in [%d, %d] ms (own-sign aligned)', P.Properties.UserData.Window));
legend({'Object (stable)', 'Location (oneDR)'}, 'Box', 'off', 'Location', 'northeast');
title(titleStr, 'Interpreter', 'none');
box off;
end
