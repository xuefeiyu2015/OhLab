function fig = PlotSignConcordance(C, titleStr, figNum)
% PlotSignConcordance  Bar plot of the fraction of neurons with the same
% preferred value (object vs location) per group, 50% chance line.
%
% Inputs
%   C        — table from ComputeSignConcordance
%   titleStr — figure title; figNum — figure number
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-05

fig = PrepareFigure(figNum, 'w', [100 100 520 380], 'Name', titleStr);
hold on;
G = height(C);
bar(1:G, C.Fraction, 0.6, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
yline(0.5, '--k');
for g = 1:G
    if isnan(C.Fraction(g)), continue; end
    text(g, C.Fraction(g) + 0.03, sprintf('%d/%d\np=%.3g', C.nConcordant(g), C.n(g), C.pBinom(g)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
end
set(gca, 'XTick', 1:G, 'XTickLabel', C.Group, 'TickLabelInterpreter', 'none');
ylim([0 1.15]); ylabel('Fraction same preferred value');
title(titleStr, 'Interpreter', 'none');
box off;
end
