function SlideAreaLegend(fig, st, names, pos)
% SlideAreaLegend  FEF / SC color key as swatch + label.
% Inputs: fig, st, names (cellstr), pos ([x y w h] normalized). Output: none.
%
% Xuefei Yu, 2026-11-05
lax = axes(fig, 'Position', pos, 'Visible', 'off'); hold(lax, 'on');
for j = 1:2
    patch(lax, [0 0.12 0.12 0], 1 - 0.35*j + [0 0 0.2 0.2], st.Area{j}, 'EdgeColor', 'none');
    text(lax, 0.18, 1 - 0.35*j + 0.1, names{j}, 'FontSize', st.LabelSize, 'Color', st.Ink, ...
        'VerticalAlignment', 'middle');
end
set(lax, 'XLim', [0 1], 'YLim', [0 1]);
end
