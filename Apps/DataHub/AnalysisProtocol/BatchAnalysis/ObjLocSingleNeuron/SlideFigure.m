function fig = SlideFigure(figNum, headline, st)
% SlideFigure  New slide figure with a bold headline (the conclusion).
% Inputs: figNum, headline (char), st (PresentationStyle). Output: fig.
%
% Xuefei Yu, 2026-11-05
fig = figure(figNum); clf(fig);
set(fig, 'Color', st.Surface, 'Position', st.FigPos, 'Name', headline);
annotation(fig, 'textbox', [0.05 0.88 0.9 0.1], 'String', headline, 'EdgeColor', 'none', ...
    'FontSize', st.HeadlineSize, 'FontWeight', 'bold', 'Color', st.Ink, 'FontName', st.FontName, ...
    'VerticalAlignment', 'middle', 'Interpreter', 'none');
end
