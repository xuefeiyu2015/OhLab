function SlideStyleAxes(ax, st)
% SlideStyleAxes  Recessive axes for slides.
% Inputs: ax, st (PresentationStyle). Output: none.
%
% Xuefei Yu, 2026-11-05
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', st.TickSize, 'FontName', st.FontName, ...
    'XColor', st.Muted, 'YColor', st.Muted, 'LineWidth', 1);
end
