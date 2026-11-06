function SlideSubtitle(fig, txt, st)
% SlideSubtitle  Gray method subtitle under the headline.
% Inputs: fig, txt, st (PresentationStyle). Output: none.
%
% Xuefei Yu, 2026-11-05
annotation(fig, 'textbox', [0.05 0.80 0.9 0.08], 'String', txt, 'EdgeColor', 'none', ...
    'FontSize', st.SubtitleSize, 'Color', st.Muted, 'FontName', st.FontName, ...
    'VerticalAlignment', 'top', 'Interpreter', 'tex');
end
