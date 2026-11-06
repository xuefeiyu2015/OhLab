function st = PresentationStyle()
% PresentationStyle  Shared look for slide figures (constants only).
% Area colors validated for CVD separation / contrast (dataviz validator):
% FEF #2a78d6 (blue), SC #eb6834 (orange).
%
% Output: struct st with colors, font sizes, figure size
%
% Xuefei Yu, 2026-11-05

hex = @(h) sscanf(h(2:end), '%2x%2x%2x', [1 3]) / 255;
st.Area    = {hex('#2a78d6'), hex('#eb6834')};   % FEF, SC
st.Ink     = hex('#1f1f1e');                     % primary text
st.Muted   = hex('#6b6a64');                     % subtitles, reference labels
st.Ref     = hex('#9a9992');                     % reference lines
st.Surface = [1 1 1];
st.FontName     = 'Helvetica';
st.HeadlineSize = 18;
st.SubtitleSize = 11;
st.LabelSize    = 14;
st.TickSize     = 12;
st.FigPos       = [100 100 1000 600];            % 16:9-ish slide panel
end
