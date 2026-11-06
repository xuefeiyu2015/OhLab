function s = SlidePString(p, names)
% SlidePString  Short area-comparison label for slides ('FEF vs SC  p = .012').
% Inputs: p (Holm-adjusted p), names (cellstr of 2). Output: s (char).
%
% Xuefei Yu, 2026-11-05
if isnan(p),      s = '';
elseif p < 0.001, s = sprintf('%s vs %s  p < .001', names{:});
elseif p < 0.05,  s = sprintf('%s vs %s  p = %.3f', names{:}, p);
else,             s = sprintf('%s vs %s  n.s.', names{:});
end
end
