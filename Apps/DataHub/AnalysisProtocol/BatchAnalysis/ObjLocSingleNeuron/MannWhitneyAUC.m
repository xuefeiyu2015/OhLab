function auc = MannWhitneyAUC(x, y)
% MannWhitneyAUC  Exact ROC area P(x > y) + 0.5 P(x == y), NaNs dropped.
% Same quantity rocN estimates (x = positive / high value class).
%
% Inputs : x, y — vectors (high, low)
% Output : auc  — scalar, NaN if either class is empty
%
% Xuefei Yu, 2026-11-05

x = x(~isnan(x)); y = y(~isnan(y));
nx = numel(x); ny = numel(y);
if nx == 0 || ny == 0
    auc = NaN; return;
end
r = tiedrank([x(:); y(:)]);
auc = (sum(r(1:nx)) - nx*(nx+1)/2) / (nx*ny);
end
