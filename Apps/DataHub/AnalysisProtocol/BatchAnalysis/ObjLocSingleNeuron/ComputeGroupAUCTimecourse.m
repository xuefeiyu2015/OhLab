function TC = ComputeGroupAUCTimecourse(A, groupMasks, alpha)
% ComputeGroupAUCTimecourse  Group mean / SEM of AUC time courses and
% per-bin significance vs 0.5 (signrank across neurons).
%
% Inputs
%   A          — N-by-T AUC matrix (already sign-aligned)
%   groupMasks — cell{G} of N-by-1 logical masks
%   alpha      — significance level (default 0.05)
% Output: struct array TC(G) with fields
%   .n, .mean (1-by-T), .sem (1-by-T), .p (1-by-T), .sig (1-by-T logical)
%
% Xuefei Yu, 2026-11-05

if nargin < 3, alpha = 0.05; end
nT = size(A,2);
TC = repmat(struct('n',0,'mean',nan(1,nT),'sem',nan(1,nT),'p',nan(1,nT),'sig',false(1,nT)), ...
            1, numel(groupMasks));
for g = 1:numel(groupMasks)
    X = A(groupMasks{g},:);
    TC(g).n = size(X,1);
    if TC(g).n == 0, continue; end
    nValid      = sum(~isnan(X),1);
    TC(g).mean  = mean(X,1,'omitnan');
    TC(g).sem   = std(X,0,1,'omitnan') ./ sqrt(nValid);
    for t = 1:nT
        v = X(~isnan(X(:,t)),t);
        if numel(v) >= 3
            TC(g).p(t) = signrank(v - 0.5);
        end
    end
    TC(g).sig = TC(g).p < alpha;
end
end
