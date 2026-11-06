function C = ComputeSignConcordance(aucObj, aucLoc, groupMasks, groupNames)
% ComputeSignConcordance  Per group, fraction of neurons whose object and
% location AUC fall on the same side of 0.5 (same preferred value).
%
% Inputs
%   aucObj, aucLoc — N-vectors of scalar AUC (100-200 ms)
%   groupMasks     — cell{G} of N-by-1 logical masks
%   groupNames     — cellstr{G}
% Output: table with one row per group
%   Group, n, nConcordant, Fraction, pBinom (two-sided exact binomial vs 0.5)
%   Neurons with either AUC exactly 0.5 or NaN are excluded.
%
% Xuefei Yu, 2026-11-05

aucObj = aucObj(:); aucLoc = aucLoc(:);
G = numel(groupMasks);
n = zeros(G,1); k = zeros(G,1); frac = nan(G,1); p = nan(G,1);
sObj = sign(aucObj - 0.5); sLoc = sign(aucLoc - 0.5);
for g = 1:G
    m = groupMasks{g}(:) & sObj ~= 0 & sLoc ~= 0 & ~isnan(sObj) & ~isnan(sLoc);
    n(g) = sum(m);
    k(g) = sum(sObj(m) == sLoc(m));
    if n(g) > 0
        frac(g) = k(g) / n(g);
        p(g)    = BinomTwoSided(k(g), n(g));
    end
end
C = table(groupNames(:), n, k, frac, p, ...
    'VariableNames', {'Group','n','nConcordant','Fraction','pBinom'});
end

function p = BinomTwoSided(k, n)
% Exact two-sided binomial test against 0.5.
p = min(1, 2 * min(binocdf(k, n, 0.5), 1 - binocdf(k - 1, n, 0.5)));
end
