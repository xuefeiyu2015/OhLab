function S = ComputeOverlapStats(aucObj, pObj, aucLoc, pLoc, masks, names, alpha)
% ComputeOverlapStats  Object vs location value relationship within subsets.
%
% Inputs
%   aucObj, pObj, aucLoc, pLoc — N-vectors (scalar AUC and ranksum p)
%   masks — cell{K} of N-by-1 logical subsets; names — cellstr{K}
%   alpha — significance level (default 0.05)
% Output: table, one row per subset
%   n, nObj, nLoc, nBoth, ExpBoth (if independent), OddsRatio, pFisher
%   Rho, pRho        — Spearman corr of (AUC_obj-0.5) vs (AUC_loc-0.5)
%   Concord, pConcord — fraction with the same preferred value (binomial vs 0.5)
%   OppBoth          — # "Both" neurons with opposite preferred value
%
% Xuefei Yu, 2026-11-05

if nargin < 7, alpha = 0.05; end
K = numel(masks);
z = nan(K,1);
[n, nO, nL, nB, eB, OR, pF, rho, pR, cc, pC, opp] = deal(z,z,z,z,z,z,z,z,z,z,z,z);
for k = 1:K
    m  = masks{k}(:) & ~isnan(aucObj(:)) & ~isnan(aucLoc(:)) & ~isnan(pObj(:)) & ~isnan(pLoc(:));
    ao = aucObj(m); al = aucLoc(m); so = pObj(m) < alpha; sl = pLoc(m) < alpha;
    n(k) = sum(m);
    if n(k) < 3, continue; end
    nO(k) = sum(so); nL(k) = sum(sl); nB(k) = sum(so & sl);
    eB(k) = n(k) * mean(so) * mean(sl);
    [~, pF(k), st] = fishertest([sum(so&sl) sum(so&~sl); sum(~so&sl) sum(~so&~sl)]);
    OR(k) = st.OddsRatio;
    [rho(k), pR(k)] = corr(ao - 0.5, al - 0.5, 'type', 'Spearman');
    C = ComputeSignConcordance(ao, al, {true(n(k),1)}, {'x'});
    cc(k) = C.Fraction; pC(k) = C.pBinom;
    b = so & sl;
    opp(k) = sum(sign(ao(b) - 0.5) ~= sign(al(b) - 0.5));
end
S = table(names(:), n, nO, nL, nB, eB, OR, pF, rho, pR, cc, pC, opp, 'VariableNames', ...
    {'Subset','n','nObj','nLoc','nBoth','ExpBoth','OddsRatio','pFisher','Rho','pRho','Concord','pConcord','OppBoth'});
end
