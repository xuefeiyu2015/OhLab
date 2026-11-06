function O = ComputeOverlapVsChance(aucObj, pObj, aucLoc, pLoc, alpha)
% ComputeOverlapVsChance  Venn counts (object only / both / location only /
% neither) versus what independence predicts given the marginal counts.
%
% Inputs
%   aucObj, pObj, aucLoc, pLoc — n-vectors for the neurons to analyze
%   alpha                      — coding threshold (p < alpha)
% Output: struct O
%   .N, .nObj, .nLoc
%   .Labels  {'Object only','Both','Location only','Neither'}
%   .Obs     1-by-4 observed counts
%   .Exp     1-by-4 expected under independence (E[both] = nObj*nLoc/N)
%   .NullLo, .NullHi  1-by-4 central 95% range under independence
%            (hypergeometric for "both"; others follow from the fixed margins)
%   .OR, .pFisher, .Verdict ('above' | 'chance' | 'below': Fisher p < alpha
%            and direction of "both" vs expected)
%   .Rho, .pRho  Spearman (AUC_obj-0.5, AUC_loc-0.5)
%   .HighLow  counts of significant high- / low-preferring neurons:
%             [objHigh objLow locHigh locLow]
%
% Xuefei Yu, 2026-11-05

ok = ~isnan(aucObj) & ~isnan(pObj) & ~isnan(aucLoc) & ~isnan(pLoc);
ao = aucObj(ok); al = aucLoc(ok); so = pObj(ok) < alpha; sl = pLoc(ok) < alpha;
N = numel(ao); nO = sum(so); nL = sum(sl); b = sum(so & sl);

O.N = N; O.nObj = nO; O.nLoc = nL;
O.Labels = {'Object only', 'Both', 'Location only', 'Neither'};
O.Obs = [nO - b, b, nL - b, N - nO - nL + b];
eb = nO * nL / N;
O.Exp = [nO - eb, eb, nL - eb, N - nO - nL + eb];
bLo = hygeinv(0.025, N, nO, nL); bHi = hygeinv(0.975, N, nO, nL);
O.NullLo = [nO - bHi, bLo, nL - bHi, N - nO - nL + bLo];
O.NullHi = [nO - bLo, bHi, nL - bLo, N - nO - nL + bHi];

[~, O.pFisher, st] = fishertest([b, nO - b; nL - b, N - nO - nL + b]);
O.OR = st.OddsRatio;
% Verdict follows the Fisher test (so headline and reported p always agree)
if O.pFisher < alpha && b > eb,     O.Verdict = 'above';
elseif O.pFisher < alpha && b < eb, O.Verdict = 'below';
else,                               O.Verdict = 'chance';
end
[O.Rho, O.pRho] = corr(ao - 0.5, al - 0.5, 'type', 'Spearman');
O.HighLow = [sum(so & ao > 0.5), sum(so & ao < 0.5), sum(sl & al > 0.5), sum(sl & al < 0.5)];

% log OR and Woolf SE (for area comparisons)
x = [b, nO - b, nL - b, N - nO - nL + b];
if any(x == 0), x = x + 0.5; end
O.LogOR   = log(x(1)*x(4) / (x(2)*x(3)));
O.LogORSE = sqrt(sum(1 ./ x));
end
