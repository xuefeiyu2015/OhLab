function C = ComputeAreaComparison(SA, SB, names)
% ComputeAreaComparison  Formal tests of FEF vs SC (or any two areas) on the
% single-unit summary metrics, with Holm correction across the family.
%
% Inputs
%   SA, SB — structs from ComputeSingleUnitSummary (area A, area B)
%   names  — {'FEF','SC'}
% Output: table C, one row per measure
%   Measure, A, B (values), Diff (B - A), Test, p, pHolm
%   Tests: fractions -> Fisher exact; pre-target d -> ranksum;
%          log odds ratio -> z (Woolf SE, +0.5 if a zero cell);
%          Spearman rho -> Fisher z (SE = sqrt(1.06/(n-3)));
%          concordance -> Fisher exact; noise-corrected r -> bootstrap.
%
% Xuefei Yu, 2026-11-05

rows = {};
% --- coding / baseline-difference fractions
for i = 1:height(SA.Coding)
    kA = SA.Coding.k(i); nA = SA.Coding.n(i); kB = SB.Coding.k(i); nB = SB.Coding.n(i);
    [~, p] = fishertest([kA nA-kA; kB nB-kB]);
    rows(end+1,:) = {['Frac ' SA.Coding.Measure{i}], kA/nA, kB/nB, 'Fisher', p}; %#ok<AGROW>
end
% --- pre-target aligned baseline AUC
for tk = {'Obj', 'Loc'}
    p = ranksum(SA.PreTgtD.(tk{1}), SB.PreTgtD.(tk{1}));
    rows(end+1,:) = {[tk{1} ' pre-target aligned AUC-0.5'], mean(SA.PreTgtD.(tk{1})), mean(SB.PreTgtD.(tk{1})), 'ranksum', p}; %#ok<AGROW>
end
% --- relationship, raw and corrected
for r = 1:2
    v = SA.Overlap.Subset{r};
    [lA, sA] = LogOR(SA.Overlap(r,:)); [lB, sB] = LogOR(SB.Overlap(r,:));
    p = 2*normcdf(-abs(lA - lB) / sqrt(sA^2 + sB^2));
    rows(end+1,:) = {['Overlap log OR (' v ')'], lA, lB, 'z (log OR)', p}; %#ok<AGROW>

    rA = SA.Overlap.Rho(r); rB = SB.Overlap.Rho(r); nA = SA.Overlap.n(r); nB = SB.Overlap.n(r);
    z = (atanh(rA) - atanh(rB)) / sqrt(1.06/(nA-3) + 1.06/(nB-3));
    rows(end+1,:) = {['Spearman rho (' v ')'], rA, rB, 'Fisher z', 2*normcdf(-abs(z))}; %#ok<AGROW>

    cA = SA.Concord.(v); cB = SB.Concord.(v);
    [~, p] = fishertest([cA.nConcordant cA.n-cA.nConcordant; cB.nConcordant cB.n-cB.nConcordant]);
    rows(end+1,:) = {['Concordance (' v ')'], cA.Fraction, cB.Fraction, 'Fisher', p}; %#ok<AGROW>

    RA = SA.Rel.(v); RB = SB.Rel.(v);
    d = RB.rDisBoot - RA.rDisBoot; d = d(~isnan(d));
    p = min(1, 2*min(mean(d <= 0), mean(d >= 0)));
    rows(end+1,:) = {['Noise-corrected r (' v ')'], RA.rDis, RB.rDis, 'bootstrap', p}; %#ok<AGROW>
end

C = cell2table(rows, 'VariableNames', {'Measure', names{1}, names{2}, 'Test', 'p'});
C.Diff  = C.(names{2}) - C.(names{1});
C.pHolm = Holm(C.p);
C = C(:, {'Measure', names{1}, names{2}, 'Diff', 'Test', 'p', 'pHolm'});
end

function [l, s] = LogOR(row)
a = row.nBoth; b = row.nObj - row.nBoth; c = row.nLoc - row.nBoth; d = row.n - row.nObj - row.nLoc + row.nBoth;
x = [a b c d];
if any(x == 0), x = x + 0.5; end
l = log(x(1)*x(4) / (x(2)*x(3)));
s = sqrt(sum(1 ./ x));
end

function q = Holm(p)
[ps, o] = sort(p(:)); m = numel(ps);
adj = min(1, cummax((m - (1:m)' + 1) .* ps));
q = nan(size(p)); q(o) = adj;
end
