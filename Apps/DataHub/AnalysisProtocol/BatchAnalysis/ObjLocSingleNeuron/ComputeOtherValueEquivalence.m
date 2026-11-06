function E = ComputeOtherValueEquivalence(aucObj, aucLoc, G, margin, nBoot, seed)
% ComputeOtherValueEquivalence  For ObjOnly / LocOnly neurons, test whether the
% "other" value AUC is really at chance or just below the significance threshold.
%
% Inputs
%   aucObj, aucLoc — N-vectors, scalar 100-200 ms AUC
%   G              — struct from ComputeValueGroups
%   margin         — equivalence margin on |AUC-0.5| (e.g. 0.05)
%   nBoot          — bootstrap resamples (e.g. 2000)
%   seed           — RNG seed (local stream; global RNG is untouched)
% Output: table, rows = {'ObjOnly (other=Loc)','LocOnly (other=Obj)'}
%   n
%   AlignedMean, AlignedCI  — mean of d = (AUC_other-0.5)*sign(AUC_ref-0.5),
%                             d > 0 means same preferred value as the reference task
%   pAligned               — signrank of d vs 0
%   AbsMean, AbsCI          — mean |AUC_other-0.5| (unsigned strength)
%   FracWithinMargin       — fraction with |AUC_other-0.5| < margin
%   pTOST                  — TOST of d within [-margin, margin] (p < .05 => equivalent to chance)
%   E.Properties.UserData.d = {d_ObjOnly, d_LocOnly}
%
% Xuefei Yu, 2026-11-05

aucObj = aucObj(:); aucLoc = aucLoc(:);
rs = RandStream('mt19937ar', 'Seed', seed);

rows  = {'ObjOnly (other=Loc)'; 'LocOnly (other=Obj)'};
masks = {G.Mask.ObjOnly, G.Mask.LocOnly};
other = {aucLoc, aucObj};
refer = {aucObj, aucLoc};

[n, am, pA, bm, fw, pT] = deal(nan(2,1));
[aCI, bCI] = deal(nan(2,2));
dAll = cell(2,1);
for r = 1:2
    m = masks{r};
    d = (other{r}(m) - 0.5) .* sign(refer{r}(m) - 0.5);
    a = abs(other{r}(m) - 0.5);
    dAll{r} = d;
    n(r) = numel(d);
    if n(r) < 3, continue; end
    am(r) = mean(d);  bm(r) = mean(a);
    aCI(r,:) = BootCI(d, nBoot, rs);
    bCI(r,:) = BootCI(a, nBoot, rs);
    pA(r) = signrank(d);
    fw(r) = mean(a < margin);
    [~, p1] = ttest(d, -margin, 'Tail', 'right');
    [~, p2] = ttest(d,  margin, 'Tail', 'left');
    pT(r) = max(p1, p2);
end
E = table(rows, n, am, aCI, pA, bm, bCI, fw, pT, 'VariableNames', ...
    {'Group','n','AlignedMean','AlignedCI','pAligned','AbsMean','AbsCI','FracWithinMargin','pTOST'});
E.Properties.UserData.d = dAll;
E.Properties.UserData.Margin = margin;
end

function ci = BootCI(x, nBoot, rs)
nx = numel(x);
idx = randi(rs, nx, nx, nBoot);
mu = mean(x(idx), 1);
ci = prctile(mu, [2.5 97.5]);
end
