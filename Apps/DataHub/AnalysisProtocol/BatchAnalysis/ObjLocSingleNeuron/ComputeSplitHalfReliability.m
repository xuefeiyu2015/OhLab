function R = ComputeSplitHalfReliability(TrialsObj, tObj, TrialsLoc, tLoc, mask, Params)
% ComputeSplitHalfReliability  Noise ceiling for the object-location AUC
% correlation. Splits each neuron's RF trials in half (stratified by value),
% computes the AUC in each half, and corrects the across-neuron correlation
% of object vs location AUC for measurement noise:
%   r_dis = r(AUC_obj, AUC_loc) / sqrt(rel_obj * rel_loc)
%   rel   = Spearman-Brown( r(AUC_half1, AUC_half2) ), averaged over splits
% Pearson correlations of (AUC - 0.5) are used throughout.
%
% Inputs
%   TrialsObj/TrialsLoc — N-by-1 structs (.Rate nTrial-by-T, .IsHigh)
%   tObj/tLoc           — 1-by-T time vectors (ms)
%   mask                — N-by-1 logical neurons to include
%   Params — .RespWin, .BaseWin, .Variant ('corrected' | 'raw'),
%            .nSplit (default 100), .nBoot (default 2000), .Seed
% Output: struct R
%   .n, .rObs, .relObj, .relLoc, .rDis, .rDisCI (95% bootstrap over neurons)
%   .rDisBoot (nBoot-by-1), .rObsBoot, .relObjBoot, .relLocBoot
%
% Xuefei Yu, 2026-11-05

if ~isfield(Params,'Variant'), Params.Variant = 'corrected'; end
if ~isfield(Params,'nSplit'),  Params.nSplit = 100;  end
if ~isfield(Params,'nBoot'),   Params.nBoot = 2000;  end
if ~isfield(Params,'Seed'),    Params.Seed = 1;      end
rs = RandStream('mt19937ar', 'Seed', Params.Seed);

idx = find(mask(:));
[vO, hO] = TrialMetric(TrialsObj(idx), tObj, Params);
[vL, hL] = TrialMetric(TrialsLoc(idx), tLoc, Params);
ok = cellfun(@(h) sum(h) >= 2 && sum(~h) >= 2, hO) & cellfun(@(h) sum(h) >= 2 && sum(~h) >= 2, hL);
vO = vO(ok); hO = hO(ok); vL = vL(ok); hL = hL(ok);
n = numel(vO);

fullO = cellfun(@(v,h) MannWhitneyAUC(v(h), v(~h)), vO, hO) - 0.5;
fullL = cellfun(@(v,h) MannWhitneyAUC(v(h), v(~h)), vL, hL) - 0.5;

S = Params.nSplit;
[aO, bO, aL, bL] = deal(nan(n, S));
for s = 1:S
    for k = 1:n
        [aO(k,s), bO(k,s)] = HalfAUC(vO{k}, hO{k}, rs);
        [aL(k,s), bL(k,s)] = HalfAUC(vL{k}, hL{k}, rs);
    end
end
aO = aO - 0.5; bO = bO - 0.5; aL = aL - 0.5; bL = bL - 0.5;

R.n      = n;
R.rObs   = PCorr(fullO, fullL);
R.relObj = mean(arrayfun(@(s) SB(PCorr(aO(:,s), bO(:,s))), 1:S));
R.relLoc = mean(arrayfun(@(s) SB(PCorr(aL(:,s), bL(:,s))), 1:S));
R.rDis   = R.rObs / sqrt(max(R.relObj, eps) * max(R.relLoc, eps));

nB = Params.nBoot;
[R.rDisBoot, R.rObsBoot, R.relObjBoot, R.relLocBoot] = deal(nan(nB,1));
for b = 1:nB
    ii = randi(rs, n, n, 1);
    s  = randi(rs, S);
    R.rObsBoot(b)   = PCorr(fullO(ii), fullL(ii));
    R.relObjBoot(b) = SB(PCorr(aO(ii,s), bO(ii,s)));
    R.relLocBoot(b) = SB(PCorr(aL(ii,s), bL(ii,s)));
    if R.relObjBoot(b) > 0 && R.relLocBoot(b) > 0
        R.rDisBoot(b) = R.rObsBoot(b) / sqrt(R.relObjBoot(b) * R.relLocBoot(b));
    end
end
R.rDisCI = prctile(R.rDisBoot(~isnan(R.rDisBoot)), [2.5 97.5]);
R.Params = Params;
end

% ---------------------------------------------------------------------------
function [v, h] = TrialMetric(Trials, t, P)
iR = t >= P.RespWin(1) & t <= P.RespWin(2);
iB = t >= P.BaseWin(1) & t <= P.BaseWin(2);
v = cell(numel(Trials),1); h = cell(numel(Trials),1);
for k = 1:numel(Trials)
    X = Trials(k).Rate;
    if isempty(X), v{k} = []; h{k} = false(0,1); continue; end
    r = mean(X(:,iR), 2, 'omitnan');
    if strcmpi(P.Variant, 'corrected')
        r = r - mean(X(:,iB), 2, 'omitnan');
    end
    keep = ~isnan(r);
    v{k} = r(keep); h{k} = Trials(k).IsHigh(keep);
end
end

function [a, b] = HalfAUC(v, h, rs)
iH = find(h); iL = find(~h);
iH = iH(randperm(rs, numel(iH))); iL = iL(randperm(rs, numel(iL)));
nh = floor(numel(iH)/2); nl = floor(numel(iL)/2);
a = MannWhitneyAUC(v(iH(1:nh)),     v(iL(1:nl)));
b = MannWhitneyAUC(v(iH(nh+1:end)), v(iL(nl+1:end)));
end

function r = PCorr(x, y)
ok = ~isnan(x) & ~isnan(y);
r = corr(x(ok), y(ok));
end

function r = SB(r)
r = 2*r / (1 + r);
end
