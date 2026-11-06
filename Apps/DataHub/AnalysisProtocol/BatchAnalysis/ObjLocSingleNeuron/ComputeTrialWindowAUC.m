function C = ComputeTrialWindowAUC(Trials, t, mask, Params)
% ComputeTrialWindowAUC  High- vs low-value ROC from single RF trials, with and
% without trial-wise pre-target baseline subtraction.
%
% Inputs
%   Trials — N-by-1 struct with .Rate (nTrial-by-T FR) and .IsHigh (logical)
%   t      — 1-by-T time vector (ms from target onset)
%   mask   — N-by-1 logical; neurons not in mask return NaN
%   Params — .RespWin (default [100 200]) target-response window (ms)
%            .BaseWin (default [-300 0])  pre-target baseline window (ms)
%            .TCSmoothMs (default 50) boxcar width for the per-bin AUC time
%                     courses only (10 ms bins are too sparse); scalar
%                     window AUCs are never smoothed
% Output: struct C (all N-by-1 unless noted)
%   .AUCRaw,  .pRaw   — response-window AUC / ranksum p, no correction
%   .AUCCorr, .pCorr  — response minus the same trial's baseline
%   .AUCBase, .pBase  — baseline window alone (the "baseline difference")
%   .AUCPSTHRaw, .AUCPSTHCorr — N-by-T per-bin AUC (corrected: bin minus trial baseline)
%   .Params
%
% Xuefei Yu, 2026-11-05

if ~isfield(Params,'RespWin'), Params.RespWin = [100 200]; end
if ~isfield(Params,'BaseWin'), Params.BaseWin = [-300 0];  end
if ~isfield(Params,'TCSmoothMs'), Params.TCSmoothMs = 50;  end
nSm = max(1, round(Params.TCSmoothMs / median(diff(t))));

N  = numel(Trials); nT = numel(t);
iR = t >= Params.RespWin(1) & t <= Params.RespWin(2);
iB = t >= Params.BaseWin(1) & t <= Params.BaseWin(2);

C.AUCRaw = nan(N,1);  C.pRaw  = nan(N,1);
C.AUCCorr = nan(N,1); C.pCorr = nan(N,1);
C.AUCBase = nan(N,1); C.pBase = nan(N,1);
C.AUCPSTHRaw  = nan(N,nT);
C.AUCPSTHCorr = nan(N,nT);

for n = find(mask(:))'
    X  = Trials(n).Rate;  hi = Trials(n).IsHigh;
    if isempty(X) || ~any(hi) || all(hi), continue; end
    base = mean(X(:,iB), 2, 'omitnan');
    resp = mean(X(:,iR), 2, 'omitnan');
    corr = resp - base;

    [C.AUCRaw(n),  C.pRaw(n)]  = AucP(resp(hi), resp(~hi));
    [C.AUCCorr(n), C.pCorr(n)] = AucP(corr(hi), corr(~hi));
    [C.AUCBase(n), C.pBase(n)] = AucP(base(hi), base(~hi));

    Xs = movmean(X, nSm, 2, 'omitnan');
    Xc = Xs - base;
    for k = 1:nT
        C.AUCPSTHRaw(n,k)  = MannWhitneyAUC(Xs(hi,k), Xs(~hi,k));
        C.AUCPSTHCorr(n,k) = MannWhitneyAUC(Xc(hi,k), Xc(~hi,k));
    end
end
C.Params = Params;
end

function [a, p] = AucP(x, y)
x = x(~isnan(x)); y = y(~isnan(y));
a = MannWhitneyAUC(x, y); p = NaN;
if numel(x) >= 3 && numel(y) >= 3
    p = ranksum(x, y);
end
end
