function L = ComputeValueLatency(Trials, t, mask, Params)
% ComputeValueLatency  Single-neuron onset latency of the value signal (high vs
% low) and of the visual response, from single RF trials.
%
% Firing rates are smoothed with a CAUSAL boxcar (each bin = mean of the
% current and preceding bins), so no post-onset activity leaks backward and
% latencies are not biased early. Value latency = first bin in LatSearchWin
% that starts a run of >= LatConsecBins bins with high-vs-low Mann-Whitney
% p < LatAlpha AND the same AUC sign as the neuron's response-window AUC.
% Two variants: Raw (rates as recorded) and Corr (each trial minus its own
% BaseWin mean — removes the oneDR pre-target anticipation offset).
% Visual latency = first such run where (rate - trial baseline) > 0 across
% all RF trials (signrank).
%
% Inputs
%   Trials — N-by-1 struct with .Rate (nTrial-by-T FR, 10 ms bins) and .IsHigh
%   t      — 1-by-T time (ms from target onset)
%   mask   — N-by-1 logical neurons to compute (others NaN)
%   Params — .RespWin [100 200], .BaseWin [-300 0], .LatSmoothMs (30),
%            .LatConsecBins (5), .LatAlpha (0.05), .LatSearchWin ([0 400])
% Output: struct L (N-by-1 unless noted)
%   .LatRaw, .LatCorr   value onset latency (ms), NaN if none / not computed
%   .SignRaw, .SignCorr sign of response-window AUC - 0.5
%   .AUCTCRaw, .AUCTCCorr  N-by-T per-bin AUC (causal smoothing)
%   .VisLat             visual response latency (ms)
%   .VisTC              N-by-T mean (rate - baseline), normalized to its peak
%                       in LatSearchWin (NaN if the peak is <= 0)
%   .t, .Params
%
% Xuefei Yu, 2026-11-06

P = Params;
if ~isfield(P,'RespWin'),       P.RespWin = [100 200]; end
if ~isfield(P,'BaseWin'),       P.BaseWin = [-300 0];  end
if ~isfield(P,'LatSmoothMs'),   P.LatSmoothMs = 30;    end
if ~isfield(P,'LatConsecBins'), P.LatConsecBins = 5;   end
if ~isfield(P,'LatAlpha'),      P.LatAlpha = 0.05;     end
if ~isfield(P,'LatSearchWin'),  P.LatSearchWin = [0 400]; end

t = t(:)'; nT = numel(t); N = numel(Trials);
nSm = max(1, round(P.LatSmoothMs / median(diff(t))));
iB = t >= P.BaseWin(1) & t <= P.BaseWin(2);
iR = t >= P.RespWin(1) & t <= P.RespWin(2);
kS = find(t >= P.LatSearchWin(1) & t <= P.LatSearchWin(2));
tS = t(kS);

L.LatRaw  = nan(N,1); L.LatCorr  = nan(N,1);
L.SignRaw = nan(N,1); L.SignCorr = nan(N,1);
L.AUCTCRaw = nan(N,nT); L.AUCTCCorr = nan(N,nT);
L.VisLat = nan(N,1); L.VisTC = nan(N,nT);

for n = find(mask(:))'
    X = Trials(n).Rate; hi = logical(Trials(n).IsHigh(:));
    if isempty(X) || sum(hi) < 3 || sum(~hi) < 3, continue; end
    base = mean(X(:,iB), 2, 'omitnan');
    resp = mean(X(:,iR), 2, 'omitnan');
    L.SignRaw(n)  = sign(MannWhitneyAUC(resp(hi), resp(~hi)) - 0.5);
    L.SignCorr(n) = sign(MannWhitneyAUC(resp(hi) - base(hi), resp(~hi) - base(~hi)) - 0.5);

    Xs = movmean(X, [nSm-1 0], 2, 'omitnan');     % causal
    Xc = Xs - base;
    pR = nan(1,nT); pC = nan(1,nT);
    for k = 1:nT
        [L.AUCTCRaw(n,k),  pR(k)] = MWAucP(Xs(hi,k), Xs(~hi,k));
        [L.AUCTCCorr(n,k), pC(k)] = MWAucP(Xc(hi,k), Xc(~hi,k));
    end
    sigR = pR(kS) < P.LatAlpha & sign(L.AUCTCRaw(n,kS)  - 0.5) == L.SignRaw(n);
    sigC = pC(kS) < P.LatAlpha & sign(L.AUCTCCorr(n,kS) - 0.5) == L.SignCorr(n);
    L.LatRaw(n)  = FirstRunTime(sigR, tS, P.LatConsecBins);
    L.LatCorr(n) = FirstRunTime(sigC, tS, P.LatConsecBins);

    % visual response: rate above the trial's own baseline, all RF trials
    pv = nan(1, numel(kS)); mv = nan(1, numel(kS));
    for j = 1:numel(kS)
        d = Xc(:, kS(j)); d = d(~isnan(d));
        if numel(d) >= 5 && any(d ~= 0), pv(j) = signrank(d); mv(j) = mean(d); end
    end
    L.VisLat(n) = FirstRunTime(pv < P.LatAlpha & mv > 0, tS, P.LatConsecBins);
    m = mean(Xc, 1, 'omitnan'); pk = max(m(kS));
    if pk > 0, L.VisTC(n,:) = m / pk; end
end
L.t = t; L.Params = P;
end

% ---------------------------------------------------------------------------
function lat = FirstRunTime(sig, tS, K)
% first time starting a run of >= K consecutive true bins
lat = NaN; run = 0;
for j = 1:numel(sig)
    if sig(j), run = run + 1; else, run = 0; end
    if run >= K, lat = tS(j - K + 1); return; end
end
end

function [auc, p] = MWAucP(x, y)
% Mann-Whitney AUC and two-sided normal-approx p with tie correction
% (same as ranksum 'approximate', without the per-call overhead).
x = x(~isnan(x)); y = y(~isnan(y));
nx = numel(x); ny = numel(y); auc = NaN; p = NaN;
if nx < 3 || ny < 3, return; end
v = [x(:); y(:)]; r = tiedrank(v); n = nx + ny;
U = sum(r(1:nx)) - nx*(nx+1)/2;
auc = U / (nx*ny);
[~,~,g] = unique(v); tc = accumarray(g, 1);
s2 = nx*ny/12 * ((n+1) - sum(tc.^3 - tc) / (n*(n-1)));
if s2 <= 0, p = 1; return; end
z = (U - nx*ny/2 - 0.5*sign(U - nx*ny/2)) / sqrt(s2);   % continuity correction
p = 2 * normcdf(-abs(z));
end
