function P = ComputePopulationOnset(TC, t, Params)
% ComputePopulationOnset  Population signal onset = half-max time of the mean
% time course, with a bootstrap over neurons.
%
% Half-max: b = mean of the curve in LatBaseWin; pk = max within
% LatSearchWin; onset = first time (linearly interpolated) in LatSearchWin,
% before the peak, where the curve reaches b + 0.5*(pk - b).
%
% Inputs
%   TC     — n-by-T per-neuron time courses, already aligned so that the
%            signal is positive (e.g. sign-aligned AUC - 0.5, normalized FR)
%   t      — 1-by-T time (ms)
%   Params — .LatBaseWin ([-100 0]), .LatSearchWin ([0 400]), .nBoot, .Seed
% Output: struct P
%   .n, .Mean (1-by-T), .SEM, .Onset (ms), .CI (1x2, 95% percentile),
%   .Boot (nBoot-by-1 bootstrap onsets; same seed -> same resamples)
%
% Xuefei Yu, 2026-11-06

if ~isfield(Params,'LatBaseWin'),   Params.LatBaseWin = [-100 0];  end
if ~isfield(Params,'LatSearchWin'), Params.LatSearchWin = [0 400]; end
TC = TC(~all(isnan(TC), 2), :);
t = t(:)'; P.n = size(TC, 1);
P.Mean = mean(TC, 1, 'omitnan');
P.SEM  = std(TC, 0, 1, 'omitnan') / sqrt(max(P.n, 1));
iB = t >= Params.LatBaseWin(1) & t <= Params.LatBaseWin(2);
kS = find(t >= Params.LatSearchWin(1) & t <= Params.LatSearchWin(2));
P.Onset = HalfMax(P.Mean, t, iB, kS);
P.Boot = nan(Params.nBoot, 1); P.CI = [NaN NaN];
if P.n < 3, return; end
s = RandStream('mt19937ar', 'Seed', Params.Seed);
for b = 1:Params.nBoot
    idx = randi(s, P.n, P.n, 1);
    P.Boot(b) = HalfMax(mean(TC(idx,:), 1, 'omitnan'), t, iB, kS);
end
P.CI = prctile(P.Boot, [2.5 97.5]);
end

function on = HalfMax(c, t, iB, kS)
on = NaN;
b = mean(c(iB), 'omitnan');
[pk, ip] = max(c(kS)); ip = kS(ip);
if ~(pk > b), return; end
thr = b + 0.5 * (pk - b);
for k = kS(1):ip
    if c(k) >= thr
        if k == kS(1) || isnan(c(k-1)), on = t(k);
        else, on = t(k-1) + (thr - c(k-1)) / (c(k) - c(k-1)) * (t(k) - t(k-1));
        end
        return;
    end
end
end
