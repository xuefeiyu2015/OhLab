function K = ComputeLatencyAreaComparison(AA, AB, names, Params)
% ComputeLatencyAreaComparison  FEF vs SC timing tests and data-driven headlines.
%
% Population onset: difference of half-max onsets (B - A), 95% CI and
%   two-sided p from the independent bootstrap distributions
%   ((count+1)/(nBoot+1), so never exactly 0).
% Single-neuron latency: medians, ranksum.
% Holm correction within each signal family (value corrected / value raw /
% visual / value - visual; object + location, population + single neuron).
%
% Inputs
%   AA, AB — ComputeAreaLatencyData outputs for area A (FEF) and B (SC)
%   names  — {'FEF','SC'}; Params — .Alpha
% Output: struct K
%   .Table   Task, Signal, Measure, A, B, Diff (B-A, ms), DiffLo, DiffHi,
%            nA, nB, Test, p, pHolm
%   .Headline / .Subtitle  struct: Value, Control, Visual, Relative
%
% Xuefei Yu, 2026-11-06

tasks = {'Obj', 'Loc'}; taskLb = {'Object', 'Location'};
sigs  = {'Corrected', 'RawNoAntic', 'Visual', 'ValueMinusVisual'};
sigLb = {'Value (baseline-corrected)', 'Value (raw, no anticipation)', 'Visual response', 'Value - visual'};
rows = {};
for i = 1:2
    for s = 1:numel(sigs)
        SA = AA.(tasks{i}).(sigs{s}); SB = AB.(tasks{i}).(sigs{s});
        if isfield(SA, 'Pop')
            d = SB.Pop.Boot - SA.Pop.Boot; d = d(~isnan(d));
            if isempty(d), ci = [NaN NaN]; p = NaN;
            else, ci = prctile(d, [2.5 97.5]); nb = numel(d);
                p = min(1, 2*min(sum(d <= 0) + 1, sum(d >= 0) + 1) / (nb + 1));
            end
            rows(end+1,:) = {taskLb{i}, sigLb{s}, 'Population half-max onset', SA.Pop.Onset, SB.Pop.Onset, ...
                SB.Pop.Onset - SA.Pop.Onset, ci(1), ci(2), SA.Pop.n, SB.Pop.n, 'bootstrap', p}; %#ok<AGROW>
        end
        la = SA.Lat(~isnan(SA.Lat)); lb = SB.Lat(~isnan(SB.Lat));
        p = NaN; if numel(la) >= 3 && numel(lb) >= 3, p = ranksum(la, lb); end
        rows(end+1,:) = {taskLb{i}, sigLb{s}, 'Single-neuron latency (median)', median(la), median(lb), ...
            median(lb) - median(la), NaN, NaN, numel(la), numel(lb), 'ranksum', p}; %#ok<AGROW>
    end
end
Tb = cell2table(rows, 'VariableNames', {'Task', 'Signal', 'Measure', names{1}, names{2}, ...
    'Diff', 'DiffLo', 'DiffHi', ['n' names{1}], ['n' names{2}], 'Test', 'p'});
Tb.pHolm = nan(height(Tb), 1);
for s = 1:numel(sigLb)
    r = strcmp(Tb.Signal, sigLb{s});
    Tb.pHolm(r) = HolmAdjust(Tb.p(r));
end
K.Table = Tb;

a = Params.Alpha;
K.Headline.Value   = TimingHeadline(Tb, 'Value (baseline-corrected)', names, a, 'value signal');
K.Headline.Control = TimingHeadline(Tb, 'Value (raw, no anticipation)', names, a, 'value signal (no-anticipation neurons)');
K.Headline.Visual  = TimingHeadline(Tb, 'Visual response', names, a, 'visual response');
K.Headline.Relative = TimingHeadline(Tb, 'Value - visual', names, a, 'value lag behind the visual response');
K.Subtitle.Value = sprintf(['Top: sign-aligned AUC (high vs low value), response minus each trial''s %d to %d ms baseline; ' ...
    'dot + bar = population half-max onset [95%% CI]. Bottom: single-neuron onset (first %d consecutive 10 ms bins p < %.2f, ' ...
    'neurons with an onset). p = Holm-corrected FEF vs SC'], Params.BaseWin, Params.LatConsecBins, Params.LatAlpha);
K.Subtitle.Control = 'Same as the main timing figure, raw responses, only neurons without a pre-target value difference';
K.Subtitle.Visual  = 'Top: visual response (rate minus trial baseline, normalized to peak), all selected neurons. Bottom: single-neuron visual latency';
K.Subtitle.Relative = 'Per neuron: value onset minus that neuron''s visual onset (same task, baseline-corrected value latency)';
end

% ---------------------------------------------------------------------------
function h = TimingHeadline(Tb, sig, names, a, what)
% Phrase from the population onset when available, else the single-neuron median.
parts = {};
for tk = {'Object', 'Location'}
    r = strcmp(Tb.Task, tk{1}) & strcmp(Tb.Signal, sig);
    rp = find(r & strcmp(Tb.Measure, 'Population half-max onset'), 1);
    if isempty(rp), rp = find(r, 1); end
    d = Tb.Diff(rp); p = Tb.pHolm(rp);
    if isnan(d)
        parts{end+1} = sprintf('%s: n/a', tk{1}); %#ok<AGROW>
    elseif p < a
        first = names{1 + (d < 0)};                     % d = B - A; d < 0 -> B earlier / shorter
        if contains(what, 'lag')
            parts{end+1} = sprintf('%s: %s shorter by %.0f ms', tk{1}, first, abs(d)); %#ok<AGROW>
        else
            parts{end+1} = sprintf('%s: %s first by %.0f ms', tk{1}, first, abs(d)); %#ok<AGROW>
        end
    else
        parts{end+1} = sprintf('%s: n.s. (%s %+.0f ms)', tk{1}, names{2}, d); %#ok<AGROW>
    end
end
if contains(what, 'lag'), lb = ''; else, lb = ' onset'; end
h = sprintf('%s%s%s — %s', upper(what(1)), what(2:end), lb, strjoin(parts, '; '));
end

function q = HolmAdjust(p)
q = nan(size(p)); ok = ~isnan(p);
[ps, o] = sort(p(ok)); m = numel(ps);
adj = min(1, cummax((m - (1:m)' + 1) .* ps(:)));
tmp = nan(m,1); tmp(o) = adj; q(ok) = tmp;
end
