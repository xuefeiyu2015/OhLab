function K = ComputeSlideDeckComparison(DA, DB, names, Params)
% ComputeSlideDeckComparison  FEF vs SC tests for the slide deck, and the
% slide headlines (worded from the results, never assumed).
%
% Inputs
%   DA, DB — ComputeSlideDeckData outputs for area A (FEF) and B (SC)
%   names  — {'FEF','SC'}
%   Params — .Alpha
% Output: struct K
%   .Table     Measure, A, B, Test, p, pHolm (Holm over all rows)
%   .ORDots    struct for PlotPresentationEffectDots (OR per condition)
%   .FracBars  struct for PlotPresentationFractionBars
%   .Headline  struct with fields Overlap, Split, Corrected, AreaOR, AreaFrac
%   .Subtitle  struct with the same fields
%
% Xuefei Yu, 2026-11-05

a = Params.Alpha;
cond   = {'All', 'WithBase', 'WithoutBase', 'Corrected'};
condLb = {'All neurons', 'With anticipation', 'Without anticipation', 'Baseline-corrected'};

rows = {};
for i = 1:numel(DA.Frac.Labels)
    kA = DA.Frac.k(i); nA = DA.Frac.n(i); kB = DB.Frac.k(i); nB = DB.Frac.n(i);
    [~, p] = fishertest([kA nA-kA; kB nB-kB]);
    rows(end+1,:) = {['% ' DA.Frac.Labels{i}], kA/nA, kB/nB, 'Fisher', p}; %#ok<AGROW>
end
[est, lo, hi] = deal(nan(numel(cond), 2));
for c = 1:numel(cond)
    OA = DA.Overlap.(cond{c}); OB = DB.Overlap.(cond{c});
    z = (OA.LogOR - OB.LogOR) / sqrt(OA.LogORSE^2 + OB.LogORSE^2);
    rows(end+1,:) = {['Overlap OR: ' condLb{c}], exp(OA.LogOR), exp(OB.LogOR), 'z (log OR)', 2*normcdf(-abs(z))}; %#ok<AGROW>
    Os = {OA, OB};
    for j = 1:2
        est(c,j) = exp(Os{j}.LogOR);
        lo(c,j)  = exp(Os{j}.LogOR - 1.96*Os{j}.LogORSE);
        hi(c,j)  = exp(Os{j}.LogOR + 1.96*Os{j}.LogORSE);
    end
end
Tb = cell2table(rows, 'VariableNames', {'Measure', names{1}, names{2}, 'Test', 'p'});
Tb.pHolm = HolmAdjust(Tb.p);
K.Table = Tb;
nF = numel(DA.Frac.Labels);

%% Headlines (data-driven)
vA = cellfun(@(c) DA.Overlap.(c).Verdict, cond, 'UniformOutput', false);
vB = cellfun(@(c) DB.Overlap.(c).Verdict, cond, 'UniformOutput', false);
K.Headline.Overlap   = VerdictHeadline({vA{1}, vB{1}});
K.Subtitle.Overlap   = 'Bars = observed neurons; black line = expected if object and location coding were independent; gray = 95% chance range';
vs = [vA(2:3), vB(2:3)];
if all(strcmp(vs, vs{1}))
    K.Headline.Split = ['With or without location anticipation: ' lower(VerdictPhrase(vs{1}))];
else
    K.Headline.Split = 'Overlap differs with vs without location anticipation';
end
K.Subtitle.Split = sprintf(['Anticipation = high vs low location difference before target onset (%d-%d ms, ranksum p < %.2f). ' ...
    'Black line = chance'], Params.BaseWin, a);
K.Headline.Corrected = ['After removing the pre-target baseline: ' lower(VerdictPhrase(vA{4}, vB{4}))];
K.Subtitle.Corrected = sprintf('Response %d-%d ms minus the same trial''s %d-%d ms baseline (both tasks). Black line = chance', ...
    Params.RespWin, Params.BaseWin);

pOR = Tb.pHolm(nF+1:end);
if all(pOR >= a)
    K.Headline.AreaOR = sprintf('%s and %s show the same overlap pattern in every analysis', names{:});
else
    K.Headline.AreaOR = sprintf('%s and %s differ in overlap for some analyses', names{:});
end
K.Subtitle.AreaOR = 'Odds ratio of coding both value types (95% CI); 1 = independent. p = Holm-corrected area difference';
pF = Tb.pHolm(1:nF);
sigF = find(pF < a);
if isempty(sigF)
    K.Headline.AreaFrac = sprintf('%s and %s have similar proportions of value-coding neurons', names{:});
elseif numel(sigF) == 1
    i = sigF;
    K.Headline.AreaFrac = sprintf('%s and %s differ only in %s (%s %.0f%% vs %s %.0f%%)', names{:}, ...
        lower(DA.Frac.Labels{i}), names{2}, 100*DB.Frac.Frac(i), names{1}, 100*DA.Frac.Frac(i));
else
    K.Headline.AreaFrac = sprintf('%s and %s differ in: %s', names{:}, strjoin(lower(DA.Frac.Labels(sigF)), ', '));
end
K.Subtitle.AreaFrac = sprintf('%% RF neurons (100-200 ms high vs low, p < %.2f; anticipation = pre-target difference). Dashed = chance', a);

%% Plot structs
K.ORDots = struct('Est', est, 'Lo', lo, 'Hi', hi, 'RowLabels', {condLb}, ...
    'Refs', struct('Value', {1}, 'Label', {'independent'}), 'XLabel', 'Odds ratio of coding both (log scale)', ...
    'XLim', [min(0.25, min(lo(:))*0.8), max(hi(:))*1.4], 'LogX', true, 'pArea', pOR, ...
    'LegendPos', [0.84 0.70 0.14 0.12], 'Headline', K.Headline.AreaOR, 'Subtitle', K.Subtitle.AreaOR);
K.FracBars = struct('Frac', [DA.Frac.Frac, DB.Frac.Frac], 'k', [DA.Frac.k, DB.Frac.k], 'n', [DA.Frac.n, DB.Frac.n], ...
    'pArea', pF', 'Categories', {DA.Frac.Labels}, 'Chance', a, ...
    'Headline', K.Headline.AreaFrac, 'Subtitle', K.Subtitle.AreaFrac);
end

% ---------------------------------------------------------------------------
function h = VerdictHeadline(v)
if all(strcmp(v, 'chance'))
    h = 'Overlap matches chance: object and location value are coded independently';
elseif all(strcmp(v, 'above'))
    h = 'More neurons code both values than expected by chance';
elseif all(strcmp(v, 'below'))
    h = 'Fewer neurons code both values than chance: separate populations';
else
    h = 'Overlap relative to chance differs between areas';
end
end

function s = VerdictPhrase(varargin)
v = varargin;
if all(strcmp(v, 'chance')),    s = 'Overlap stays at chance';
elseif all(strcmp(v, 'above')), s = 'Overlap stays above chance';
elseif all(strcmp(v, 'below')), s = 'Overlap stays below chance';
else,                           s = 'Overlap differs between areas';
end
end

function q = HolmAdjust(p)
[ps, o] = sort(p(:)); m = numel(ps);
adj = min(1, cummax((m - (1:m)' + 1) .* ps));
q = nan(size(p)); q(o) = adj;
end
