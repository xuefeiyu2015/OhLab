function PD = ComputePresentationData(SA, SB, Cmp, CA, CB, names, Params)
% ComputePresentationData  Everything the presentation figures need, already
% computed (plot functions only render). Headlines are worded from the
% results, so they stay truthful if the data change.
%
% Inputs
%   SA, SB  — ComputeSingleUnitSummary outputs (area A = FEF, B = SC)
%   Cmp     — ComputeAreaComparison table
%   CA, CB  — struct('Obj',CObj,'Loc',CLoc,'Mask',mask) per area (for scatter)
%   names   — {'FEF','SC'}
%   Params  — .Alpha
% Output: struct PD with fields .Fig1 ... .Fig5 (data + .Headline/.Subtitle)
%
% Xuefei Yu, 2026-11-05

a = Params.Alpha;
pc = @(m) Cmp.pHolm(strcmp(Cmp.Measure, m));
S = {SA, SB};

%% Fig 1 — coding fractions (baseline-corrected target response)
F = zeros(2,2); K = zeros(2,2); N = zeros(2,2);   % rows: Obj, Loc; cols: areas
rowsC = {'Obj response (corrected)', 'Loc response (corrected)'};
for i = 1:2
    for j = 1:2
        r = strcmp(S{j}.Coding.Measure, rowsC{i});
        K(i,j) = S{j}.Coding.k(r); N(i,j) = S{j}.Coding.n(r); F(i,j) = K(i,j)/N(i,j);
    end
end
PD.Fig1.Frac = F; PD.Fig1.k = K; PD.Fig1.n = N;
PD.Fig1.pArea = [pc(['Frac ' rowsC{1}]), pc(['Frac ' rowsC{2}])];
PD.Fig1.Categories = {'Object value', 'Location value'};
PD.Fig1.Chance = a;
PD.Fig1.Headline = 'Both FEF and SC encode object and location value';
PD.Fig1.Subtitle = sprintf(['%% RF neurons with high vs low value difference, target response ' ...
    '(%d-%d ms, baseline-subtracted); dashed = chance (%.0f%%)'], Params.RespWin, 100*a);

%% Fig 2 — pre-target (baseline) value difference
rowsB = {'Obj baseline diff', 'Loc baseline diff'};
for i = 1:2
    for j = 1:2
        r = strcmp(S{j}.Coding.Measure, rowsB{i});
        K(i,j) = S{j}.Coding.k(r); N(i,j) = S{j}.Coding.n(r); F(i,j) = K(i,j)/N(i,j);
    end
end
PD.Fig2.Frac = F; PD.Fig2.k = K; PD.Fig2.n = N;
PD.Fig2.pArea = [pc(['Frac ' rowsB{1}]), pc(['Frac ' rowsB{2}])];
PD.Fig2.Categories = {'Object task', 'Location task'};
PD.Fig2.Chance = a;
agree = cellfun(@(s) s.BaseRespAgree.Fraction(1), S);
locHigher = all(F(2,:) > 2*F(1,:));
if locHigher
    PD.Fig2.Headline = 'Location value appears before the target; object value does not';
else
    PD.Fig2.Headline = 'Pre-target value differences in both tasks';
end
PD.Fig2.Subtitle = sprintf(['%% neurons with high vs low difference before target onset (%d-%d ms). ' ...
    'Same sign as target response: %s %.0f%%, %s %.0f%%'], Params.BaseWin, names{1}, 100*agree(1), names{2}, 100*agree(2));

%% Fig 3 — overlap odds ratio (corrected)
[est, lo, hi] = deal(zeros(2,1));
for j = 1:2
    row = S{j}.Overlap(2,:);
    [l, s] = LogOR(row);
    est(j) = exp(l); lo(j) = exp(l - 1.96*s); hi(j) = exp(l + 1.96*s);
    PD.Fig3.Obs(j) = row.nBoth; PD.Fig3.Exp(j) = row.ExpBoth; PD.Fig3.n(j) = row.n;
end
PD.Fig3.Est = est; PD.Fig3.Lo = lo; PD.Fig3.Hi = hi;
PD.Fig3.pArea = pc('Overlap log OR (Corrected)');
if all(lo > 1)
    PD.Fig3.Headline = 'Neurons coding one value type are more likely to code the other';
elseif all(hi < 1)
    PD.Fig3.Headline = 'Object and location value use separate neurons';
else
    PD.Fig3.Headline = 'Overlap of object- and location-coding neurons';
end
PD.Fig3.Subtitle = sprintf('Odds ratio of coding both (95%% CI); 1 = independent, <1 = separate populations. Both: %s %d vs %.0f expected, %s %d vs %.0f', ...
    names{1}, PD.Fig3.Obs(1), PD.Fig3.Exp(1), names{2}, PD.Fig3.Obs(2), PD.Fig3.Exp(2));
PD.Fig3.E = struct('Est', est', 'Lo', lo', 'Hi', hi', 'RowLabels', {{'Coding both value types'}}, ...
    'Refs', struct('Value', {1}, 'Label', {'independent'}), 'XLabel', 'Odds ratio (log scale)', ...
    'XLim', [min(0.5, min(lo)*0.7), max(hi)*1.6], 'LogX', true, 'pArea', PD.Fig3.pArea, 'LegendPos', [0.84 0.15 0.14 0.18], ...
    'Headline', PD.Fig3.Headline, 'Subtitle', PD.Fig3.Subtitle);

%% Fig 4 — AUC scatter per area (corrected)
Cs = {CA, CB};
for j = 1:2
    m = Cs{j}.Mask(:);
    PD.Fig4.AUCObj{j} = Cs{j}.Obj.AUCCorr(m);
    PD.Fig4.AUCLoc{j} = Cs{j}.Loc.AUCCorr(m);
    PD.Fig4.Both{j}   = Cs{j}.Obj.pCorr(m) < a & Cs{j}.Loc.pCorr(m) < a;
    PD.Fig4.Either{j} = xor(Cs{j}.Obj.pCorr(m) < a, Cs{j}.Loc.pCorr(m) < a);
    PD.Fig4.rObs(j)   = S{j}.Rel.Corrected.rObs;
    PD.Fig4.rDis(j)   = S{j}.Rel.Corrected.rDis;
end
PD.Fig4.Headline = 'Each neuron''s object and location value preferences are only loosely related';
PD.Fig4.Subtitle = sprintf('Each dot = one RF neuron; AUC > 0.5 = prefers high value (%d-%d ms, baseline-subtracted). Dotted = identical tuning', Params.RespWin);

%% Fig 5 — observed vs noise-corrected correlation
for j = 1:2
    R = S{j}.Rel.Corrected;
    PD.Fig5.rObs(j) = R.rObs; PD.Fig5.rObsCI(j,:) = prctile(R.rObsBoot, [2.5 97.5]);
    PD.Fig5.rDis(j) = R.rDis; PD.Fig5.rDisCI(j,:) = R.rDisCI;
    PD.Fig5.relObj(j) = R.relObj; PD.Fig5.relLoc(j) = R.relLoc;
end
PD.Fig5.pArea = pc('Noise-corrected r (Corrected)');
if all(PD.Fig5.rDisCI(:,2) < 1) && all(PD.Fig5.rDisCI(:,1) > 0)
    PD.Fig5.Headline = 'Tuning is shared but far from identical, even after noise correction';
elseif all(PD.Fig5.rDisCI(:,2) < 1)
    PD.Fig5.Headline = 'Tuning is far from identical, even after noise correction';
else
    PD.Fig5.Headline = 'Once noise is accounted for, tuning could be identical';
end
PD.Fig5.Subtitle = sprintf(['Pearson r of object vs location AUC across neurons (95%% bootstrap CI). ' ...
    'Noise-corrected = r / sqrt(split-half reliability_{obj} x reliability_{loc}); 1 = identical tuning']);
PD.Fig5.E = struct('Est', [PD.Fig5.rObs; PD.Fig5.rDis], ...
    'Lo', [PD.Fig5.rObsCI(:,1)'; PD.Fig5.rDisCI(:,1)'], 'Hi', [PD.Fig5.rObsCI(:,2)'; PD.Fig5.rDisCI(:,2)'], ...
    'RowLabels', {{'Observed r', 'Noise-corrected r'}}, ...
    'Refs', struct('Value', {0, 1}, 'Label', {'independent', 'identical tuning'}), ...
    'XLabel', 'Correlation of object vs location AUC across neurons', ...
    'XLim', [min(-0.2, min(PD.Fig5.rObsCI(:))-0.05), max(1.2, max(PD.Fig5.rDisCI(:))+0.05)], 'LogX', false, ...
    'pArea', [NaN; PD.Fig5.pArea], 'LegendPos', [0.84 0.40 0.14 0.16], 'Headline', PD.Fig5.Headline, 'Subtitle', PD.Fig5.Subtitle);
PD.Names = names;
PD.AreaDiffSummary = Cmp(Cmp.pHolm < a, :);
end

function [l, s] = LogOR(row)
x = [row.nBoth, row.nObj - row.nBoth, row.nLoc - row.nBoth, row.n - row.nObj - row.nLoc + row.nBoth];
if any(x == 0), x = x + 0.5; end
l = log(x(1)*x(4) / (x(2)*x(3)));
s = sqrt(sum(1 ./ x));
end
