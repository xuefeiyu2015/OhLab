function P = ComputePreTargetAUC(aucPsthObj, tObj, aucPsthLoc, tLoc, aucObj, aucLoc, win, groupMasks, groupNames)
% ComputePreTargetAUC  Mean AUC in a pre-target window per task and group.
% Each task is sign-aligned by its OWN 100-200 ms AUC, so a pre-target
% deviation from 0.5 means pre-target activity already differs in the
% direction of that task's later value preference. This is expected under the
% OneDR block-structure (anticipation/drift) confound, and should not happen in
% the interleaved object task.
%
% Inputs
%   aucPsthObj/aucPsthLoc — N-by-T AUC time courses; tObj/tLoc — 1-by-T (ms)
%   aucObj/aucLoc         — N-vectors, scalar 100-200 ms AUCs
%   win                   — [t1 t2] pre-target window (ms), e.g. [-300 0]
%   groupMasks/groupNames — cell{G} masks and names
% Output: table, rows = groups, columns per task (Obj*, Loc*):
%   n, ObjMean, ObjSEM, ObjP, LocMean, LocSEM, LocP   (signrank vs 0.5)
%   P.Properties.UserData.PerNeuron = struct('Obj',N-by-1,'Loc',N-by-1)
%
% Xuefei Yu, 2026-11-05

wObj = WindowMeanAligned(aucPsthObj, tObj, aucObj, win);
wLoc = WindowMeanAligned(aucPsthLoc, tLoc, aucLoc, win);

G = numel(groupMasks);
[n, mO, sO, pO, mL, sL, pL] = deal(nan(G,1));
for g = 1:G
    m = groupMasks{g}(:);
    n(g) = sum(m);
    [mO(g), sO(g), pO(g)] = Summ(wObj(m));
    [mL(g), sL(g), pL(g)] = Summ(wLoc(m));
end
P = table(groupNames(:), n, mO, sO, pO, mL, sL, pL, 'VariableNames', ...
    {'Group','n','ObjMean','ObjSEM','ObjP','LocMean','LocSEM','LocP'});
P.Properties.UserData.PerNeuron = struct('Obj', wObj, 'Loc', wLoc);
P.Properties.UserData.Window = win;
end

function w = WindowMeanAligned(A, t, aucScalar, win)
idx = t >= win(1) & t <= win(2);
w = mean(A(:,idx), 2, 'omitnan');
flip = aucScalar(:) < 0.5;
w(flip) = 1 - w(flip);
end

function [mu, se, p] = Summ(v)
v = v(~isnan(v));
mu = mean(v); se = std(v) / sqrt(numel(v)); p = NaN;
if numel(v) >= 3, p = signrank(v - 0.5); end
if isempty(v), mu = NaN; se = NaN; end
end
