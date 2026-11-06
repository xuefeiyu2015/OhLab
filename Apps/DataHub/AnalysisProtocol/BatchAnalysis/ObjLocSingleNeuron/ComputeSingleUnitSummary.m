function S = ComputeSingleUnitSummary(T, CObj, CLoc, mask, Params)
% ComputeSingleUnitSummary  All single-unit object vs location value metrics
% for one area / subpopulation.
%
% Inputs
%   T          — struct from LoadObjLocNeuronTable (LoadTrials = true)
%   CObj, CLoc — structs from ComputeTrialWindowAUC
%   mask       — N-by-1 logical neurons to include
%   Params     — .Alpha, .RespWin, .BaseWin, .nSplit, .nBoot, .Seed
% Output: struct S
%   .n
%   .Coding     table: Measure, k, n, Frac  (obj/loc x raw/corrected/baseline)
%   .PreTgt     table: per task, baseline AUC aligned by own response
%               preference, d = (AUCbase-0.5)*sign(AUCresp-0.5): mean, SEM, p
%   .PreTgtD    struct .Obj/.Loc per-neuron d (for area comparison)
%   .BaseRespAgree  table (loc baseline-diff neurons: baseline sign == response sign)
%   .Overlap    ComputeOverlapStats table, rows {'Raw','Corrected'}
%   .Concord    ComputeSignConcordance tables .Raw/.Corrected (n, nConcordant)
%   .Split      ComputeBaselineSplitResults output
%   .Rel        split-half noise ceiling, .Corrected / .Raw
%
% Xuefei Yu, 2026-11-05

a = Params.Alpha; m = mask(:);
S.n = sum(m);

% --- coding fractions
meas = {'Obj response (raw)'; 'Obj response (corrected)'; 'Loc response (raw)'; ...
        'Loc response (corrected)'; 'Obj baseline diff'; 'Loc baseline diff'};
pv   = {CObj.pRaw, CObj.pCorr, CLoc.pRaw, CLoc.pCorr, CObj.pBase, CLoc.pBase};
k = zeros(6,1); nn = zeros(6,1);
for i = 1:6
    valid = m & ~isnan(pv{i});
    k(i) = sum(valid & pv{i} < a); nn(i) = sum(valid);
end
S.Coding = table(meas, k, nn, k./nn, 'VariableNames', {'Measure','k','n','Frac'});

% --- pre-target (baseline) AUC aligned by own response preference
dO = (CObj.AUCBase(m) - 0.5) .* sign(CObj.AUCRaw(m) - 0.5);
dL = (CLoc.AUCBase(m) - 0.5) .* sign(CLoc.AUCRaw(m) - 0.5);
dO = dO(~isnan(dO)); dL = dL(~isnan(dL));
S.PreTgtD = struct('Obj', dO, 'Loc', dL);
S.PreTgt = table({'Obj'; 'Loc'}, [numel(dO); numel(dL)], [mean(dO); mean(dL)], ...
    [std(dO)/sqrt(numel(dO)); std(dL)/sqrt(numel(dL))], [signrank(dO); signrank(dL)], ...
    'VariableNames', {'Task','n','MeanAlignedBaselineAUCminus05','SEM','pSignrank'});

% --- does the oneDR baseline preference match the response preference?
S.BaseRespAgree = ComputeSignConcordance(CLoc.AUCBase, CLoc.AUCRaw, ...
    {m & CLoc.pBase < a, m & CLoc.pBase >= a}, {'LocBaselineDiff', 'NoLocBaselineDiff'});

% --- object vs location relationship
S.Overlap = [ ...
    ComputeOverlapStats(CObj.AUCRaw,  CObj.pRaw,  CLoc.AUCRaw,  CLoc.pRaw,  {m}, {'Raw'},       a); ...
    ComputeOverlapStats(CObj.AUCCorr, CObj.pCorr, CLoc.AUCCorr, CLoc.pCorr, {m}, {'Corrected'}, a)];
S.Concord.Raw       = ComputeSignConcordance(CObj.AUCRaw,  CLoc.AUCRaw,  {m}, {'Raw'});
S.Concord.Corrected = ComputeSignConcordance(CObj.AUCCorr, CLoc.AUCCorr, {m}, {'Corrected'});

% --- baseline split
S.Split = ComputeBaselineSplitResults(CObj, CLoc, m, Params);

% --- noise ceiling
P = Params; P.Variant = 'corrected';
S.Rel.Corrected = ComputeSplitHalfReliability(T.Obj.Trials, T.Obj.TrialTime, T.Loc.Trials, T.Loc.TrialTime, m, P);
P.Variant = 'raw';
S.Rel.Raw       = ComputeSplitHalfReliability(T.Obj.Trials, T.Obj.TrialTime, T.Loc.Trials, T.Loc.TrialTime, m, P);
end
