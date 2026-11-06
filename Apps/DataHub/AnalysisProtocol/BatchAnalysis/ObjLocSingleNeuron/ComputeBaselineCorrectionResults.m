function B = ComputeBaselineCorrectionResults(T, CObj, CLoc, mask, Params)
% ComputeBaselineCorrectionResults  Compare object vs location value coding
% before and after trial-wise baseline subtraction, for one subpopulation.
%
% Inputs
%   T          — struct from LoadObjLocNeuronTable (LoadTrials = true)
%   CObj, CLoc — structs from ComputeTrialWindowAUC
%   mask       — N-by-1 logical subpopulation
%   Params     — as for ComputeObjLocSingleNeuronResults (+ Alpha)
% Output: struct B
%   .Raw, .Corrected — ComputeObjLocSingleNeuronResults output for each variant
%   .Overlap         — ComputeOverlapStats table, rows {'Raw','Corrected'}
%   .BaselineSig     — table: # neurons with a baseline difference per task
%   .Validation      — Spearman rho / median |diff| of trial-level raw AUC vs
%                      the protocol's saved TgtAUC_RF (smoothed FR), per task
%
% Xuefei Yu, 2026-11-05

a = Params.Alpha; m = mask(:);
Tr = ComputeAUCVariantTable(T, CObj, CLoc, 'raw');
Tc = ComputeAUCVariantTable(T, CObj, CLoc, 'corrected');
B.Raw       = ComputeObjLocSingleNeuronResults(Tr, m, Params);
B.Corrected = ComputeObjLocSingleNeuronResults(Tc, m, Params);

B.Overlap = [ ...
    ComputeOverlapStats(CObj.AUCRaw,  CObj.pRaw,  CLoc.AUCRaw,  CLoc.pRaw,  {m}, {'Raw'},       a); ...
    ComputeOverlapStats(CObj.AUCCorr, CObj.pCorr, CLoc.AUCCorr, CLoc.pCorr, {m}, {'Corrected'}, a)];

task = {'Obj'; 'Loc'};
nBase = [sum(m & CObj.pBase < a); sum(m & CLoc.pBase < a)];
nRaw  = [sum(m & CObj.pRaw  < a); sum(m & CLoc.pRaw  < a)];
nCorr = [sum(m & CObj.pCorr < a); sum(m & CLoc.pCorr < a)];
B.BaselineSig = table(task, repmat(sum(m),2,1), nBase, nRaw, nCorr, ...
    'VariableNames', {'Task','n','nBaselineDiff','nRespSigRaw','nRespSigCorrected'});

rho = nan(2,1); md = nan(2,1);
Cs = {CObj, CLoc}; Ts = {T.Obj, T.Loc};
for k = 1:2
    x = Cs{k}.AUCRaw(m); y = Ts{k}.TgtAUC(m); ok = ~isnan(x) & ~isnan(y);
    if sum(ok) > 3
        rho(k) = corr(x(ok), y(ok), 'type', 'Spearman');
        md(k)  = median(abs(x(ok) - y(ok)));
    end
end
B.Validation = table(task, rho, md, 'VariableNames', {'Task','RhoVsSavedAUC','MedianAbsDiff'});
end
