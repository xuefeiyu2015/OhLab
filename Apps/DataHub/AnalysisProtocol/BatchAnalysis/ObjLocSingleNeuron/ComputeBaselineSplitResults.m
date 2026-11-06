function S = ComputeBaselineSplitResults(CObj, CLoc, mask, Params)
% ComputeBaselineSplitResults  Split neurons by a pre-target high/low
% difference in the location (oneDR) task, and compare the object-location
% value relationship in each subset.
%
% Inputs
%   CObj, CLoc — structs from ComputeTrialWindowAUC
%   mask       — N-by-1 logical subpopulation
%   Params     — .Alpha (baseline split and coding threshold, default 0.05)
% Output: struct S
%   .Names       {'All','BaselineDiff','NoBaselineDiff'}
%   .Masks       cell of N-by-1 logical
%   .OverlapRaw  ComputeOverlapStats table, uncorrected response AUC
%   .OverlapCorr ComputeOverlapStats table, baseline-corrected response AUC
%   .BaseVsResp  table: does the oneDR baseline preference match the oneDR
%                (raw) response preference? / the object-task preference?
%   .AlignedObj, .AlignedLoc  N-by-T AUC time courses, each task aligned by
%                its own raw response sign (for subset time-course plots)
%                NOTE: aligning by the response-window sign inflates AUC in that
%                window even for null neurons; interpret the pre-target part only.
%   .TCObj, .TCLoc  ComputeGroupAUCTimecourse over S.Masks
%
% Xuefei Yu, 2026-11-05

if ~isfield(Params,'Alpha'), Params.Alpha = 0.05; end
a = Params.Alpha; m = mask(:);
hasBase = ~isnan(CLoc.pBase);
S.Names = {'All', 'BaselineDiff', 'NoBaselineDiff'};
S.Masks = {m & hasBase, m & CLoc.pBase < a, m & CLoc.pBase >= a};

S.OverlapRaw  = ComputeOverlapStats(CObj.AUCRaw,  CObj.pRaw,  CLoc.AUCRaw,  CLoc.pRaw,  S.Masks, S.Names, a);
S.OverlapCorr = ComputeOverlapStats(CObj.AUCCorr, CObj.pCorr, CLoc.AUCCorr, CLoc.pCorr, S.Masks, S.Names, a);

c1 = ComputeSignConcordance(CLoc.AUCBase, CLoc.AUCRaw, S.Masks, S.Names);
c2 = ComputeSignConcordance(CLoc.AUCBase, CObj.AUCRaw, S.Masks, S.Names);
S.BaseVsResp = table(S.Names(:), c1.n, c1.Fraction, c1.pBinom, c2.Fraction, c2.pBinom, 'VariableNames', ...
    {'Subset','n','LocBase_vs_LocResp','p1','LocBase_vs_ObjResp','p2'});

ref = 3 * ones(numel(m), 1);   % each task by its own raw response sign
[S.AlignedObj, S.AlignedLoc] = ComputeSignAlignedAUC(CObj.AUCPSTHRaw, CLoc.AUCPSTHRaw, CObj.AUCRaw, CLoc.AUCRaw, ref);
S.TCObj = ComputeGroupAUCTimecourse(S.AlignedObj, S.Masks, a);
S.TCLoc = ComputeGroupAUCTimecourse(S.AlignedLoc, S.Masks, a);
S.Params = Params;
end
