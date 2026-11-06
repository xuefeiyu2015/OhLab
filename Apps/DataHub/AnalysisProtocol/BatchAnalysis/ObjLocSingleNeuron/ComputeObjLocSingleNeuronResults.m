function R = ComputeObjLocSingleNeuronResults(T, mask, Params)
% ComputeObjLocSingleNeuronResults  Run all single-neuron computations for one
% subpopulation (e.g. Sel.All / Sel.Fast / Sel.Others).
%
% Inputs
%   T      — struct from LoadObjLocNeuronTable
%   mask   — N-by-1 logical, neurons to include
%   Params — struct: Alpha, NoCodingCriterion, EffectMargin, PreTgtWin,
%            HeatmapSortMode, nBoot, Seed
% Output: struct R
%   .Index        indices into T of the included neurons (n-by-1)
%   .AUCObj/.AUCLoc, .pObj/.pLoc   scalar 100-200 ms AUC / ranksum p (n-by-1)
%   .Groups       from ComputeValueGroups
%   .AlignedObj/.AlignedLoc        sign-aligned AUC time courses (n-by-T)
%   .TCObj/.TCLoc group time courses (ComputeGroupAUCTimecourse)
%   .GroupNamesTC groups used for time courses
%   .Concordance  table (ComputeSignConcordance)
%   .PreTarget    table (ComputePreTargetAUC)
%   .Equivalence  table (ComputeOtherValueEquivalence)
%   .SortIdx / .SortIdxLoc            heatmap row order by object / location AUC (into R.Index)
%   .SortBoundary / .SortBoundaryLoc  row where the sorted AUC crosses 0.5 (signed mode)
%
% Xuefei Yu, 2026-11-05

R.Index  = find(mask(:));
R.AUCObj = T.Obj.TgtAUC(R.Index);   R.pObj = T.Obj.pRanksum(R.Index);
R.AUCLoc = T.Loc.TgtAUC(R.Index);   R.pLoc = T.Loc.pRanksum(R.Index);
psthObj  = T.Obj.AUCPSTH(R.Index,:);
psthLoc  = T.Loc.AUCPSTH(R.Index,:);

R.Groups = ComputeValueGroups(R.AUCObj, R.pObj, R.AUCLoc, R.pLoc, Params);
G = R.Groups;

[R.AlignedObj, R.AlignedLoc] = ComputeSignAlignedAUC(psthObj, psthLoc, R.AUCObj, R.AUCLoc, G.Ref);

R.GroupNamesTC = {'ObjOnly','LocOnly','Both','None'};
masksTC = cellfun(@(nm) G.Mask.(nm), R.GroupNamesTC, 'UniformOutput', false);
R.TCObj = ComputeGroupAUCTimecourse(R.AlignedObj, masksTC, Params.Alpha);
R.TCLoc = ComputeGroupAUCTimecourse(R.AlignedLoc, masksTC, Params.Alpha);

concNames = [R.GroupNamesTC, {'All'}];
concMasks = [masksTC, {true(numel(R.Index),1)}];
R.Concordance = ComputeSignConcordance(R.AUCObj, R.AUCLoc, concMasks, concNames);

R.PreTarget = ComputePreTargetAUC(psthObj, T.Obj.Time, psthLoc, T.Loc.Time, ...
    R.AUCObj, R.AUCLoc, Params.PreTgtWin, concMasks, concNames);

R.Equivalence = ComputeOtherValueEquivalence(R.AUCObj, R.AUCLoc, G, ...
    Params.EffectMargin, Params.nBoot, Params.Seed);

R.SortIdx    = ComputeAUCSortOrder(R.AUCObj, Params.HeatmapSortMode);
R.SortIdxLoc = ComputeAUCSortOrder(R.AUCLoc, Params.HeatmapSortMode);
if strcmpi(Params.HeatmapSortMode, 'signed')
    R.SortBoundary    = sum(R.AUCObj > 0.5) + 0.5;
    R.SortBoundaryLoc = sum(R.AUCLoc > 0.5) + 0.5;
else
    R.SortBoundary    = [];
    R.SortBoundaryLoc = [];
end
R.Params = Params;
end
