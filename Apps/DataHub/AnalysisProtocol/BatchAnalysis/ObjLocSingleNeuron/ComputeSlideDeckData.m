function D = ComputeSlideDeckData(T, CObj, CLoc, mask, Params)
% ComputeSlideDeckData  All per-area numbers for the object vs location value
% slide deck (ObjLocValue_SlideDeck). Pure computation.
%
% Inputs
%   T          — LoadObjLocNeuronTable output (LoadTrials = true)
%   CObj, CLoc — ComputeTrialWindowAUC outputs (trial-level, for baseline
%                split and baseline-corrected AUCs)
%   mask       — N-by-1 logical neurons to include
%   Params     — .Alpha, .NoCodingCriterion, .EffectMargin, .PreTgtWin,
%                .HeatmapSortMode, .nBoot, .Seed
% Output: struct D
%   .R            ComputeObjLocSingleNeuronResults on the protocol's saved AUCs
%                 (scatter, heatmap, 4-group population AUC time courses)
%   .Overlap.All          ComputeOverlapVsChance, saved AUCs, all neurons
%   .Overlap.WithBase     ... neurons WITH a pre-target value difference (oneDR)
%   .Overlap.WithoutBase  ... neurons WITHOUT it
%   .Overlap.Corrected    ... trial-level baseline-subtracted AUCs, all neurons
%   .Frac         coding / anticipation fractions with counts
%
% Xuefei Yu, 2026-11-05

a = Params.Alpha;
D.R = ComputeObjLocSingleNeuronResults(T, mask, Params);
idx = D.R.Index;

ao = D.R.AUCObj; po = D.R.pObj; al = D.R.AUCLoc; pl = D.R.pLoc;
pBase = CLoc.pBase(idx);
withB = pBase < a; withoutB = pBase >= a;   % NaN baseline p -> in neither

D.Overlap.All         = ComputeOverlapVsChance(ao, po, al, pl, a);
D.Overlap.WithBase    = ComputeOverlapVsChance(ao(withB), po(withB), al(withB), pl(withB), a);
D.Overlap.WithoutBase = ComputeOverlapVsChance(ao(withoutB), po(withoutB), al(withoutB), pl(withoutB), a);
D.Overlap.Corrected   = ComputeOverlapVsChance(CObj.AUCCorr(idx), CObj.pCorr(idx), ...
                                               CLoc.AUCCorr(idx), CLoc.pCorr(idx), a);
D.BaseMask = struct('With', withB, 'Without', withoutB);

valid = ~isnan(pBase);
D.Frac.Labels = {'Object value', 'Location value', 'Location anticipation'};
D.Frac.k = [sum(po < a); sum(pl < a); sum(withB)];
D.Frac.n = [sum(~isnan(po)); sum(~isnan(pl)); sum(valid)];
D.Frac.Frac = D.Frac.k ./ D.Frac.n;
end
