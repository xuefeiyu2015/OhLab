function Tv = ComputeAUCVariantTable(T, CObj, CLoc, variant)
% ComputeAUCVariantTable  Copy of T with the AUC fields replaced by a
% trial-level variant, so ComputeObjLocSingleNeuronResults can be reused.
%
% Inputs
%   T          — struct from LoadObjLocNeuronTable (LoadTrials = true)
%   CObj, CLoc — structs from ComputeTrialWindowAUC
%   variant    — 'raw' (no correction) or 'corrected' (baseline subtracted)
% Output
%   Tv — T with Obj/Loc .TgtAUC, .pRanksum, .AUCPSTH, .Time replaced;
%        trial data dropped to keep it light
%
% Xuefei Yu, 2026-11-05

switch lower(variant)
    case 'raw',       f = {'AUCRaw',  'pRaw',  'AUCPSTHRaw'};
    case 'corrected', f = {'AUCCorr', 'pCorr', 'AUCPSTHCorr'};
    otherwise, error('ComputeAUCVariantTable:BadVariant', 'Unknown variant "%s".', variant);
end
Tv = T;
tasks = {'Obj', 'Loc'}; Cs = {CObj, CLoc};
for k = 1:2
    tk = tasks{k}; C = Cs{k};
    Tv.(tk).TgtAUC   = C.(f{1});
    Tv.(tk).pRanksum = C.(f{2});
    Tv.(tk).AUCPSTH  = C.(f{3});
    Tv.(tk).Time     = T.(tk).TrialTime;
    if isfield(Tv.(tk), 'Trials'), Tv.(tk) = rmfield(Tv.(tk), 'Trials'); end
end
end
