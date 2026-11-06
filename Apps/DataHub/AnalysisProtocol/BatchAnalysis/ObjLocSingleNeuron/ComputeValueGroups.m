function G = ComputeValueGroups(aucObj, pObj, aucLoc, pLoc, Params)
% ComputeValueGroups  Classify neurons by object / location value coding.
%
% Inputs
%   aucObj, pObj — N-vectors: object-task AUC (100-200 ms) and ranksum p
%   aucLoc, pLoc — N-vectors: location-task AUC and ranksum p
%   Params       — struct with fields
%       .Alpha             coding threshold, p < Alpha (default 0.05)
%       .NoCodingCriterion 'pvalue' (p >= Alpha, default)
%                          'effect' (|AUC-0.5| < EffectMargin)
%                          'both'   (p >= Alpha AND |AUC-0.5| < EffectMargin)
%       .EffectMargin      default 0.05
% Output: struct G
%   .Label  N-by-1 integer code, see .Names
%   .Names  {'ObjOnly','LocOnly','Both','None','Ambiguous','Unclassified'}
%   .Mask   struct of N-by-1 logical per name
%   .Ref    N-by-1 reference task for sign alignment: 1 = obj, 2 = loc,
%           3 = own sign (Both), 0 = no flip
%   .Counts 1-by-6 counts (sum == N)
%
% Xuefei Yu, 2026-11-05

if ~isfield(Params,'Alpha'),             Params.Alpha = 0.05;               end
if ~isfield(Params,'NoCodingCriterion'), Params.NoCodingCriterion = 'pvalue'; end
if ~isfield(Params,'EffectMargin'),      Params.EffectMargin = 0.05;        end

aucObj = aucObj(:); pObj = pObj(:); aucLoc = aucLoc(:); pLoc = pLoc(:);

valid = ~isnan(aucObj) & ~isnan(pObj) & ~isnan(aucLoc) & ~isnan(pLoc);

codeObj = pObj < Params.Alpha;
codeLoc = pLoc < Params.Alpha;
noObj   = NoCoding(aucObj, pObj, Params) & ~codeObj;
noLoc   = NoCoding(aucLoc, pLoc, Params) & ~codeLoc;

G.Names = {'ObjOnly','LocOnly','Both','None','Ambiguous','Unclassified'};
M.ObjOnly      = valid & codeObj & noLoc;
M.LocOnly      = valid & codeLoc & noObj;
M.Both         = valid & codeObj & codeLoc;
M.None         = valid & noObj & noLoc;
M.Unclassified = ~valid;
M.Ambiguous    = valid & ~(M.ObjOnly | M.LocOnly | M.Both | M.None);

G.Label = zeros(numel(aucObj),1);
G.Counts = zeros(1, numel(G.Names));
for k = 1:numel(G.Names)
    G.Label(M.(G.Names{k})) = k;
    G.Counts(k) = sum(M.(G.Names{k}));
end
G.Mask = M;

G.Ref = zeros(numel(aucObj),1);
G.Ref(M.ObjOnly) = 1;
G.Ref(M.LocOnly) = 2;
G.Ref(M.Both)    = 3;
G.Params = Params;
end

function tf = NoCoding(auc, p, Params)
switch lower(Params.NoCodingCriterion)
    case 'pvalue', tf = p >= Params.Alpha;
    case 'effect', tf = abs(auc - 0.5) < Params.EffectMargin;
    case 'both',   tf = p >= Params.Alpha & abs(auc - 0.5) < Params.EffectMargin;
    otherwise
        error('ComputeValueGroups:BadCriterion', 'Unknown NoCodingCriterion "%s".', Params.NoCodingCriterion);
end
end
