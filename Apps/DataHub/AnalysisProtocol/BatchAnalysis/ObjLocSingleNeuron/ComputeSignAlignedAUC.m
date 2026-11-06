function [AObj, ALoc] = ComputeSignAlignedAUC(aucPsthObj, aucPsthLoc, aucObj, aucLoc, ref)
% ComputeSignAlignedAUC  Flip AUC time courses so the reference task's
% preferred value maps to AUC > 0.5.
%
% Inputs
%   aucPsthObj, aucPsthLoc — N-by-T AUC time courses (object / location task)
%   aucObj, aucLoc         — N-vectors, scalar 100-200 ms AUCs (define the sign)
%   ref                    — N-vector: 1 = flip both by object sign (ObjOnly)
%                                      2 = flip both by location sign (LocOnly)
%                                      3 = flip each task by its own sign (Both)
%                                      0 = no flip
% Outputs
%   AObj, ALoc — N-by-T aligned AUC time courses
%
% Note: the non-significant task is never flipped by its own sign (that would
% push its mean above 0.5 by construction).
%
% Xuefei Yu, 2026-11-05

ref = ref(:); aucObj = aucObj(:); aucLoc = aucLoc(:);
flipObj = (ref == 1 & aucObj < 0.5) | (ref == 2 & aucLoc < 0.5) | (ref == 3 & aucObj < 0.5);
flipLoc = (ref == 1 & aucObj < 0.5) | (ref == 2 & aucLoc < 0.5) | (ref == 3 & aucLoc < 0.5);

AObj = aucPsthObj;  AObj(flipObj,:) = 1 - AObj(flipObj,:);
ALoc = aucPsthLoc;  ALoc(flipLoc,:) = 1 - ALoc(flipLoc,:);
end
