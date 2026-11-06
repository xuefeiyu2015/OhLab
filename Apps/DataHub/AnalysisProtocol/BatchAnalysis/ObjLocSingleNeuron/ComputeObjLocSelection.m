function Sel = ComputeObjLocSelection(T, Params)
% ComputeObjLocSelection  Neuron inclusion masks + opto ModuType.
% Matches the selection used in Batch_ObjLocValue_Organize (Sel_All/Fast/Others).
%
% Inputs
%   T      — struct from LoadObjLocNeuronTable
%   Params — struct with fields
%            .MinFR        (default 1 sp/s, both tasks)
%            .MinTrials    (default 10; per RF condition: high/low x obj/loc)
%            .RequireOpto  (default true; as in Organize's Sel_All)
%            .FastLatency  (default 5 ms)
%            .OptoAlpha    (default 0.05)
% Output: struct Sel (all N-by-1 logical unless noted)
%   .ModuType  N-by-1: 1 Fast Ex, 2 Slow Ex, 3 Inhibition, 4 No effect, NaN no opto
%   .All, .Fast, .Others
%
% Xuefei Yu, 2026-11-05

if ~isfield(Params,'MinFR'),       Params.MinFR = 1;        end
if ~isfield(Params,'MinTrials'),   Params.MinTrials = 10;   end
if ~isfield(Params,'RequireOpto'), Params.RequireOpto = true; end
if ~isfield(Params,'FastLatency'), Params.FastLatency = 5;  end
if ~isfield(Params,'OptoAlpha'),   Params.OptoAlpha = 0.05; end

lat = T.Opto.Latency; sm = T.Opto.SM; smp = T.Opto.SM_p;
sigOpto = smp < Params.OptoAlpha & ~isnan(lat);

ModuType = nan(size(lat));
ModuType(T.HasOpto & ~sigOpto)                              = 4;
ModuType(sigOpto & sm <= 0)                                 = 3;
ModuType(sigOpto & sm > 0 & lat >= Params.FastLatency)      = 2;
ModuType(sigOpto & sm > 0 & lat <  Params.FastLatency)      = 1;

responseValid = T.Loc.FRMean > Params.MinFR & T.Obj.FRMean > Params.MinFR;
visualResp    = T.Obj.pVis < 0.05 & T.Obj.VisDiff > 0;   % stable task, as in Organize

trialOK = T.Obj.nHighRF >= Params.MinTrials & T.Obj.nLowRF >= Params.MinTrials & ...
          T.Loc.nHighRF >= Params.MinTrials & T.Loc.nLowRF >= Params.MinTrials;

base = T.HasStable & T.HasOneDR & responseValid & visualResp & trialOK;
if Params.RequireOpto
    base = base & T.HasOpto;
end

Sel.ModuType = ModuType;
Sel.All      = base;
Sel.Fast     = base & ModuType == 1;
Sel.Others   = base & ModuType == 4;
end
