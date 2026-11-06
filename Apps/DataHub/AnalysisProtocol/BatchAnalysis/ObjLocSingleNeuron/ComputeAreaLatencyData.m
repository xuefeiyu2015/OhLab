function A = ComputeAreaLatencyData(LObj, LLoc, CObj, CLoc, mask, Params)
% ComputeAreaLatencyData  Value and visual timing for one area / subgroup.
%
% Value-coding neurons per task and variant:
%   Corrected  — response-minus-baseline p < Alpha (all neurons; removes the
%                oneDR pre-target anticipation offset)            [main]
%   RawNoAntic — raw response p < Alpha AND no pre-target high/low
%                difference (baseline p >= Alpha)                 [control]
% Population curve = sign-aligned per-bin AUC - 0.5 (sign of the neuron's
% response-window AUC), onset by ComputePopulationOnset.
%
% Inputs
%   LObj, LLoc — ComputeValueLatency outputs (object / location task)
%   CObj, CLoc — ComputeTrialWindowAUC outputs (window AUC + p)
%   mask       — N-by-1 logical neurons in the subgroup
%   Params     — .Alpha + ComputePopulationOnset params
% Output: struct A with A.(Obj|Loc).(Corrected|RawNoAntic|Visual|ValueMinusVisual)
%   value variants: .Index, .n, .Lat (n-by-1, NaN = no onset found),
%                   .nFound, .Pop (ComputePopulationOnset)
%   Visual: .Lat (all mask neurons), .Pop (normalized visual response)
%   ValueMinusVisual: .Lat  (Corrected value latency - visual latency, same neurons)
%   A.(task).t
%
% Xuefei Yu, 2026-11-06

a = Params.Alpha;
tasks = {'Obj', 'Loc'}; Ls = {LObj, LLoc}; Cs = {CObj, CLoc};
for i = 1:2
    L = Ls{i}; C = Cs{i}; tk = tasks{i};
    sel.Corrected  = mask & C.pCorr < a;
    sel.RawNoAntic = mask & C.pRaw < a & C.pBase >= a;
    tcs.Corrected  = (L.AUCTCCorr - 0.5) .* sign(C.AUCCorr - 0.5);
    tcs.RawNoAntic = (L.AUCTCRaw  - 0.5) .* sign(C.AUCRaw  - 0.5);
    lats.Corrected = L.LatCorr; lats.RawNoAntic = L.LatRaw;
    for v = {'Corrected', 'RawNoAntic'}
        vv = v{1}; idx = find(sel.(vv));
        S.Index = idx; S.n = numel(idx); S.Lat = lats.(vv)(idx);
        S.nFound = sum(~isnan(S.Lat));
        S.Pop = ComputePopulationOnset(tcs.(vv)(idx,:), L.t, Params);
        A.(tk).(vv) = S;
    end
    iv = find(mask);
    A.(tk).Visual = struct('Index', iv, 'n', numel(iv), 'Lat', L.VisLat(iv), ...
        'nFound', sum(~isnan(L.VisLat(iv))), 'Pop', ComputePopulationOnset(L.VisTC(iv,:), L.t, Params));
    ic = A.(tk).Corrected.Index;
    A.(tk).ValueMinusVisual = struct('Index', ic, 'n', numel(ic), 'Lat', L.LatCorr(ic) - L.VisLat(ic), ...
        'nFound', sum(~isnan(L.LatCorr(ic) - L.VisLat(ic))));
    A.(tk).t = L.t;
end
end
