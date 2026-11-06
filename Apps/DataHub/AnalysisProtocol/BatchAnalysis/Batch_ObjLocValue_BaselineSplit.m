function Batch_ObjLocValue_BaselineSplit(Data, StartFile, EndFile, ShowFigureFlag, OutputFlag)
% Batch_ObjLocValue_BaselineSplit  Split neurons by whether they show a
% pre-target high vs low difference in the location (oneDR) task, then compare
% the object (stable) vs location value relationship in each subset.
% Run on FEF_BothRA_StableValue or SC_BothRA_StableValue.
%
% Baseline difference: ranksum p < Alpha between high- and low-value RF trials
% on mean FR in BaseWin (oneDR task). Relationship metrics per subset: overlap
% of significant neurons (Fisher), Spearman rho of AUCs, sign concordance;
% computed for raw and baseline-corrected response AUCs.
%
% Inputs : Data (batch struct), StartFile, EndFile, ShowFigureFlag, OutputFlag
% Outputs: figures; .mat -> <SaveRoot>/<Area>/BaselineSplit_<Area>.mat
% Neuron selection: identical to Batch_ObjLocValue_SingleNeuron / Organize.
%
% Xuefei Yu, 2026-11-05

%% ---------------- Parameters ----------------
SaveRoot = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/ObjectLocationPSTH';

Params.RespWin          = [100 200];   % ms
Params.BaseWin          = [-300 0];    % ms, baseline window for the split
Params.TCSmoothMs       = 50;          % ms boxcar, AUC time courses/heatmaps only
Params.Alpha            = 0.05;
Params.DispWin          = [-500 300];
Params.MinFR            = 1;
Params.MinTrials        = 10;
Params.RequireOpto      = true;
Params.ExcludeOptoFiles = {'Robin060923FRACSAC1'};
Params.Subgroups        = {'All', 'Fast', 'Others'};
FigureStartNum          = 500;

if isfield(Data, 'SingleNeuronParams')          % optional overrides
    ov = fieldnames(Data.SingleNeuronParams);
    for k = 1:numel(ov), Params.(ov{k}) = Data.SingleNeuronParams.(ov{k}); end
end

HelperDir = fullfile(fileparts(mfilename('fullpath')), 'ObjLocSingleNeuron');
if ~contains([pathsep, path, pathsep], [pathsep, HelperDir, pathsep]), addpath(HelperDir); end

%% ---------------- Load + select ----------------
T = LoadObjLocNeuronTable(Data.ResultFilePath(StartFile:EndFile), Data.FileName(StartFile:EndFile), ...
    Data.ChannelNumber(StartFile:EndFile,:), Params.ExcludeOptoFiles, true);
Sel = ComputeObjLocSelection(T, Params);
if contains(Data.BatchFileName, 'FEF'),    Area = 'FEF';
elseif contains(Data.BatchFileName, 'SC'), Area = 'SC';
else,                                      Area = 'Unknown';
end

%% ---------------- Compute ----------------
CObj = ComputeTrialWindowAUC(T.Obj.Trials, T.Obj.TrialTime, Sel.All, Params);
CLoc = ComputeTrialWindowAUC(T.Loc.Trials, T.Loc.TrialTime, Sel.All, Params);

Results = struct();
for s = 1:numel(Params.Subgroups)
    sg = Params.Subgroups{s};
    Results.(sg) = ComputeBaselineSplitResults(CObj, CLoc, Sel.(sg), Params);
    S = Results.(sg);
    fprintf('\n=== %s %s — split by oneDR baseline difference ===\n', Area, sg);
    fprintf('-- raw response AUC --\n');       disp(S.OverlapRaw);
    fprintf('-- baseline-corrected AUC --\n'); disp(S.OverlapCorr);
    fprintf('-- baseline preference vs response preference --\n'); disp(S.BaseVsResp);
end

%% ---------------- Plot ----------------
if ShowFigureFlag
    fn = FigureStartNum;
    for s = 1:numel(Params.Subgroups)
        sg = Params.Subgroups{s}; S = Results.(sg);
        if ~any(S.Masks{1}), continue; end
        tag = sprintf('%s %s', Area, sg);

        fn = fn + 1;
        PlotObjLocAUCScatterBySubset(CObj.AUCRaw, CObj.pRaw, CLoc.AUCRaw, CLoc.pRaw, S.Masks(2:3), S.Names(2:3), ...
            Params.Alpha, ['Raw AUC by oneDR baseline difference — ' tag], fn);
        fn = fn + 1;
        PlotObjLocAUCScatterBySubset(CObj.AUCCorr, CObj.pCorr, CLoc.AUCCorr, CLoc.pCorr, S.Masks(2:3), S.Names(2:3), ...
            Params.Alpha, ['Corrected AUC by oneDR baseline difference — ' tag], fn);
        fn = fn + 1;
        PlotOverlapStatsComparison(S.OverlapRaw,  ['Object vs location (raw) by baseline subset — ' tag], fn);
        fn = fn + 1;
        PlotOverlapStatsComparison(S.OverlapCorr, ['Object vs location (corrected) by baseline subset — ' tag], fn);
        fn = fn + 1;
        PlotGroupAUCTimecourse(S.TCObj, T.Obj.TrialTime, S.TCLoc, T.Loc.TrialTime, S.Names, Params.DispWin, ...
            ['AUC time course by baseline subset (own-sign aligned: 100-200 ms inflated by construction; read pre-target) — ' tag], fn);
    end
end

%% ---------------- Save ----------------
if OutputFlag
    SaveDir = fullfile(SaveRoot, Area);
    if ~exist(SaveDir, 'dir'), mkdir(SaveDir); end
    SaveFile = fullfile(SaveDir, sprintf('BaselineSplit_%s.mat', Area));
    save(SaveFile, 'CObj', 'CLoc', 'Sel', 'Results', 'Params', '-v7.3');
    fprintf('Saved -> %s\n', SaveFile);
end
end
