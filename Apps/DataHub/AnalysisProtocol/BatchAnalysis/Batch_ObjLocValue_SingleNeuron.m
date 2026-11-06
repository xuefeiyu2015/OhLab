function Batch_ObjLocValue_SingleNeuron(Data, StartFile, EndFile, ShowFigureFlag, OutputFlag)
% Batch_ObjLocValue_SingleNeuron  Single-neuron comparison of object value
% (StableValue) vs location value (OneDR) coding inside the RF.
% Run on FEF_BothRA_StableValue or SC_BothRA_StableValue; each result file
% holds both OutputData.StableValue and OutputData.OneDR.
%
% Inputs
%   Data           — batch struct (BatchFileName, ResultFilePath, FileName, ChannelNumber)
%   StartFile/EndFile — file range within the batch
%   ShowFigureFlag — 1 to draw figures
%   OutputFlag     — 1 to save results
% Outputs: none (figures; .mat saved to <SaveRoot>/<Area>/SingleNeuron_<Area>_<criterion>.mat)
%
% Helpers: BatchAnalysis/ObjLocSingleNeuron/ (Compute* = pure, Plot* = render only)
% Plan: DataHub/plans/v2.md
% Xuefei Yu, 2026-11-05

%% ---------------- Parameters ----------------
SaveRoot = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/ObjectLocationPSTH';

Params.Alpha             = 0.05;      % coding: ranksum p < Alpha (100-200 ms window)
Params.NoCodingCriterion = 'pvalue';  % 'pvalue' | 'effect' | 'both'
Params.EffectMargin      = 0.05;      % |AUC-0.5| < margin (for 'effect'/'both' and TOST)
Params.PreTgtWin         = [-300 0];  % ms, pre-target window
Params.HeatmapSortMode   = 'signed';  % 'signed' | 'abs'
Params.HeatmapCLim       = [0.2 0.8];
Params.DispWin           = [-500 300];
Params.HeatmapWin        = [-200 400];  % ms, heatmap x-range
Params.nBoot             = 2000;
Params.Seed              = 1;
Params.MinFR             = 1;
Params.MinTrials         = 10;        % per RF condition (high/low x obj/loc), as in Organize
Params.RequireOpto       = true;      % match Organize's Sel_All
Params.ExcludeOptoFiles  = {'Robin060923FRACSAC1'};
Subgroups                = {'All', 'Fast', 'Others'};
FigureStartNum           = 300;

% Optional overrides: Data.SingleNeuronParams.<field> replaces Params.<field>
% e.g. Data.SingleNeuronParams.NoCodingCriterion = 'effect';
if isfield(Data, 'SingleNeuronParams')
    ov = fieldnames(Data.SingleNeuronParams);
    for k = 1:numel(ov)
        Params.(ov{k}) = Data.SingleNeuronParams.(ov{k});
    end
end

%% ---------------- Helper path ----------------
% Helpers live next to this file; add them so GUI runs work without genpath.
HelperDir = fullfile(fileparts(mfilename('fullpath')), 'ObjLocSingleNeuron');
if ~contains([pathsep, path, pathsep], [pathsep, HelperDir, pathsep])
    addpath(HelperDir);
end

%% ---------------- Load ----------------
FilesName = Data.ResultFilePath(StartFile:EndFile);
T = LoadObjLocNeuronTable(FilesName, Data.FileName(StartFile:EndFile), ...
    Data.ChannelNumber(StartFile:EndFile,:), Params.ExcludeOptoFiles);

if contains(Data.BatchFileName, 'FEF'),    Area = 'FEF';
elseif contains(Data.BatchFileName, 'SC'), Area = 'SC';
else,                                      Area = 'Unknown';
end

Sel = ComputeObjLocSelection(T, Params);

%% ---------------- Compute per subgroup ----------------
Results = struct();
for s = 1:numel(Subgroups)
    sg = Subgroups{s};
    Results.(sg) = ComputeObjLocSingleNeuronResults(T, Sel.(sg), Params);
    G = Results.(sg).Groups;
    fprintf('\n=== %s %s (n=%d, no-coding = %s) ===\n', Area, sg, numel(Results.(sg).Index), Params.NoCodingCriterion);
    for k = 1:numel(G.Names)
        fprintf('  %-12s : %d\n', G.Names{k}, G.Counts(k));
    end
    disp(Results.(sg).Concordance);
    disp(Results.(sg).PreTarget);
    disp(Results.(sg).Equivalence);
end

%% ---------------- Plot ----------------
if ShowFigureFlag
    fn = FigureStartNum;
    for s = 1:numel(Subgroups)
        sg = Subgroups{s};
        R  = Results.(sg);
        if isempty(R.Index), continue; end
        tag = sprintf('%s %s [%s]', Area, sg, Params.NoCodingCriterion);

        fn = fn + 1;
        rows = R.Index(R.SortIdx);   % raw (unflipped) AUC, sorted by object AUC
        PlotAUCHeatmapPair(T.Obj.AUCPSTH(rows,:), T.Obj.Time, ...
                           T.Loc.AUCPSTH(rows,:), T.Loc.Time, ...
                           Params.HeatmapWin, Params.HeatmapCLim, ...
                           ['AUC heatmap sorted by object AUC — ' tag], fn, R.SortBoundary, 'object AUC');
        fn = fn + 1;
        rowsLoc = R.Index(R.SortIdxLoc);   % same, sorted by location AUC
        PlotAUCHeatmapPair(T.Obj.AUCPSTH(rowsLoc,:), T.Obj.Time, ...
                           T.Loc.AUCPSTH(rowsLoc,:), T.Loc.Time, ...
                           Params.HeatmapWin, Params.HeatmapCLim, ...
                           ['AUC heatmap sorted by location AUC — ' tag], fn, R.SortBoundaryLoc, 'location AUC');
        fn = fn + 1;
        PlotGroupAUCTimecourse(R.TCObj, T.Obj.Time, R.TCLoc, T.Loc.Time, R.GroupNamesTC, ...
                               Params.DispWin, ['Sign-aligned AUC by group — ' tag], fn);
        fn = fn + 1;
        PlotSignConcordance(R.Concordance, ['Preferred-value sign concordance — ' tag], fn);
        fn = fn + 1;
        PlotOtherValueAUCDistribution(R.Equivalence, ['Other-value AUC — ' tag], fn);
        fn = fn + 1;
        PlotPreTargetAUC(R.PreTarget, ['Pre-target AUC — ' tag], fn);
    end
end

%% ---------------- Save ----------------
if OutputFlag
    SaveDir = fullfile(SaveRoot, Area);
    if ~exist(SaveDir, 'dir'), mkdir(SaveDir); end
    SaveFile = fullfile(SaveDir, sprintf('SingleNeuron_%s_%s.mat', Area, Params.NoCodingCriterion));
    save(SaveFile, 'T', 'Sel', 'Results', 'Params', '-v7.3');
    fprintf('Saved -> %s\n', SaveFile);
end
end
