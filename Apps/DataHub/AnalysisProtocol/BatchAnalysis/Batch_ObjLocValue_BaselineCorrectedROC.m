function Batch_ObjLocValue_BaselineCorrectedROC(Data, StartFile, EndFile, ShowFigureFlag, OutputFlag)
% Batch_ObjLocValue_BaselineCorrectedROC  Remove the pre-target baseline
% trial by trial (both StableValue and OneDR), recompute the high- vs low-value
% ROC on the target response, and compare object vs location value coding
% before vs after correction.
% Run on FEF_BothRA_StableValue or SC_BothRA_StableValue.
%
% Method: for each RF trial, corrected response = mean FR in RespWin minus mean
% FR in BaseWin of the same trial. AUC = exact Mann-Whitney P(high > low),
% p = ranksum. Uncorrected ("raw") AUCs use the same trials and method, so
% raw vs corrected comparisons are like-for-like. Single-trial rates are the
% unsmoothed 10 ms bins in 'psth_targeton_raw', so raw AUCs differ slightly
% from the protocol's TgtAUC_RF (smoothed); see Validation table.
%
% Inputs : Data (batch struct), StartFile, EndFile, ShowFigureFlag, OutputFlag
% Outputs: figures; .mat -> <SaveRoot>/<Area>/BaselineCorrectedROC_<Area>.mat
% Neuron selection: identical to Batch_ObjLocValue_SingleNeuron / Organize.
%
% Xuefei Yu, 2026-11-05

%% ---------------- Parameters ----------------
SaveRoot = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/ObjectLocationPSTH';

Params.RespWin           = [100 200];   % ms, target response
Params.BaseWin           = [-300 0];    % ms, pre-target baseline (same trial)
Params.TCSmoothMs       = 50;          % ms boxcar, AUC time courses/heatmaps only
Params.Alpha             = 0.05;
Params.NoCodingCriterion = 'pvalue';
Params.EffectMargin      = 0.05;
Params.PreTgtWin         = [-300 0];
Params.HeatmapSortMode   = 'signed';
Params.HeatmapCLim       = [0.2 0.8];
Params.HeatmapWin        = [-200 400];
Params.DispWin           = [-500 300];
Params.nBoot             = 2000;
Params.Seed              = 1;
Params.MinFR             = 1;
Params.MinTrials         = 10;
Params.RequireOpto       = true;
Params.ExcludeOptoFiles  = {'Robin060923FRACSAC1'};
Params.Subgroups         = {'All', 'Fast', 'Others'};
FigureStartNum           = 400;

if isfield(Data, 'SingleNeuronParams')          % optional overrides
    ov = fieldnames(Data.SingleNeuronParams);
    for k = 1:numel(ov), Params.(ov{k}) = Data.SingleNeuronParams.(ov{k}); end
end

HelperDir = fullfile(fileparts(mfilename('fullpath')), 'ObjLocSingleNeuron');
if ~contains([pathsep, path, pathsep], [pathsep, HelperDir, pathsep]), addpath(HelperDir); end

%% ---------------- Load + select ----------------
T = LoadObjLocNeuronTable(Data.ResultFilePath(StartFile:EndFile), Data.FileName(StartFile:EndFile), ...
    Data.ChannelNumber(StartFile:EndFile,:), Params.ExcludeOptoFiles, true);
Sel  = ComputeObjLocSelection(T, Params);
Area = AreaFromBatchName(Data.BatchFileName);

%% ---------------- Compute ----------------
CObj = ComputeTrialWindowAUC(T.Obj.Trials, T.Obj.TrialTime, Sel.All, Params);
CLoc = ComputeTrialWindowAUC(T.Loc.Trials, T.Loc.TrialTime, Sel.All, Params);

Results = struct();
for s = 1:numel(Params.Subgroups)
    sg = Params.Subgroups{s};
    Results.(sg) = ComputeBaselineCorrectionResults(T, CObj, CLoc, Sel.(sg), Params);
    B = Results.(sg);
    fprintf('\n=== %s %s (n=%d) — baseline-corrected ROC ===\n', Area, sg, sum(Sel.(sg)));
    disp(B.Validation); disp(B.BaselineSig); disp(B.Overlap);
    fprintf('Groups ObjOnly/LocOnly/Both/None  raw: %s   corrected: %s\n', ...
        mat2str(B.Raw.Groups.Counts(1:4)), mat2str(B.Corrected.Groups.Counts(1:4)));
    fprintf('Pre-target AUC (All, own-sign aligned)  raw Loc %.3f / Obj %.3f   corrected Loc %.3f / Obj %.3f\n', ...
        B.Raw.PreTarget.LocMean(end), B.Raw.PreTarget.ObjMean(end), ...
        B.Corrected.PreTarget.LocMean(end), B.Corrected.PreTarget.ObjMean(end));
end

%% ---------------- Plot ----------------
if ShowFigureFlag
    fn = FigureStartNum;
    tObj = T.Obj.TrialTime; tLoc = T.Loc.TrialTime;
    for s = 1:numel(Params.Subgroups)
        sg = Params.Subgroups{s}; B = Results.(sg);
        if isempty(B.Corrected.Index), continue; end
        tag = sprintf('%s %s', Area, sg);

        fn = fn + 1;
        PlotRawVsCorrectedAUC(CObj, CLoc, Sel.(sg), Params.Alpha, ['Raw vs baseline-corrected AUC — ' tag], fn);
        fn = fn + 1;
        PlotOverlapStatsComparison(B.Overlap, ['Object vs location: raw vs corrected — ' tag], fn);

        Rc = B.Corrected;
        fn = fn + 1;
        PlotObjLocAUCScatterBySubset(Rc.AUCObj, Rc.pObj, Rc.AUCLoc, Rc.pLoc, {true(numel(Rc.Index),1)}, {sg}, ...
            Params.Alpha, ['Corrected AUC: object vs location — ' tag], fn);
        fn = fn + 1;
        rows = Rc.Index(Rc.SortIdx);
        PlotAUCHeatmapPair(CObj.AUCPSTHCorr(rows,:), tObj, CLoc.AUCPSTHCorr(rows,:), tLoc, ...
            Params.HeatmapWin, Params.HeatmapCLim, ['Corrected AUC heatmap sorted by object AUC — ' tag], fn, Rc.SortBoundary, 'object AUC');
        fn = fn + 1;
        rows = Rc.Index(Rc.SortIdxLoc);
        PlotAUCHeatmapPair(CObj.AUCPSTHCorr(rows,:), tObj, CLoc.AUCPSTHCorr(rows,:), tLoc, ...
            Params.HeatmapWin, Params.HeatmapCLim, ['Corrected AUC heatmap sorted by location AUC — ' tag], fn, Rc.SortBoundaryLoc, 'location AUC');
        fn = fn + 1;
        PlotGroupAUCTimecourse(Rc.TCObj, tObj, Rc.TCLoc, tLoc, Rc.GroupNamesTC, Params.DispWin, ...
            ['Corrected sign-aligned AUC by group — ' tag], fn);
    end
end

%% ---------------- Save ----------------
if OutputFlag
    SaveDir = fullfile(SaveRoot, Area);
    if ~exist(SaveDir, 'dir'), mkdir(SaveDir); end
    SaveFile = fullfile(SaveDir, sprintf('BaselineCorrectedROC_%s.mat', Area));
    save(SaveFile, 'CObj', 'CLoc', 'Sel', 'Results', 'Params', '-v7.3');
    fprintf('Saved -> %s\n', SaveFile);
end
end

function Area = AreaFromBatchName(name)
% AreaFromBatchName  'FEF' / 'SC' / 'Unknown' from the batch file name.
if contains(name, 'FEF'),    Area = 'FEF';
elseif contains(name, 'SC'), Area = 'SC';
else,                        Area = 'Unknown';
end
end
