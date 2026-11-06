function Out = ObjLocValue_SingleUnitSummary(ShowFigureFlag, OutputFlag)
% ObjLocValue_SingleUnitSummary  Summarize all single-unit object vs location
% value observations for FEF and SC, test FEF vs SC formally, and produce
% presentation-ready figures (one conclusion per figure).
%
% Runs both areas itself (not a per-batch GUI protocol). From MATLAB:
%   ObjLocValue_SingleUnitSummary(1, 1)
%
% Inputs
%   ShowFigureFlag — 1 to draw the 5 slide figures (default 1)
%   OutputFlag     — 1 to save .mat, comparison .csv, summary .txt and
%                    figures (.png 300 dpi + .fig) (default 0)
% Output
%   Out — struct: Summary.(FEF/SC), Comparison (table), Presentation, Text
%
% Metrics (trial-level ROC on RF trials, see ComputeTrialWindowAUC):
%   coding fractions, pre-target value differences, overlap odds ratio,
%   AUC correlation (Spearman), sign concordance, split-half noise-corrected r,
%   baseline-difference split. FEF vs SC: Fisher / ranksum / z / bootstrap,
%   Holm-corrected (ComputeAreaComparison).
%
% Xuefei Yu, 2026-11-05

if nargin < 1, ShowFigureFlag = 1; end
if nargin < 2, OutputFlag = 0; end

%% ---------------- Paths & parameters ----------------
DataHubPath  = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/DataHub_YXF/Apps/DataHub';
BatchFileDir = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/';
ResultDir    = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/';
SaveDir      = fullfile(ResultDir, 'ObjectLocationPSTH', 'SingleUnitSummary');

Areas      = {'FEF', 'SC'};
BatchFiles = {'FEF_BothRA_StableValue.txt', 'SC_BothRA_StableValue.txt'};

Params.RespWin          = [100 200];
Params.BaseWin          = [-300 0];
Params.TCSmoothMs       = 50;
Params.Alpha            = 0.05;
Params.MinFR            = 1;
Params.MinTrials        = 10;
Params.RequireOpto      = true;
Params.ExcludeOptoFiles = {'Robin060923FRACSAC1'};
Params.Subgroup         = 'All';      % 'All' | 'Fast' | 'Others'
Params.nSplit           = 100;        % split-half repeats
Params.nBoot            = 2000;
Params.Seed             = 1;
FigureStartNum          = 700;

addpath(genpath(DataHubPath));

%% ---------------- Compute per area ----------------
Summary = struct(); Cache = struct();
for a = 1:2
    fprintf('Loading %s ...\n', Areas{a});
    Data = BuildObjLocBatchData(BatchFileDir, ResultDir, BatchFiles{a});
    T    = LoadObjLocNeuronTable(Data.ResultFilePath, Data.FileName, Data.ChannelNumber, ...
                                 Params.ExcludeOptoFiles, true);
    Sel  = ComputeObjLocSelection(T, Params);
    m    = Sel.(Params.Subgroup);
    CObj = ComputeTrialWindowAUC(T.Obj.Trials, T.Obj.TrialTime, m, Params);
    CLoc = ComputeTrialWindowAUC(T.Loc.Trials, T.Loc.TrialTime, m, Params);
    Summary.(Areas{a}) = ComputeSingleUnitSummary(T, CObj, CLoc, m, Params);
    Cache.(Areas{a})   = struct('Obj', CObj, 'Loc', CLoc, 'Mask', m);
end

%% ---------------- Compare areas ----------------
SA = Summary.(Areas{1}); SB = Summary.(Areas{2});
Cmp = ComputeAreaComparison(SA, SB, Areas);
PD  = ComputePresentationData(SA, SB, Cmp, Cache.(Areas{1}), Cache.(Areas{2}), Areas, Params);
Txt = FormatSingleUnitSummaryText(PD, SA, SB, Cmp, Areas);

fprintf('\n%s\n', Txt{:});
fprintf('\n--- Full FEF vs SC comparison ---\n'); disp(Cmp);

%% ---------------- Figures ----------------
figs = gobjects(0);
if ShowFigureFlag
    figs(1) = PlotPresentationFractionBars(PD.Fig1, FigureStartNum + 1, Areas, '% of neurons coding value');
    figs(2) = PlotPresentationFractionBars(PD.Fig2, FigureStartNum + 2, Areas, '% of neurons, pre-target difference');
    figs(3) = PlotPresentationEffectDots(PD.Fig3.E, FigureStartNum + 3, Areas);
    figs(4) = PlotPresentationAUCScatter(PD.Fig4, FigureStartNum + 4, Areas);
    figs(5) = PlotPresentationEffectDots(PD.Fig5.E, FigureStartNum + 5, Areas);
end

%% ---------------- Save ----------------
if OutputFlag
    if ~exist(SaveDir, 'dir'), mkdir(SaveDir); end
    tag = Params.Subgroup;
    save(fullfile(SaveDir, sprintf('SingleUnitSummary_%s.mat', tag)), 'Summary', 'Cmp', 'PD', 'Params', '-v7.3');
    writetable(Cmp, fullfile(SaveDir, sprintf('AreaComparison_%s.csv', tag)));
    fid = fopen(fullfile(SaveDir, sprintf('Observations_%s.txt', tag)), 'w');
    fprintf(fid, '%s\n', Txt{:}); fclose(fid);
    names = {'1_CodingFractions', '2_PreTargetValue', '3_Overlap', '4_AUCScatter', '5_NoiseCorrectedCorrelation'};
    for k = 1:numel(figs)
        f = fullfile(SaveDir, sprintf('Slide%s_%s', names{k}, tag));
        exportgraphics(figs(k), [f '.png'], 'Resolution', 300);
        savefig(figs(k), [f '.fig']);
    end
    fprintf('Saved summary, comparison table and figures -> %s\n', SaveDir);
end

Out = struct('Summary', Summary, 'Comparison', Cmp, 'Presentation', PD, 'Text', {Txt});
end
