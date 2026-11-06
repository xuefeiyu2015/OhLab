function Out = ObjLocValue_SlideDeck(ShowFigureFlag, OutputFlag)
% ObjLocValue_SlideDeck  Presentation figures for single-unit object vs
% location value coding in FEF and SC, in talk order:
%   1  Scatter + marginal histograms + Venn            (per area)
%   2  Population AUC heatmap sorted by object value   (per area)
%   3  Population-average AUC for the 4 neuron types   (per area)
%   4  Venn as bars vs chance (independence)           (FEF | SC)
%   5  Split by location anticipation (pre-target)     (FEF/SC x with/without)
%   6  Baseline-corrected responses vs chance          (FEF | SC)
%   7  FEF vs SC: proportions, and overlap OR in every analysis
% Headlines are generated from the results.
%
% Slides 1-5 use the protocols' saved RF AUCs (TgtAUC_RF / AUCPSTH, same as
% Organize); the anticipation split and slide 6 use trial-level ROC
% (ComputeTrialWindowAUC). Neuron selection = Organize (FR, visual, trials, opto).
%
% Inputs : ShowFigureFlag (default 1), OutputFlag (default 0: save PNG/FIG,
%          comparison CSV and .mat to <ResultDir>/ObjectLocationPSTH/SlideDeck)
% Output : Out — struct with Deck.(FEF/SC), Comparison
%
% Xuefei Yu, 2026-11-05

if nargin < 1, ShowFigureFlag = 1; end
if nargin < 2, OutputFlag = 0; end

%% ---------------- Paths & parameters ----------------
DataHubPath  = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/DataHub_YXF/Apps/DataHub';
BatchFileDir = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/';
ResultDir    = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/';
SaveDir      = fullfile(ResultDir, 'ObjectLocationPSTH', 'SlideDeck');

Areas      = {'FEF', 'SC'};
BatchFiles = {'FEF_BothRA_StableValue.txt', 'SC_BothRA_StableValue.txt'};

Params.RespWin           = [100 200];
Params.BaseWin           = [-300 0];
Params.TCSmoothMs        = 50;
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
Params.Subgroup          = 'All';
FigBase                  = 800;

addpath(genpath(DataHubPath));

%% ---------------- Compute ----------------
Deck = struct(); Times = struct();
for a = 1:2
    fprintf('Loading %s ...\n', Areas{a});
    Data = BuildObjLocBatchData(BatchFileDir, ResultDir, BatchFiles{a});
    T    = LoadObjLocNeuronTable(Data.ResultFilePath, Data.FileName, Data.ChannelNumber, Params.ExcludeOptoFiles, true);
    Sel  = ComputeObjLocSelection(T, Params);
    m    = Sel.(Params.Subgroup);
    CObj = ComputeTrialWindowAUC(T.Obj.Trials, T.Obj.TrialTime, m, Params);
    CLoc = ComputeTrialWindowAUC(T.Loc.Trials, T.Loc.TrialTime, m, Params);
    Deck.(Areas{a}) = ComputeSlideDeckData(T, CObj, CLoc, m, Params);
    Deck.(Areas{a}).AUCPSTHObj = T.Obj.AUCPSTH;   % saved AUC time courses (rows = T neurons)
    Deck.(Areas{a}).AUCPSTHLoc = T.Loc.AUCPSTH;
    Times.(Areas{a}) = struct('Obj', T.Obj.Time, 'Loc', T.Loc.Time);
end
DA = Deck.(Areas{1}); DB = Deck.(Areas{2});
K  = ComputeSlideDeckComparison(DA, DB, Areas, Params);

PrintDeckSummary(Deck, K, Areas);

%% ---------------- Figures ----------------
figs = gobjects(0); fnames = {};
if ShowFigureFlag
    for a = 1:2
        A = Areas{a}; D = Deck.(A); R = D.R; O = D.Overlap.All; tt = Times.(A);
        % 1  scatter + Venn
        h = sprintf('%s: %d neurons code object value, %d location value, %d both', A, O.nObj, O.nLoc, O.Obs(2));
        figs(end+1) = PlotSlideScatterVenn(R.AUCObj, R.pObj, R.AUCLoc, R.pLoc, O, Params.Alpha, FigBase + 10*a + 1, h); %#ok<AGROW>
        fnames{end+1} = sprintf('Slide1_ScatterVenn_%s', A); %#ok<AGROW>
        % 2  heatmap sorted by object value
        rows = R.Index(R.SortIdx);
        h = sprintf('%s: location value is only loosely ordered by object value (r = %.2f)', A, O.Rho);
        figs(end+1) = PlotAUCHeatmapPair(D.AUCPSTHObj(rows,:), tt.Obj, D.AUCPSTHLoc(rows,:), tt.Loc, ...
            Params.HeatmapWin, Params.HeatmapCLim, h, FigBase + 10*a + 2, R.SortBoundary, 'object AUC'); %#ok<AGROW>
        fnames{end+1} = sprintf('Slide2_HeatmapByObject_%s', A); %#ok<AGROW>
        % 3  population AUC for the 4 neuron types
        h = sprintf('%s: population AUC by neuron type (sign-aligned to the significant value)', A);
        figs(end+1) = PlotGroupAUCTimecourse(R.TCObj, tt.Obj, R.TCLoc, tt.Loc, ...
            {'Object only', 'Location only', 'Both', 'Neither'}, Params.DispWin, h, FigBase + 10*a + 3); %#ok<AGROW>
        fnames{end+1} = sprintf('Slide3_GroupAUC_%s', A); %#ok<AGROW>
    end
    % 4  Venn as bars vs chance
    figs(end+1) = PlotSlideOverlapVsChance({DA.Overlap.All, DB.Overlap.All}, Areas, [1 2], FigBase + 41, ...
        K.Headline.Overlap, K.Subtitle.Overlap);
    fnames{end+1} = 'Slide4_OverlapVsChance';
    % 5  split by anticipation
    figs(end+1) = PlotSlideOverlapVsChance({DA.Overlap.WithBase, DA.Overlap.WithoutBase, DB.Overlap.WithBase, DB.Overlap.WithoutBase}, ...
        {[Areas{1} ' with anticipation'], [Areas{1} ' without anticipation'], [Areas{2} ' with anticipation'], [Areas{2} ' without anticipation']}, ...
        [2 2], FigBase + 51, K.Headline.Split, K.Subtitle.Split);
    fnames{end+1} = 'Slide5_SplitByAnticipation';
    % 6  baseline-corrected
    figs(end+1) = PlotSlideOverlapVsChance({DA.Overlap.Corrected, DB.Overlap.Corrected}, Areas, [1 2], FigBase + 61, ...
        K.Headline.Corrected, K.Subtitle.Corrected);
    fnames{end+1} = 'Slide6_BaselineCorrected';
    % 7  FEF vs SC
    figs(end+1) = PlotPresentationFractionBars(K.FracBars, FigBase + 71, Areas, '% of neurons');
    fnames{end+1} = 'Slide7a_AreaProportions';
    figs(end+1) = PlotPresentationEffectDots(K.ORDots, FigBase + 72, Areas);
    fnames{end+1} = 'Slide7b_AreaOverlapOR';
end

%% ---------------- Save ----------------
if OutputFlag
    if ~exist(SaveDir, 'dir'), mkdir(SaveDir); end
    for k = 1:numel(figs)
        f = fullfile(SaveDir, fnames{k});
        exportgraphics(figs(k), [f '.png'], 'Resolution', 300);
        savefig(figs(k), [f '.fig']);
    end
    writetable(K.Table, fullfile(SaveDir, 'AreaComparison.csv'));
    save(fullfile(SaveDir, 'SlideDeckData.mat'), 'Deck', 'K', 'Params', '-v7.3');
    fprintf('Saved %d figures + tables -> %s\n', numel(figs), SaveDir);
end
Out = struct('Deck', Deck, 'Comparison', K);
end

% ---------------------------------------------------------------------------
function PrintDeckSummary(Deck, K, Areas)
% PrintDeckSummary  Console overview of the numbers behind each slide.
cond = {'All', 'WithBase', 'WithoutBase', 'Corrected'};
fprintf('\n%-6s %-12s %5s %5s %5s %6s %7s %6s %8s  %s\n', 'Area', 'Analysis', 'N', 'nObj', 'nLoc', 'Both', 'Chance', 'OR', 'pFisher', 'verdict');
for a = 1:2
    for c = 1:numel(cond)
        O = Deck.(Areas{a}).Overlap.(cond{c});
        fprintf('%-6s %-12s %5d %5d %5d %6d %7.1f %6.2f %8.3g  %s\n', Areas{a}, cond{c}, O.N, O.nObj, O.nLoc, ...
            O.Obs(2), O.Exp(2), O.OR, O.pFisher, O.Verdict);
    end
end
fprintf('\nHeadlines:\n  4: %s\n  5: %s\n  6: %s\n  7a: %s\n  7b: %s\n', K.Headline.Overlap, K.Headline.Split, ...
    K.Headline.Corrected, K.Headline.AreaFrac, K.Headline.AreaOR);
fprintf('\nFEF vs SC tests (Holm-corrected):\n'); disp(K.Table);
end
