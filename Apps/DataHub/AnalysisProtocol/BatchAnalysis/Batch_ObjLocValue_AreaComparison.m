function Batch_ObjLocValue_AreaComparison(Data, StartFile, EndFile, ShowFigureFlag, OutputFlag)
% Batch_ObjLocValue_AreaComparison  FEF vs SC comparison of single-unit object
% (StableValue) and location (OneDR) value coding, runnable from the DataHub GUI.
% Main question: TIMING — in which area does the value signal rise first?
%
% Load FEFSC_BothRA_StableValue.txt in the GUI (recommended). Each row is
% assigned to FEF or SC by matching session + channel range against
% FEF_BothRA_StableValue.txt / SC_BothRA_StableValue.txt. If the loaded batch
% contains only one area, the other area's batch file is loaded automatically.
%
% Timing (per subgroup All / Fast / Others, per task):
%   value onset — per-neuron first run of LatConsecBins significant high-vs-low
%   bins (causal smoothing), and population half-max onset of the sign-aligned
%   AUC with a bootstrap over neurons. Main = baseline-corrected responses (removes
%   oneDR anticipation); control = raw responses in neurons without anticipation.
%   Visual-response onset of the same neurons as reference, and value - visual
%   lag per neuron (is an SC lead just faster visual input?).
%   FEF vs SC: bootstrap onset difference, ranksum on latencies, Holm.
% Overlap (Params.RunOverlap): object/location coding overlap vs chance and
%   coding fractions (ComputeSlideDeckComparison).
% Same neuron selection as Organize (FR, visual, trial count, opto).
%
% Inputs : Data (batch struct), StartFile, EndFile, ShowFigureFlag, OutputFlag
% Outputs: timing figures (value onset, control, visual onset, value - visual);
%          overlap figures if RunOverlap; if OutputFlag: <ResultDir>/
%          ObjectLocationPSTH/AreaComparison/ PNG+FIG, Timing_<Subgroup>.csv,
%          AreaComparison_<Subgroup>.csv, AreaComparison.mat
%
% Xuefei Yu, 2026-11-06

%% ---------------- Paths & parameters ----------------
BatchFileDir = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/';
ResultDir    = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/';
SaveDir      = fullfile(ResultDir, 'ObjectLocationPSTH', 'AreaComparison');

Areas    = {'FEF', 'SC'};
RefFiles = {'FEF_BothRA_StableValue.txt', 'SC_BothRA_StableValue.txt'};

Params.RespWin           = [100 200];
Params.BaseWin           = [-300 0];
Params.TCSmoothMs        = 50;
Params.Alpha             = 0.05;
Params.NoCodingCriterion = 'pvalue';
Params.EffectMargin      = 0.05;
Params.PreTgtWin         = [-300 0];
Params.HeatmapSortMode   = 'signed';
Params.nBoot             = 2000;
Params.Seed              = 1;
Params.MinFR             = 1;
Params.MinTrials         = 10;
Params.RequireOpto       = true;
Params.ExcludeOptoFiles  = {'Robin060923FRACSAC1'};
Params.Subgroups         = {'All', 'Fast', 'Others'};   % computed + printed
Params.PlotSubgroups     = {'All'};                     % plotted / saved as figures
Params.RunOverlap        = true;                        % also the overlap-vs-chance figures
% timing
Params.LatSmoothMs       = 30;          % causal boxcar (ms)
Params.LatConsecBins     = 5;           % consecutive significant 10 ms bins
Params.LatAlpha          = 0.05;        % per-bin threshold
Params.LatSearchWin      = [0 400];     % ms, onset / peak search
Params.LatBaseWin        = [-100 0];    % ms, baseline for population half-max
Params.LatDispWin        = [-100 400];  % ms, curve display
FigureStartNum           = 900;

if isfield(Data, 'SingleNeuronParams')          % optional overrides
    ov = fieldnames(Data.SingleNeuronParams);
    for k = 1:numel(ov), Params.(ov{k}) = Data.SingleNeuronParams.(ov{k}); end
end

HelperDir = fullfile(fileparts(mfilename('fullpath')), 'ObjLocSingleNeuron');
if ~contains([pathsep, path, pathsep], [pathsep, HelperDir, pathsep]), addpath(HelperDir); end

%% ---------------- Assign rows to areas ----------------
Rows = SubsetBatch(Data, StartFile:EndFile);
Ref  = struct('Area', {}, 'FileName', {}, 'ChannelNumber', {}, 'Data', {});
for a = 1:numel(Areas)
    D = BuildObjLocBatchData(BatchFileDir, ResultDir, RefFiles{a});
    Ref(a) = struct('Area', Areas{a}, 'FileName', {D.FileName}, 'ChannelNumber', D.ChannelNumber, 'Data', D);
end
RowArea = ComputeBatchRowArea(Rows.FileName, Rows.ChannelNumber, Ref);
if any(cellfun(@isempty, RowArea))
    fprintf('Skipping %d rows not found in %s / %s:\n', sum(cellfun(@isempty, RowArea)), RefFiles{:});
    fprintf('  %s\n', Rows.FileName{cellfun(@isempty, RowArea)});
end

AreaData = struct();
for a = 1:numel(Areas)
    idx = find(strcmp(RowArea, Areas{a}));
    if isempty(idx)
        fprintf('Loaded batch has no %s rows -> loading %s\n', Areas{a}, RefFiles{a});
        AreaData.(Areas{a}) = Ref(a).Data;
    else
        AreaData.(Areas{a}) = SubsetBatch(Rows, idx);
    end
    fprintf('%s: %d sessions\n', Areas{a}, numel(AreaData.(Areas{a}).FileName));
end

%% ---------------- Compute ----------------
Deck = struct(); Lat = struct();
for a = 1:numel(Areas)
    A = Areas{a}; D = AreaData.(A);
    fprintf('Loading %s neurons ...\n', A);
    T    = LoadObjLocNeuronTable(D.ResultFilePath, D.FileName, D.ChannelNumber, Params.ExcludeOptoFiles, true);
    Sel  = ComputeObjLocSelection(T, Params);
    CObj = ComputeTrialWindowAUC(T.Obj.Trials, T.Obj.TrialTime, Sel.All, Params);
    CLoc = ComputeTrialWindowAUC(T.Loc.Trials, T.Loc.TrialTime, Sel.All, Params);
    fprintf('Computing %s latencies ...\n', A);
    LObj = ComputeValueLatency(T.Obj.Trials, T.Obj.TrialTime, Sel.All, Params);
    LLoc = ComputeValueLatency(T.Loc.Trials, T.Loc.TrialTime, Sel.All, Params);
    for s = 1:numel(Params.Subgroups)
        sg = Params.Subgroups{s};
        Lat.(sg).(A) = ComputeAreaLatencyData(LObj, LLoc, CObj, CLoc, Sel.(sg), Params);
        if Params.RunOverlap
            Deck.(sg).(A) = ComputeSlideDeckData(T, CObj, CLoc, Sel.(sg), Params);
        end
    end
end

LatCmp = struct();
for s = 1:numel(Params.Subgroups)
    sg = Params.Subgroups{s};
    LatCmp.(sg) = ComputeLatencyAreaComparison(Lat.(sg).(Areas{1}), Lat.(sg).(Areas{2}), Areas, Params);
    fprintf('\n===== %s neurons: TIMING (Diff = %s - %s, ms) =====\n', sg, Areas{2}, Areas{1});
    disp(LatCmp.(sg).Table);
    fprintf('  %s\n', LatCmp.(sg).Headline.Value, LatCmp.(sg).Headline.Visual, LatCmp.(sg).Headline.Relative);
end

Cmp = struct();
for s = 1:numel(Params.Subgroups) * Params.RunOverlap
    sg = Params.Subgroups{s};
    Cmp.(sg) = ComputeSlideDeckComparison(Deck.(sg).(Areas{1}), Deck.(sg).(Areas{2}), Areas, Params);
    PrintComparison(Deck.(sg), Cmp.(sg), Areas, sg);
end

%% ---------------- Figures ----------------
figs = gobjects(0); fnames = {};
if ShowFigureFlag
    for s = 1:numel(Params.PlotSubgroups)
        sg = Params.PlotSubgroups{s};
        if ~isfield(LatCmp, sg), continue; end
        fb = FigureStartNum + 20*(s-1);
        pre = [sg ' neurons — '];
        % ---- timing (main)
        LA = Lat.(sg).(Areas{1}); LB = Lat.(sg).(Areas{2}); KL = LatCmp.(sg);
        sigs = {'Corrected', 'Visual', 'ValueMinusVisual', 'RawNoAntic'};
        hl   = {KL.Headline.Value, KL.Headline.Visual, KL.Headline.Relative, KL.Headline.Control};
        sb   = {KL.Subtitle.Value, KL.Subtitle.Visual, KL.Subtitle.Relative, KL.Subtitle.Control};
        fn   = {'Timing_ValueOnset', 'Timing_VisualOnset', 'Timing_ValueMinusVisual', 'Timing_ValueOnsetRawNoAntic'};
        for q = 1:numel(sigs)
            figs(end+1) = PlotLatencySlide(LatencyCols(LA, LB, KL.Table, sigs{q}, Params), ...
                Areas, fb + 10 + q, [pre hl{q}], sb{q}); %#ok<AGROW>
            fnames{end+1} = [fn{q} '_' sg]; %#ok<AGROW>
        end
        if ~isfield(Cmp, sg), continue; end
        % ---- overlap
        K  = Cmp.(sg); DA = Deck.(sg).(Areas{1}); DB = Deck.(sg).(Areas{2});
        for a = 1:numel(Areas)
            A = Areas{a}; R = Deck.(sg).(A).R; O = Deck.(sg).(A).Overlap.All;
            h = sprintf('%s%s: %d object, %d location, %d both', pre, A, O.nObj, O.nLoc, O.Obs(2));
            figs(end+1) = PlotSlideScatterVenn(R.AUCObj, R.pObj, R.AUCLoc, R.pLoc, O, Params.Alpha, fb + a, h); %#ok<AGROW>
            fnames{end+1} = sprintf('ScatterVenn_%s_%s', A, sg); %#ok<AGROW>
        end
        figs(end+1) = PlotSlideOverlapVsChance({DA.Overlap.All, DB.Overlap.All}, Areas, [1 2], fb + 4, ...
            [pre K.Headline.Overlap], K.Subtitle.Overlap); %#ok<AGROW>
        fnames{end+1} = ['OverlapVsChance_' sg]; %#ok<AGROW>
        figs(end+1) = PlotSlideOverlapVsChance({DA.Overlap.WithBase, DA.Overlap.WithoutBase, DB.Overlap.WithBase, DB.Overlap.WithoutBase}, ...
            {[Areas{1} ' with anticipation'], [Areas{1} ' without anticipation'], [Areas{2} ' with anticipation'], [Areas{2} ' without anticipation']}, ...
            [2 2], fb + 5, [pre K.Headline.Split], K.Subtitle.Split); %#ok<AGROW>
        fnames{end+1} = ['SplitByAnticipation_' sg]; %#ok<AGROW>
        figs(end+1) = PlotSlideOverlapVsChance({DA.Overlap.Corrected, DB.Overlap.Corrected}, Areas, [1 2], fb + 6, ...
            [pre K.Headline.Corrected], K.Subtitle.Corrected); %#ok<AGROW>
        fnames{end+1} = ['BaselineCorrected_' sg]; %#ok<AGROW>
        FB = K.FracBars; FB.Headline = [pre FB.Headline];
        figs(end+1) = PlotPresentationFractionBars(FB, fb + 7, Areas, '% of neurons'); %#ok<AGROW>
        fnames{end+1} = ['AreaProportions_' sg]; %#ok<AGROW>
        OD = K.ORDots; OD.Headline = [pre OD.Headline];
        figs(end+1) = PlotPresentationEffectDots(OD, fb + 8, Areas); %#ok<AGROW>
        fnames{end+1} = ['AreaOverlapOR_' sg]; %#ok<AGROW>
    end
end

%% ---------------- Save ----------------
if OutputFlag
    if ~exist(SaveDir, 'dir'), mkdir(SaveDir); end
    for k = 1:numel(figs)
        f = fullfile(SaveDir, fnames{k});
        exportgraphics(figs(k), [f '.png'], 'Resolution', 300);
        savefig(figs(k), [f '.fig']);
    end
    for s = 1:numel(Params.Subgroups)
        sg = Params.Subgroups{s};
        writetable(LatCmp.(sg).Table, fullfile(SaveDir, sprintf('Timing_%s.csv', sg)));
        if isfield(Cmp, sg)
            writetable(Cmp.(sg).Table, fullfile(SaveDir, sprintf('AreaComparison_%s.csv', sg)));
        end
    end
    save(fullfile(SaveDir, 'AreaComparison.mat'), 'Lat', 'LatCmp', 'Deck', 'Cmp', 'Params', 'RowArea', '-v7.3');
    fprintf('Saved %d figures + tables -> %s\n', numel(figs), SaveDir);
end
end

% ---------------------------------------------------------------------------
function Cols = LatencyCols(LA, LB, Tb, fld, Params)
% LatencyCols  Arrange one timing signal (object | location columns) for
% PlotLatencySlide; p values are looked up from the comparison table.
sigLb = struct('Corrected', 'Value (baseline-corrected)', 'RawNoAntic', 'Value (raw, no anticipation)', ...
               'Visual', 'Visual response', 'ValueMinusVisual', 'Value - visual');
tasks = {'Obj', 'Loc'}; taskLb = {'Object', 'Location'};
isRel = strcmp(fld, 'ValueMinusVisual');
for i = 1:2
    SA = LA.(tasks{i}).(fld); SB = LB.(tasks{i}).(fld);
    r  = strcmp(Tb.Task, taskLb{i}) & strcmp(Tb.Signal, sigLb.(fld));
    pp = Tb.pHolm(r & strcmp(Tb.Measure, 'Population half-max onset'));
    pl = Tb.pHolm(r & strcmp(Tb.Measure, 'Single-neuron latency (median)'));
    if isempty(pp), pp = NaN; end
    C.Title = sprintf('%s value task', taskLb{i});
    C.t = LA.(tasks{i}).t;
    if isRel, C.PopA = []; C.PopB = []; else, C.PopA = SA.Pop; C.PopB = SB.Pop; end
    if strcmp(fld, 'Visual'), C.YLabel = 'Normalized response'; else, C.YLabel = 'AUC - 0.5 (sign-aligned)'; end
    C.XLimCurve = Params.LatDispWin;
    C.LatA = SA.Lat; C.LatB = SB.Lat;
    if isRel, C.LatLabel = 'Value onset - visual onset (ms)'; C.XLimLat = [-100 400];
    else,     C.LatLabel = 'Onset latency (ms)';              C.XLimLat = Params.LatSearchWin;
    end
    C.pPop = pp; C.pLat = pl;
    Cols(i) = C; %#ok<AGROW>
end
end

function S = SubsetBatch(Data, idx)
% SubsetBatch  Rows idx of a DataHub batch struct (fields used by the loader).
S.BatchFileName  = Data.BatchFileName;
S.ResultFilePath = Data.ResultFilePath(idx);
S.FileName       = Data.FileName(idx);
S.ChannelNumber  = Data.ChannelNumber(idx, :);
S.FileName       = S.FileName(:);
S.ResultFilePath = S.ResultFilePath(:);
end

function PrintComparison(DeckSg, K, Areas, sg)
% PrintComparison  Console table of overlap vs chance and FEF vs SC tests.
cond = {'All', 'WithBase', 'WithoutBase', 'Corrected'};
fprintf('\n===== %s neurons =====\n', sg);
fprintf('%-5s %-12s %5s %5s %5s %5s %7s %6s %8s  %s\n', 'Area', 'Analysis', 'N', 'nObj', 'nLoc', 'Both', 'Chance', 'OR', 'pFisher', 'verdict');
for a = 1:numel(Areas)
    for c = 1:numel(cond)
        O = DeckSg.(Areas{a}).Overlap.(cond{c});
        fprintf('%-5s %-12s %5d %5d %5d %5d %7.1f %6.2f %8.3g  %s\n', Areas{a}, cond{c}, O.N, O.nObj, O.nLoc, ...
            O.Obs(2), O.Exp(2), O.OR, O.pFisher, O.Verdict);
    end
end
fprintf('%s vs %s (Holm-corrected):\n', Areas{:}); disp(K.Table);
end
