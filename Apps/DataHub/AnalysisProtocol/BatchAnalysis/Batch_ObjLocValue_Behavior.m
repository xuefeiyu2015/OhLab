function Batch_ObjLocValue_Behavior(Data, StartFile, EndFile, ShowFigureFlag, OutputFlag)
%BATCH_OBJLOCVALUE_BEHAVIOR Population comparison of saccade behavior for
%high vs low value in the object value (stable value) and location value
%(OneDR) tasks: reaction time, peak velocity, gaze duration.
%
% Load with a stable value batch file (e.g. FEF_BothRA_StableValue.txt).
% Reads OutputData.StableValueBehavior and OutputData.OneDR_Behavior saved
% by RTGoodBadAnalysis.m / OneDRRT.m (SaveValueBehaviorOutput.m format), and,
% if present, OutputData.StableValuePerformance / OneDR_Performance saved by
% ValueTaskPerformance.m (error or break-fixation rate, anticipatory saccades;
% fig 208; choose with FailureMeasure at the top).
%
% Unit = session. All channel files of one session carry identical
% behavior, so each batch-file row contributes once (first channel file
% that has BOTH tasks exported in the current format). Rows without both
% are skipped and listed.
%
% Inputs:  DataHub batch signature (Data, StartFile, EndFile,
%          ShowFigureFlag, OutputFlag)
% Outputs: figures (if ShowFigureFlag);
%          Results/ObjLocBehavior/<Area>/ObjLocBehavior_<BatchFile>.mat (if OutputFlag)
%
% Xuefei Yu, 2026-11-05

BatchFileName = Data.BatchFileName;
FilesName     = Data.ResultFilePath(StartFile:EndFile);
SessionNames  = Data.FileName(StartFile:EndFile);

Area = GetAreaLabel(BatchFileName);

% Failure measure for fig 208 / console: 'Error' (any uncompleted trial
% with the target shown) or 'BreakFix' (fixation broken before go)
FailureMeasure = 'Error';

%% Load (I/O)
Sessions = LoadSessionBehavior(FilesName, SessionNames);
fprintf('%s: %d / %d sessions have both object and location value behavior\n', ...
    BatchFileName, numel(Sessions.Name), numel(FilesName));
if ~isempty(Sessions.Skipped)
    fprintf('Skipped (missing or old-format export):\n');
    fprintf('   %s\n', Sessions.Skipped{:});
end
if isempty(Sessions.Name)
    disp('No session to analyze. Re-run RTGoodBadAnalysis and OneDR_RT with export on.');
    return
end

%% Compute
Stats = ComputeBatchValueBehavior(Sessions, FailureMeasure);
PrintStats(Stats);
PrintPerfStats(Stats);

%% Plot
if ShowFigureFlag
    TitleStr = sprintf('%s (%s, n=%d sessions)', strrep(BatchFileName, '_', ' '), Area, Stats.nSession);
    PlotPairedHighLow(Stats, TitleStr);
    PlotEffectObjVsLoc(Stats, TitleStr);
    for m = 1:numel(Stats.Metrics)
        PlotFeatureHistogram(Stats, m, TitleStr, 202 + m); % figs 203-205
    end
    PlotFeatureSpace3D(Stats, TitleStr); % fig 206
    PlotViolin(Stats, TitleStr); % fig 207
    PlotPerformance(Stats, TitleStr); % fig 208
end

%% Save
if OutputFlag
    ResultRoot = fileparts(FilesName{1}{1});
    SaveDir = fullfile(ResultRoot, 'ObjLocBehavior', Area);
    if ~exist(SaveDir, 'dir')
        mkdir(SaveDir);
    end
    SaveFile = fullfile(SaveDir, sprintf('ObjLocBehavior_%s.mat', BatchFileName));
    save(SaveFile, 'Sessions', 'Stats');
    fprintf('Saved: %s\n', SaveFile);
end

end

%% ======================= Load =======================
function Sessions = LoadSessionBehavior(FilesName, SessionNames)
% Collect per-session tables and session means for both tasks.
% Output struct:
%   .Name        {nSession x 1} session (raw file) names
%   .ResultFile  {nSession x 1} result file used
%   .Obj, .Loc   struct arrays with fields TrialTable, SessionMean, Summary,
%                Perf (ValueTaskPerformance .Summary, [] if not exported)
%   .Skipped     {k x 1} sessions without both exports
Tasks  = {'Obj', 'Loc'};
Fields = {'StableValueBehavior', 'OneDR_Behavior'};
PerfFields = {'StableValuePerformance', 'OneDR_Performance'};

Sessions.Name = {}; Sessions.ResultFile = {}; Sessions.Skipped = {};
Sessions.Obj = struct('TrialTable', {}, 'SessionMean', {}, 'Summary', {}, 'Perf', {});
Sessions.Loc = Sessions.Obj;

for i = 1:numel(FilesName)
    Found = false;
    for j = 1:numel(FilesName{i})
        f = FilesName{i}{j};
        if ~exist(f, 'file')
            continue
        end
        L = load(f, 'OutputData');
        if ~isfield(L, 'OutputData') || ~all(isfield(L.OutputData, Fields))
            continue
        end
        Entries = cellfun(@(x) L.OutputData.(x), Fields, 'UniformOutput', false);
        if ~all(cellfun(@(e) isfield(e, 'TrialTable') && isfield(e, 'SessionMean'), Entries))
            continue % old-format export
        end
        n = numel(Sessions.Name) + 1;
        Sessions.Name{n,1} = SessionNames{i};
        Sessions.ResultFile{n,1} = f;
        for t = 1:2
            Sessions.(Tasks{t})(n).TrialTable  = Entries{t}.TrialTable;
            Sessions.(Tasks{t})(n).SessionMean = Entries{t}.SessionMean;
            Sessions.(Tasks{t})(n).Summary     = Entries{t}.Behavior.Summary;
            Sessions.(Tasks{t})(n).Perf        = [];
            if isfield(L.OutputData, PerfFields{t})
                Sessions.(Tasks{t})(n).Perf = L.OutputData.(PerfFields{t}).Summary;
            end
        end
        Found = true;
        break
    end
    if ~Found
        Sessions.Skipped{end+1,1} = SessionNames{i};
    end
end
end

%% ======================= Compute =======================
function Stats = ComputeBatchValueBehavior(Sessions, FailureMeasure)
% Pure computation on loaded sessions.
% Stats.(Task).(Metric):
%   .Good, .Bad   [nSession x 1] session means
%   .Diff         Good - Bad
%   .SessionP     within-session rank-sum p
%   .MeanGood/.SemGood/.MeanBad/.SemBad, .P (signrank, paired across sessions)
%   .PooledGood, .PooledBad  all trials pooled (for reference)
% Stats.Compare.(Metric): obj vs loc value effect (Diff)
%   .P (signrank on Diff_obj vs Diff_loc), .R/.PCorr (Spearman)
% Stats.Trials.(Task): table of all valid trials pooled across sessions
%   (Monkey, Session, IsGood, RT, PeakV, GazeDur); trials with any
%   missing metric removed
% Stats.Monkeys: monkey names present ('Robin', 'Adams')
Tasks   = {'Obj', 'Loc'};
Metrics = {'RT', 'PeakV', 'GazeDur'};
nS = numel(Sessions.Name);

Stats.Tasks = Tasks;
Stats.TaskLabels = {'Object value', 'Location value'};
Stats.Metrics = Metrics;
Stats.MetricLabels = {'Reaction time (ms)', 'Peak velocity (deg/s)', 'Gaze duration (ms)'};
Stats.nSession = nS;
Stats.SessionNames = Sessions.Name;

for t = 1:numel(Tasks)
    S = Sessions.(Tasks{t});
    for m = 1:numel(Metrics)
        M = Metrics{m};
        Good = arrayfun(@(x) x.SessionMean.(M).Good.Mean, S)';
        Bad  = arrayfun(@(x) x.SessionMean.(M).Bad.Mean,  S)';
        R.Good = Good;
        R.Bad  = Bad;
        R.Diff = Good - Bad;
        R.SessionP = arrayfun(@(x) x.Summary.(M).P, S)';
        R.MeanGood = mean(Good, 'omitnan'); R.SemGood = SemOmit(Good);
        R.MeanBad  = mean(Bad,  'omitnan'); R.SemBad  = SemOmit(Bad);
        Valid = ~isnan(Good) & ~isnan(Bad);
        R.N = sum(Valid);
        if R.N > 1
            R.P = signrank(Good(Valid), Bad(Valid));
        else
            R.P = NaN;
        end
        PooledG = cell(nS,1); PooledB = cell(nS,1);
        for s = 1:nS
            T = S(s).TrialTable;
            PooledG{s} = T.(M)(T.IsGood);
            PooledB{s} = T.(M)(~T.IsGood);
        end
        R.PooledGood = cell2mat(PooledG);
        R.PooledBad  = cell2mat(PooledB);
        Stats.(Tasks{t}).(M) = R;
    end
end

% Pooled trial tables with monkey label (session name starting with 'A' = Adams)
MonkeyOfSession = repmat({'Robin'}, nS, 1);
MonkeyOfSession(startsWith(Sessions.Name, 'A')) = {'Adams'};
Stats.MonkeyOfSession = MonkeyOfSession;
Stats.Monkeys = intersect({'Robin', 'Adams'}, MonkeyOfSession, 'stable');
for t = 1:numel(Tasks)
    Parts = cell(nS, 1);
    for s = 1:nS
        T = Sessions.(Tasks{t})(s).TrialTable(:, [{'IsGood'}, Metrics]);
        T.Monkey  = repmat(MonkeyOfSession(s), height(T), 1);
        T.Session = repmat(Sessions.Name(s), height(T), 1);
        Parts{s} = T;
    end
    T = vertcat(Parts{:});
    Stats.Trials.(Tasks{t}) = T(all(~isnan(T{:, Metrics}), 2), :);
end

% Failure rate (FailureMeasure: 'Error' or 'BreakFix') and anticipation
% (ValueTaskPerformance). Per task, measure and
% monkey: session rates High/Low (NaN if not exported), paired signrank
% across sessions, and Fisher test on counts pooled over sessions.
PerfMeasures = {FailureMeasure, 'Ant'};
Stats.PerfMeasures = PerfMeasures;
FailLabels = struct('Error', 'Error rate (%)', 'BreakFix', 'Break-fixation rate (%)');
Stats.PerfLabels = {FailLabels.(FailureMeasure), 'Anticipatory saccades (%)'};
for t = 1:numel(Tasks)
    S = Sessions.(Tasks{t});
    HasPerf = arrayfun(@(x) ~isempty(x.Perf), S)';
    for q = 1:numel(PerfMeasures)
        Q = PerfMeasures{q};
        RateH = nan(nS,1); RateL = nan(nS,1); NH = zeros(nS,1); NL = zeros(nS,1);
        for s = find(HasPerf)'
            Pf = S(s).Perf;
            if any(strcmp(Q, {'Error', 'BreakFix'}))
                if ~isfield(Pf, [Q 'Rate_High'])
                    continue % exported before this measure existed: re-run ValueTaskPerformance
                end
                RateH(s) = Pf.([Q 'Rate_High']); RateL(s) = Pf.([Q 'Rate_Low']);
                NH(s) = Pf.([Q 'N_High']);       NL(s) = Pf.([Q 'N_Low']);
            else
                RateH(s) = Pf.AntRate_High;   RateL(s) = Pf.AntRate_Low;
                NH(s) = Pf.AntN_High;         NL(s) = Pf.AntN_Low;
            end
        end
        R = struct('High', RateH, 'Low', RateL, 'NHigh', NH, 'NLow', NL, 'HasPerf', HasPerf);
        for k = 1:numel(Stats.Monkeys)
            Mk = strcmp(MonkeyOfSession, Stats.Monkeys{k}) & HasPerf & ~isnan(RateH) & ~isnan(RateL);
            Z.N = sum(Mk);
            Z.MeanHigh = mean(RateH(Mk)); Z.SemHigh = std(RateH(Mk)) / sqrt(Z.N);
            Z.MeanLow  = mean(RateL(Mk)); Z.SemLow  = std(RateL(Mk)) / sqrt(Z.N);
            Z.P = NaN;
            if Z.N > 1
                Z.P = signrank(RateH(Mk), RateL(Mk));
            end
            kH = round(sum(RateH(Mk) .* NH(Mk))); nH = sum(NH(Mk));
            kL = round(sum(RateL(Mk) .* NL(Mk))); nL = sum(NL(Mk));
            Z.PooledHigh = kH / nH; Z.PooledLow = kL / nL;
            Z.PPooled = NaN;
            if nH > 0 && nL > 0
                [~, Z.PPooled] = fishertest([kH, nH - kH; kL, nL - kL]);
            end
            R.(Stats.Monkeys{k}) = Z;
        end
        Stats.Perf.(Tasks{t}).(Q) = R;
    end
end

for m = 1:numel(Metrics)
    M = Metrics{m};
    DObj = Stats.Obj.(M).Diff;
    DLoc = Stats.Loc.(M).Diff;
    Valid = ~isnan(DObj) & ~isnan(DLoc);
    C.N = sum(Valid);
    if C.N > 2
        C.P = signrank(DObj(Valid), DLoc(Valid));
        [C.R, C.PCorr] = corr(DObj(Valid), DLoc(Valid), 'Type', 'Spearman');
    else
        C.P = NaN; C.R = NaN; C.PCorr = NaN;
    end
    Stats.Compare.(M) = C;
end
end

function s = SemOmit(x)
x = x(~isnan(x));
s = std(x) / sqrt(numel(x));
end

function Area = GetAreaLabel(BatchFileName)
if contains(BatchFileName, 'FEFSC')
    Area = 'FEFSC';
elseif contains(BatchFileName, 'FEF')
    Area = 'FEF';
elseif contains(BatchFileName, 'SC')
    Area = 'SC';
else
    Area = 'Unknown';
end
end

%% ======================= Report =======================
function PrintStats(Stats)
fprintf('\n%-16s %-8s %18s %18s %8s\n', 'Task', 'Metric', 'High (mean±sem)', 'Low (mean±sem)', 'p');
for t = 1:numel(Stats.Tasks)
    for m = 1:numel(Stats.Metrics)
        R = Stats.(Stats.Tasks{t}).(Stats.Metrics{m});
        fprintf('%-16s %-8s %10.1f ± %-5.1f %10.1f ± %-5.1f %8.3g\n', Stats.TaskLabels{t}, ...
            Stats.Metrics{m}, R.MeanGood, R.SemGood, R.MeanBad, R.SemBad, R.P);
    end
end
end

%% ======================= Plot =======================
function Fig = PlotPairedHighLow(Stats, TitleStr)
% Rows: task; columns: metric. Each line = one session (high vs low mean).
% Filled markers: session with significant within-session difference (p<0.05).
Fig = PrepareFigure(201, 'w', [50, 80, 1200, 800], 'Name', 'ObjLoc_Behavior_HighLow');
nT = numel(Stats.Tasks); nM = numel(Stats.Metrics);
for t = 1:nT
    [ColorGood, ColorBad] = ValueColors(Stats.Tasks{t}); % obj red/blue, loc purple/green
    for m = 1:nM
        R  = Stats.(Stats.Tasks{t}).(Stats.Metrics{m});
        ax = subplot(nT, nM, (t-1)*nM + m); hold(ax, 'on');
        Sig = R.SessionP < 0.05;
        plot(ax, [1 2], [R.Good R.Bad]', '-', 'Color', [0.7 0.7 0.7]);
        scatter(ax, ones(sum(~Sig),1),   R.Good(~Sig), 25, ColorGood);
        scatter(ax, 2*ones(sum(~Sig),1), R.Bad(~Sig),  25, ColorBad);
        scatter(ax, ones(sum(Sig),1),    R.Good(Sig),  25, ColorGood, 'filled');
        scatter(ax, 2*ones(sum(Sig),1),  R.Bad(Sig),   25, ColorBad,  'filled');
        errorbar(ax, [0.75 2.25], [R.MeanGood R.MeanBad], [R.SemGood R.SemBad], 'ks-', ...
            'MarkerFaceColor', 'k', 'LineWidth', 1.5);
        set(ax, 'XTick', [1 2], 'XTickLabel', {'High', 'Low'}, 'XLim', [0.5 2.5], 'TickDir', 'out');
        box(ax, 'off');
        ylabel(ax, Stats.MetricLabels{m});
        title(ax, sprintf('%s: p=%.2g (n=%d)', Stats.TaskLabels{t}, R.P, R.N));
    end
end
sgtitle(Fig, {TitleStr, 'lines: sessions; filled: within-session p<0.05; black: mean ± sem'});
end

function Fig = PlotEffectObjVsLoc(Stats, TitleStr)
% Value effect (High - Low) per session: object value (x) vs location value (y).
Fig = PrepareFigure(202, 'w', [80, 120, 1200, 420], 'Name', 'ObjLoc_Behavior_EffectCompare');
for m = 1:numel(Stats.Metrics)
    M  = Stats.Metrics{m};
    x  = Stats.Obj.(M).Diff;
    y  = Stats.Loc.(M).Diff;
    C  = Stats.Compare.(M);
    ax = subplot(1, numel(Stats.Metrics), m); hold(ax, 'on');
    Lim = max(abs([x; y]), [], 'omitnan') * 1.1;
    if isempty(Lim) || isnan(Lim) || Lim == 0
        Lim = 1;
    end
    plot(ax, [-Lim Lim], [0 0], 'k:'); plot(ax, [0 0], [-Lim Lim], 'k:');
    plot(ax, [-Lim Lim], [-Lim Lim], '--', 'Color', [0.6 0.6 0.6]);
    scatter(ax, x, y, 30, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.7);
    axis(ax, 'square'); set(ax, 'XLim', [-Lim Lim], 'YLim', [-Lim Lim], 'TickDir', 'out');
    xlabel(ax, 'Object value: High - Low');
    ylabel(ax, 'Location value: High - Low');
    title(ax, {Stats.MetricLabels{m}, sprintf('obj vs loc p=%.2g; rho=%.2f (p=%.2g)', C.P, C.R, C.PCorr)});
end
sgtitle(Fig, TitleStr);
end

function Fig = PlotFeatureHistogram(Stats, m, TitleStr, FigNum)
% Distribution of one metric (Stats.Metrics{m}), high vs low value.
% Rows: monkey; columns: task. Dashed lines: mean of each distribution.
% RT uses fixed bins (-150:20:500 ms re: FP off); other metrics use
% 25-unit bins spanning each panel's data range. Inset (top right): bar
% plot of high vs low mean +/- sem (same trials).
Metric = Stats.Metrics{m};
Label  = Stats.MetricLabels{m};
Unit   = regexp(Label, '\((.*)\)', 'tokens', 'once');
nMk = numel(Stats.Monkeys); nT = numel(Stats.Tasks);
Fig = PrepareFigure(FigNum, 'w', [100 + 20*m, 100, 1000, 380*nMk], ...
    'Name', ['ObjLoc_Behavior_Hist_' Metric]);
MainAx = gobjects(nMk, nT);
BarMean = nan(nMk, nT, 2); BarSem = nan(nMk, nT, 2);
for k = 1:nMk
    for t = 1:nT
        [ColorGood, ColorBad] = ValueColors(Stats.Tasks{t});
        T  = Stats.Trials.(Stats.Tasks{t});
        T  = T(strcmp(T.Monkey, Stats.Monkeys{k}), :);
        G  = T.(Metric)(T.IsGood);
        B  = T.(Metric)(~T.IsGood);
        if strcmp(Metric, 'RT')
            Edges = -150:20:500;
        else
            BinW = 25;
            All  = [G; B];
            if isempty(All), All = 0; end
            Edges = floor(min(All)/BinW)*BinW : BinW : ceil(max(All)/BinW)*BinW + BinW;
        end
        ax = subplot(nMk, nT, (k-1)*nT + t); hold(ax, 'on');
        histogram(ax, G, Edges, 'Normalization', 'probability', 'FaceColor', ColorGood, ...
            'FaceAlpha', 0.4, 'EdgeColor', 'none');
        histogram(ax, B, Edges, 'Normalization', 'probability', 'FaceColor', ColorBad, ...
            'FaceAlpha', 0.4, 'EdgeColor', 'none');
        YL = ylim(ax);
        plot(ax, mean(G) * [1 1], YL, '--', 'Color', ColorGood, 'LineWidth', 1.5);
        plot(ax, mean(B) * [1 1], YL, '--', 'Color', ColorBad,  'LineWidth', 1.5);
        BarMean(k, t, :) = [mean(G), mean(B)];
        BarSem(k, t, :)  = [std(G) / sqrt(numel(G)), std(B) / sqrt(numel(B))];
        if numel(G) > 1 && numel(B) > 1
            P = ranksum(G, B);
        else
            P = NaN;
        end
        title(ax, {sprintf('%s - %s', Stats.Monkeys{k}, Stats.TaskLabels{t}), ...
            sprintf('mean High %.0f / Low %.0f %s, p=%.2g', mean(G), mean(B), Unit{1}, P)});
        xlabel(ax, Label); ylabel(ax, 'Fraction of trials');
        legend(ax, {sprintf('High (n=%d)', numel(G)), sprintf('Low (n=%d)', numel(B))}, ...
            'Box', 'off', 'Location', 'northwest');
        set(ax, 'TickDir', 'out'); box(ax, 'off');
        % Headroom at the top for the inset; no ticks in the headroom
        Tk = yticks(ax);
        ylim(ax, [YL(1), YL(2) + 0.6 * diff(YL)]);
        yticks(ax, Tk(Tk <= YL(2)));
        MainAx(k, t) = ax;
    end
end
sgtitle(Fig, {TitleStr, ['inset: mean ', char(177), ' sem']}, 'Interpreter', 'none');

% Insets: drawn after sgtitle so panel positions are final
for k = 1:nMk
    for t = 1:nT
        [ColorGood, ColorBad] = ValueColors(Stats.Tasks{t});
        AxPos = MainAx(k, t).Position;
        In = axes(Fig, 'Position', [AxPos(1) + 0.70*AxPos(3), AxPos(2) + 0.60*AxPos(4), ...
            0.26*AxPos(3), 0.27*AxPos(4)]);
        hold(In, 'on');
        Mu = squeeze(BarMean(k, t, :))'; Se = squeeze(BarSem(k, t, :))';
        b = bar(In, [1 2], Mu, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'none');
        b.CData = [ColorGood; ColorBad];
        errorbar(In, [1 2], Mu, Se, 'k.', 'LineWidth', 1, 'CapSize', 3);
        Lo = min([0, Mu - Se]); Hi = max([0, Mu + Se]);
        set(In, 'XTick', [1 2], 'XTickLabel', {'High', 'Low'}, 'XLim', [0.4 2.6], ...
            'YLim', [Lo, Hi] + [-0.05 0.15] * (Hi - Lo + eps), ...
            'TickDir', 'out', 'FontSize', 7, 'Color', 'none');
        ylabel(In, Unit{1}, 'FontSize', 7);
        box(In, 'off');
    end
end
end

function Fig = PlotFeatureSpace3D(Stats, TitleStr)
% Trials in (RT, peak velocity, gaze duration) space, high vs low value.
% Rows: monkey; columns: task. Small dots: trials; large markers: means.
nMk = numel(Stats.Monkeys); nT = numel(Stats.Tasks);
Fig = PrepareFigure(206, 'w', [120, 80, 1100, 480*nMk], 'Name', 'ObjLoc_Behavior_3D');
for k = 1:nMk
    for t = 1:nT
        [ColorGood, ColorBad] = ValueColors(Stats.Tasks{t});
        T  = Stats.Trials.(Stats.Tasks{t});
        T  = T(strcmp(T.Monkey, Stats.Monkeys{k}), :);
        G  = T{T.IsGood,  {'RT', 'PeakV', 'GazeDur'}};
        B  = T{~T.IsGood, {'RT', 'PeakV', 'GazeDur'}};
        ax = subplot(nMk, nT, (k-1)*nT + t); hold(ax, 'on');
        scatter3(ax, G(:,1), G(:,2), G(:,3), 8, ColorGood, 'filled', 'MarkerFaceAlpha', 0.35);
        scatter3(ax, B(:,1), B(:,2), B(:,3), 8, ColorBad,  'filled', 'MarkerFaceAlpha', 0.35);
        mG = mean(G, 1); mB = mean(B, 1);
        scatter3(ax, mG(1), mG(2), mG(3), 140, ColorGood, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        scatter3(ax, mB(1), mB(2), mB(3), 140, ColorBad,  'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        plot3(ax, [mG(1) mB(1)], [mG(2) mB(2)], [mG(3) mB(3)], 'k-', 'LineWidth', 1.5);
        xlabel(ax, 'RT (ms)'); ylabel(ax, 'Peak velocity (deg/s)'); zlabel(ax, 'Gaze duration (ms)');
        title(ax, sprintf('%s - %s (High n=%d, Low n=%d)', Stats.Monkeys{k}, ...
            Stats.TaskLabels{t}, size(G,1), size(B,1)));
        legend(ax, {'High', 'Low'}, 'Box', 'off', 'Location', 'northeast');
        grid(ax, 'on'); view(ax, 45, 25); set(ax, 'TickDir', 'out');
    end
end
sgtitle(Fig, {TitleStr, 'dots: trials; large markers: mean'});
end

function Fig = PlotViolin(Stats, TitleStr)
% Violin plots of trial-level saccade parameters, high vs low value, for
% both tasks. Rows: monkey; columns: metric. Each panel: object high/low,
% location high/low. White dot: median; black bar: interquartile range.
% Title: rank-sum p (high vs low) per task.
Pos   = [1 2 4 5];                         % obj high, obj low, loc high, loc low
nMk = numel(Stats.Monkeys); nM = numel(Stats.Metrics);
Fig = PrepareFigure(207, 'w', [140, 60, 1300, 420*nMk], 'Name', 'ObjLoc_Behavior_Violin');
for k = 1:nMk
    for m = 1:nM
        M  = Stats.Metrics{m};
        ax = subplot(nMk, nM, (k-1)*nM + m); hold(ax, 'on');
        P  = nan(1, numel(Stats.Tasks));
        for t = 1:numel(Stats.Tasks)
            [ColorGood, ColorBad] = ValueColors(Stats.Tasks{t});
            T = Stats.Trials.(Stats.Tasks{t});
            T = T(strcmp(T.Monkey, Stats.Monkeys{k}), :);
            Vals = {T.(M)(T.IsGood), T.(M)(~T.IsGood)};
            Cols = {ColorGood, ColorBad};
            for g = 1:2
                V = Vals{g};
                x = Pos(2*(t-1) + g);
                if numel(V) > 1
                    violinplot(ax, x * ones(size(V)), V, 'FaceColor', Cols{g}, ...
                        'FaceAlpha', 0.5, 'EdgeColor', 'none');
                end
                if ~isempty(V)
                    Q = prctile(V, [25 50 75]);
                    plot(ax, [x x], Q([1 3]), 'k-', 'LineWidth', 3);
                    plot(ax, x, Q(2), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
                end
            end
            if numel(Vals{1}) > 1 && numel(Vals{2}) > 1
                P(t) = ranksum(Vals{1}, Vals{2});
            end
        end
        set(ax, 'XTick', Pos, 'XTickLabel', {'High\newlineObj', 'Low\newlineObj', ...
            'High\newlineLoc', 'Low\newlineLoc'}, ...
            'TickLabelInterpreter', 'tex', 'XLim', [0 6], 'TickDir', 'out');
        box(ax, 'off');
        ylabel(ax, Stats.MetricLabels{m});
        title(ax, {sprintf('%s', Stats.Monkeys{k}), ...
            sprintf('Object p=%.2g | Location p=%.2g', P(1), P(2))});
    end
end
sgtitle(Fig, {TitleStr, 'trials pooled across sessions; white dot: median; bar: IQR'});
end

function PrintPerfStats(Stats)
% Console summary of error rate / anticipation per monkey and task
Fm = Stats.PerfMeasures{1};
nHas = sum(Stats.Perf.Obj.(Fm).HasPerf & Stats.Perf.Loc.(Fm).HasPerf);
fprintf('\n%s / anticipation (ValueTaskPerformance exported for both tasks in %d / %d sessions)\n', ...
    Stats.PerfLabels{1}, ...
    nHas, Stats.nSession);
for k = 1:numel(Stats.Monkeys)
    for t = 1:numel(Stats.Tasks)
        for q = 1:numel(Stats.PerfMeasures)
            Z = Stats.Perf.(Stats.Tasks{t}).(Stats.PerfMeasures{q}).(Stats.Monkeys{k});
            fprintf('%-6s %-15s %-6s High %5.1f%% Low %5.1f%% (n=%d sessions, signrank p=%.2g; pooled Fisher p=%.2g)\n', ...
                Stats.Monkeys{k}, Stats.TaskLabels{t}, Stats.PerfMeasures{q}, ...
                100*Z.MeanHigh, 100*Z.MeanLow, Z.N, Z.P, Z.PPooled);
        end
    end
end
end

function Fig = PlotPerformance(Stats, TitleStr)
% Failure rate (error or break-fixation, Stats.PerfMeasures{1}) and
% anticipatory saccade rate, high vs low value.
% Rows: monkey; columns: Failure-Object, Failure-Location, Anticipation-Object,
% Anticipation-Location. Bars (full colour): mean +/- sem across sessions;
% white dots: sessions (paired by grey lines); bracket: paired signrank across sessions
% (*** p<0.001, ** p<0.01, * p<0.05, n.s.). Title also gives Fisher on
% counts pooled over sessions.
nMk = numel(Stats.Monkeys);
Fm = Stats.PerfMeasures{1};
Cols = {{'Obj', Fm}, {'Loc', Fm}, {'Obj', 'Ant'}, {'Loc', 'Ant'}};
Fig = PrepareFigure(208, 'w', [160, 80, 1400, 420*nMk], 'Name', 'ObjLoc_Behavior_Performance');
for k = 1:nMk
    Mk = strcmp(Stats.MonkeyOfSession, Stats.Monkeys{k});
    for c = 1:numel(Cols)
        Task = Cols{c}{1}; Q = Cols{c}{2};
        [ColorGood, ColorBad] = ValueColors(Task);
        R  = Stats.Perf.(Task).(Q);
        Z  = R.(Stats.Monkeys{k});
        Sel = Mk & R.HasPerf & ~isnan(R.High) & ~isnan(R.Low);
        ax = subplot(nMk, numel(Cols), (k-1)*numel(Cols) + c); hold(ax, 'on');
        H = 100 * R.High(Sel); L = 100 * R.Low(Sel);
        Top = 1;
        if any(Sel)
            % Bars (full colour): mean across sessions
            b = bar(ax, [1 2], 100 * [Z.MeanHigh Z.MeanLow], 0.6, 'FaceColor', 'flat', ...
                'EdgeColor', 'none');
            b.CData = [ColorGood; ColorBad];
            % Session dots on top (white, black edge), paired; small jitter
            Rs = RandStream('mt19937ar', 'Seed', k*10 + c);
            Jx = 0.08 * (rand(Rs, numel(H), 1) - 0.5) * 2;
            plot(ax, [1 + Jx, 2 + Jx]', [H L]', '-', 'Color', [0.35 0.35 0.35 0.6]);
            scatter(ax, [1 + Jx; 2 + Jx], [H; L], 22, 'w', 'filled', 'MarkerEdgeColor', 'k');
            % Sem error bars last, so they stay visible
            errorbar(ax, [1 2], 100 * [Z.MeanHigh Z.MeanLow], 100 * [Z.SemHigh Z.SemLow], ...
                'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
            Top = max([H; L; 100 * [Z.MeanHigh + Z.SemHigh; Z.MeanLow + Z.SemLow]; 1]);
            % Significance bracket
            Yb = Top * 1.08; Tick = Top * 0.03;
            plot(ax, [1 1 2 2], [Yb - Tick, Yb, Yb, Yb - Tick], 'k-', 'LineWidth', 1);
            text(ax, 1.5, Yb, SigStars(Z.P), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
        end
        set(ax, 'XTick', [1 2], 'XTickLabel', {'High', 'Low'}, 'XLim', [0.4 2.6], 'TickDir', 'out');
        ylim(ax, [0, Top * 1.25]);                    % rates: axis from 0
        box(ax, 'off');
        ylabel(ax, Stats.PerfLabels{strcmp(Stats.PerfMeasures, Q)});
        TaskLabel = Stats.TaskLabels{strcmp(Stats.Tasks, Task)};
        title(ax, {sprintf('%s - %s', Stats.Monkeys{k}, TaskLabel), ...
            sprintf('signrank p=%.2g (n=%d) | pooled p=%.2g', Z.P, Z.N, Z.PPooled)});
    end
end
sgtitle(Fig, {TitleStr, ['bars: mean ', char(177), ' sem across sessions; white dots: sessions; ', ...
    'bracket: paired signrank (* p<0.05, ** p<0.01, *** p<0.001)']}, 'Interpreter', 'none');
end

function S = SigStars(p)
% Significance label for a p-value
if isnan(p)
    S = 'n/a';
elseif p < 0.001
    S = '***';
elseif p < 0.01
    S = '**';
elseif p < 0.05
    S = '*';
else
    S = 'n.s.';
end
end

