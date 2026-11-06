function ValueTaskPerformance(Data, StartTrial, EndTrial, ShowFigureFlag, OutputFlag)
%VALUETASKPERFORMANCE Per-session break-fixation rate (plus overall error
%rate) and anticipatory saccade rate
%for high vs low value, in the object value (stable value) and location
%value (OneDR) tasks.
%
% Needs both Good and All Data of the session. Whichever DataHub passes,
% the other is loaded from the packed file (Data.Path/<session>.mat).
% StartTrial/EndTrial are interpreted in the trial numbering of the data
% DataHub passed (Good Data by default).
%
% Inputs:  DataHub protocol signature (Data, StartTrial, EndTrial,
%          ShowFigureFlag, OutputFlag)
% Outputs: figure (if ShowFigureFlag); OutputData.StableValuePerformance /
%          OneDR_Performance / ObjLocationValue_Performance /
%          FlexibleValuePerformance in each channel's result file (if OutputFlag)
%
% Computation: ComputeValueTaskPerformance
% Plotting:    PlotValueTaskPerformance
% Saving:      SaveProtocolEntry
%
% Xuefei Yu, 2026-11-05

ECode;

%% Task-specific settings
Params = struct();
switch Data.TaskCode
    case ONEDR_TASK
        ProtocolMarker = 'OneDR_Performance';
        Params.NumSwitchExclude = 5;          % match OneDRAnalysis / OneDRRT
    case OBJ_LOCATIONVALUE_TASK
        ProtocolMarker = 'ObjLocationValue_Performance';
        Params.NumSwitchExclude = 5;
    case FLEXIBLE_VALUE_TASK
        ProtocolMarker = 'FlexibleValuePerformance';
    otherwise
        ProtocolMarker = 'StableValuePerformance';
end

%% Good + All Data of this session
[GoodData, AllData, PassedIsAll] = GetGoodAndAllData(Data);
if PassedIsAll
    % Range was given in All Data numbering: convert to Good numbering
    Correct = find(any(AllData.EventChannel == CORRECTCD, 2));
    InRange = Correct(Correct >= StartTrial & Correct <= EndTrial);
    GoodRange = [find(Correct == InRange(1)), find(Correct == InRange(end))];
else
    GoodRange = [StartTrial, EndTrial];
end

Perf = ComputeValueTaskPerformance(GoodData, AllData, GoodRange, Params);

fprintf('%s | break-fixation rate High %.1f%% (%d/%d), Low %.1f%% (%d/%d), p=%.2g | anticipation (RT<%dms) High %.1f%%, Low %.1f%%, p=%.2g\n', ...
    Data.FileName, 100*Perf.BreakFix.High.Rate, Perf.BreakFix.High.NBreak, Perf.BreakFix.High.N, ...
    100*Perf.BreakFix.Low.Rate, Perf.BreakFix.Low.NBreak, Perf.BreakFix.Low.N, Perf.BreakFix.P, ...
    Perf.Anticipation.ThresholdMs, 100*Perf.Anticipation.High.Rate, ...
    100*Perf.Anticipation.Low.Rate, Perf.Anticipation.P);

if ShowFigureFlag
    PlotValueTaskPerformance(Perf, ['ValuePerformance_', Data.FileName]);
end

if OutputFlag
    Entry.Task      = Data.TaskType;
    Entry.TaskCode  = Data.TaskCode;
    Entry.Perf      = rmfield(Perf, 'Behavior');   % behavior is exported by RTGoodBadAnalysis / OneDR_RT
    Entry.Summary   = struct( ...
        'ErrorRate_High', Perf.Error.High.Rate, 'ErrorRate_Low', Perf.Error.Low.Rate, ...
        'ErrorN_High', Perf.Error.High.N, 'ErrorN_Low', Perf.Error.Low.N, 'Error_P', Perf.Error.P, ...
        'BreakFixRate_High', Perf.BreakFix.High.Rate, 'BreakFixRate_Low', Perf.BreakFix.Low.Rate, ...
        'BreakFixN_High', Perf.BreakFix.High.N, 'BreakFixN_Low', Perf.BreakFix.Low.N, ...
        'BreakFix_P', Perf.BreakFix.P, ...
        'AntRate_High', Perf.Anticipation.High.Rate, 'AntRate_Low', Perf.Anticipation.Low.Rate, ...
        'AntN_High', Perf.Anticipation.High.N, 'AntN_Low', Perf.Anticipation.Low.N, ...
        'Ant_P', Perf.Anticipation.P, ...
        'PreGoRate_High', Perf.Anticipation.High.RatePreGo, ...
        'PreGoRate_Low', Perf.Anticipation.Low.RatePreGo, ...
        'NumTargetNoValue', Perf.NumTargetNoValue, 'NumValueInferred', Perf.NumValueInferred);
    OutputFiles = SaveProtocolEntry(Data, ProtocolMarker, Entry);
    fprintf('Data Exported (%d channel files): %s ... %s\n', numel(OutputFiles), ...
        OutputFiles{1}, OutputFiles{end});
else
    disp('No data export');
end
end

%% ======================= Local functions =======================
function [GoodData, AllData, PassedIsAll] = GetGoodAndAllData(Data)
% Return GoodData and AllData for this session. The one DataHub passed is
% used as is; the other is loaded from the packed file.
PassedIsAll = isfield(Data, 'DataType') && any(strcmp(Data.DataType, 'All Data'));
Need = 'AllData';
if PassedIsAll
    Need = 'GoodData';
end
PackedFile = '';
if isfield(Data, 'Path') && ~isempty(Data.Path)
    PackedFile = fullfile(Data.Path, [erase(Data.FileName, '.rst'), '.mat']);
end
if isempty(PackedFile) || ~exist(PackedFile, 'file')
    error('ValueTaskPerformance:NoPackedFile', ...
        'Cannot find packed file for %s (Data.Path missing or file not found): %s', ...
        Data.FileName, PackedFile);
end
L = load(PackedFile, Need);
if PassedIsAll
    AllData = Data;  GoodData = L.GoodData;
else
    GoodData = Data; AllData = L.AllData;
end
end
