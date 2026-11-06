function RTGoodBadAnalysis(Data, StartTrial, EndTrial, ShowFigureFlag, OutputFlag)
%RTGOODBADANALYSIS Per-session behavior for the object value (stable value)
%task: reaction time, peak velocity and post-saccade gaze duration for
%saccades to high vs low value objects.
%
% Inputs:  DataHub protocol signature (Data, StartTrial, EndTrial,
%          ShowFigureFlag, OutputFlag)
% Outputs: figure (if ShowFigureFlag); OutputData.StableValueBehavior (or
%          FlexibleValueBehavior) in the per-neuron result file (if OutputFlag)
%
% Computation: ComputeValueSaccadeBehavior
% Plotting:    PlotValueSaccadeBehavior
% Saving:      SaveValueBehaviorOutput
%
% Xuefei Yu 05122020; refactored 2026-11-05

ECode;

Params = struct();
Params.NumSwitchExclude = 0; % stable object-value associations, no switches

Behavior = ComputeValueSaccadeBehavior(Data, StartTrial, EndTrial, Params);

if ShowFigureFlag
    PlotValueSaccadeBehavior(Behavior, ['ObjectValue_Behavior_', Data.FileName]);
end

if Data.TaskCode == FLEXIBLE_VALUE_TASK
    ProtocolMarker = 'FlexibleValueBehavior';
else
    ProtocolMarker = 'StableValueBehavior';
end

if OutputFlag
    OutputFiles = SaveValueBehaviorOutput(Data, ProtocolMarker, Behavior);
    fprintf('Data Exported (%d channel files): %s ... %s\n', numel(OutputFiles), ...
        OutputFiles{1}, OutputFiles{end});
else
    disp('No data export');
end

end
