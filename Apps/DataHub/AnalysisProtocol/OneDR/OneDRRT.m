function OneDRRT(Data, StartTrial, EndTrial, ShowFigureFlag, OutputFlag)
%ONEDRRT Per-session behavior for the location value (one-direction reward)
%task: reaction time, peak velocity and post-saccade gaze duration for
%saccades to high vs low value locations.
%
% The first NumSwitchExclude trials after each value-location reversal are
% dropped (same convention as OneDRAnalysis). Summary is also computed
% separately for each target side (Behavior.SummaryBySide) so that value
% effects can be checked independent of direction.
%
% Inputs:  DataHub protocol signature (Data, StartTrial, EndTrial,
%          ShowFigureFlag, OutputFlag)
% Outputs: figure (if ShowFigureFlag); OutputData.OneDR_Behavior (or
%          ObjLocationValue_Behavior) in the per-neuron result file (if OutputFlag)
%
% Computation: ComputeValueSaccadeBehavior
% Plotting:    PlotValueSaccadeBehavior
% Saving:      SaveValueBehaviorOutput
%
% Xuefei Yu; refactored 2026-11-05

ECode;

Params = struct();
Params.NumSwitchExclude = 5; % match OneDRAnalysis nSwitchExclude

Behavior = ComputeValueSaccadeBehavior(Data, StartTrial, EndTrial, Params);

if ShowFigureFlag
    PlotValueSaccadeBehavior(Behavior, ['LocationValue_Behavior_', Data.FileName]);
end

if Data.TaskCode == OBJ_LOCATIONVALUE_TASK
    ProtocolMarker = 'ObjLocationValue_Behavior';
else
    ProtocolMarker = 'OneDR_Behavior';
end

if OutputFlag
    OutputFiles = SaveValueBehaviorOutput(Data, ProtocolMarker, Behavior);
    fprintf('Data Exported (%d channel files): %s ... %s\n', numel(OutputFiles), ...
        OutputFiles{1}, OutputFiles{end});
else
    disp('No data export');
end

end
