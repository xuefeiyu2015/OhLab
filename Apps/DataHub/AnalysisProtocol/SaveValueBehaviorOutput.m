function OutputFiles = SaveValueBehaviorOutput(Data, ProtocolMarker, Behavior)
%SAVEVALUEBEHAVIOROUTPUT Write value-task behavior results into the
%per-neuron result files Results/[Monkey]_[Date]_N#_C#.mat as
%OutputData.(ProtocolMarker), one file per channel in Data.SpikeChannel
%(batch mode passes a channel range, e.g. 4:42; behavior is identical for
%all channels of a session). Existing fields in each file are preserved;
%the ProtocolMarker field is replaced.
%
% Stored entry:
%   .Task, .TaskCode, .TrialType
%   .TrialTable  table, one row per trial: TrialIdx, IsGood, TargetAngle,
%                RT, PeakV, GazeDur, Amp, Artifact, LandedOutside, GazeCensored
%   .SessionMean struct: SessionMean.(RT|PeakV|GazeDur) with fields
%                All/Good/Bad mean, sem, n (NaN trials excluded)
%   .DataStamp   containers.Map:
%                per trial  : RT, PeakV, GazeDur, IsGood, TrialIdx (all trials)
%                             <Metric>_Good, <Metric>_Bad (split by value)
%                per session: <Metric>_Mean, <Metric>_Good_Mean, <Metric>_Bad_Mean
%                (RT_Good, RT_Bad, RT_Good_Mean, RT_Bad_Mean kept for the
%                legacy batch scripts)
%   .Behavior    full struct from ComputeValueSaccadeBehavior
%
% Inputs:
%   Data           : DataHub Data struct
%   ProtocolMarker : e.g. 'StableValueBehavior', 'OneDR_Behavior'
%   Behavior       : struct from ComputeValueSaccadeBehavior
% Output:
%   OutputFiles    : {nChannel x 1} full paths of the saved .mat files
%
% Xuefei Yu

%% Build entry
B = Behavior;
G = B.IsGood;
Metrics = {'RT', 'PeakV', 'GazeDur'};

% Per-trial table
Entry.TrialTable = table(B.TrialIdx, B.IsGood, B.TargetAngle, B.RT, B.PeakV, ...
    B.GazeDur, B.Amp, B.Artifact, B.LandedOutside, B.GazeCensored, ...
    'VariableNames', {'TrialIdx','IsGood','TargetAngle','RT','PeakV', ...
    'GazeDur','Amp','Artifact','LandedOutside','GazeCensored'});

% Session means and DataStamp
NameStr = {'TrialIdx', 'IsGood'};
DataLib = {B.TrialIdx, B.IsGood};
for m = 1:numel(Metrics)
    M = Metrics{m};
    V = B.(M);
    Entry.SessionMean.(M).All  = MeanSem(V);
    Entry.SessionMean.(M).Good = MeanSem(V(G));
    Entry.SessionMean.(M).Bad  = MeanSem(V(~G));
    NameStr = [NameStr, {M, [M '_Good'], [M '_Bad'], ...
        [M '_Mean'], [M '_Good_Mean'], [M '_Bad_Mean']}]; %#ok<AGROW>
    DataLib = [DataLib, {V, V(G), V(~G), ...
        Entry.SessionMean.(M).All.Mean, Entry.SessionMean.(M).Good.Mean, ...
        Entry.SessionMean.(M).Bad.Mean}]; %#ok<AGROW>
end

Entry.TrialType = 'GoodTrials';
Entry.Task      = Data.TaskType;
Entry.TaskCode  = Data.TaskCode;
Entry.DataStamp = containers.Map(NameStr, DataLib);
Entry.Behavior  = B;

%% Write to every channel's result file
OutputFiles = SaveProtocolEntry(Data, ProtocolMarker, Entry);
end

function S = MeanSem(V)
% Mean, sem and n of non-NaN values
V = V(~isnan(V));
S.Mean = mean(V);
S.Sem  = std(V) / sqrt(numel(V));
S.N    = numel(V);
end
