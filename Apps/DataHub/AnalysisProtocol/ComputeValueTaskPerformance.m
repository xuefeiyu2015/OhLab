function Perf = ComputeValueTaskPerformance(GoodData, AllData, GoodRange, Params)
%COMPUTEVALUETASKPERFORMANCE Error rate and anticipatory saccade rate for
%high vs low value trials in one session (object value or OneDR task).
%
% Pure computation: no plotting, printing or file I/O.
%
% Trials analysed: All Data trials in which the TARGET WAS SHOWN (TGTCD)
% and the value is known. Breaks before target onset are never counted.
% Value comes from the GOODOBJ/BADOBJ code. That code is sent at a
% task-dependent time (object task: usually before target onset; OneDR:
% 200-400 ms after), so a fixation break between target onset and the code
% leaves the trial uncoded. Its value is then recovered (InferMissingValue):
%   location task : target side (LEFTHEMI/RIGHTHEMI, or target-ID code
%                   mapped to side) vs the current high-value side, taken
%                   from the neighbouring coded trials; left unknown if a
%                   reversal falls between them
%   object task   : target-ID code (1000-1999) -> value learned from coded
%                   trials
% A rule is used only if it reproduces the coded trials' values with
% accuracy >= MinInferenceAccuracy. Unresolved trials: NumTargetNoValue.
%
% Break-fixation rate: fixation broken after target onset and before the
% go signal ('DelayBreak'). Error rate: not completed (no CORRECTCD).
% Error stage:
%   'DelayBreak' : BROKEFIX logged before any FIXOFF (broke fixation
%                  before the go signal; a FIXOFF logged after the break
%                  is task cleanup, as in the object task)
%   'NoSaccade'  : FAILED_SAC, or FIXOFF without SAC_OCCURRED
%   'GazeBreak'  : any other error after go (FAILED_GAZ: left target early)
%
% Anticipation: among completed trials, saccade onset (DetectTargetSaccade,
% via ComputeValueSaccadeBehavior) earlier than AnticipationMs after the
% go signal. Also reports onset before the go signal (RT < 0).
%
% Reversal exclusion (OneDR): the first NumSwitchExclude trials after each
% value reversal are excluded from both rates. Reversals are taken from
% the Good Data trials (same convention as ComputeValueSaccadeBehavior)
% and placed right after the last completed pre-reversal trial in All Data.
%
% Inputs:
%   GoodData  : DataHub GoodData struct (completed trials)
%   AllData   : DataHub AllData struct (all trials) of the same session
%   GoodRange : [StartTrial EndTrial] in Good Data trial numbers
%   Params    : struct, fields override DefaultParams() (and are passed to
%               ComputeValueSaccadeBehavior)
%
% Output Perf struct:
%   .Params
%   .TrialTable   one row per All Data trial in range with target shown
%                 and known value: AllTrialIdx, IsGood, ValueInferred,
%                 Completed, ErrorStage, TrialsSinceSwitch, Included
%   .NumTargetNoValue  target shown but value unknown (not in rates)
%   .NumValueInferred  trials whose value was recovered
%   .ValueInference    struct: Method, Accuracy (on coded trials), NInferred
%   .Error        struct: High/Low .N .NErr .Rate, .ByStage (table), .P (Fisher)
%   .BreakFix     struct: High/Low .N .NBreak .Rate, .P (Fisher); fixation
%                 broken before the go signal (ErrorStage 'DelayBreak')
%   .Anticipation struct: High/Low .N .NAnt .Rate .NPreGo .RatePreGo, .P, .PPreGo
%   .Behavior     output of ComputeValueSaccadeBehavior on the Good trials
%                 (run without reversal exclusion)
%
% Xuefei Yu

ECode;
if nargin < 4 || isempty(Params)
    Params = struct();
end
Def = DefaultParams();
f = fieldnames(Params);
for i = 1:numel(f)
    Def.(f{i}) = Params.(f{i});
end
Params = Def;

%% ---- Map the Good-trial range onto All Data ----
EC        = AllData.EventChannel;
Completed = any(EC == CORRECTCD, 2);
CorrectIdx = find(Completed);                          % All index of each Good trial
GoodRange = [max(1, GoodRange(1)), min(numel(CorrectIdx), GoodRange(2))];
AllFirst = 1;
if GoodRange(1) > 1
    AllFirst = CorrectIdx(GoodRange(1) - 1) + 1;       % errors before first Good trial
end
AllLast = AllData.TrialNum;
if GoodRange(2) < numel(CorrectIdx)
    AllLast = CorrectIdx(GoodRange(2));
end
AllIdx = (AllFirst:AllLast)';

%% ---- Value label and error stage (All Data) ----
HasGood   = any(EC(AllIdx,:) == GOODOBJ, 2);
HasBad    = any(EC(AllIdx,:) == BADOBJ, 2);
HasTarget = any(EC(AllIdx,:) == TGTCD, 2);
Known     = HasGood | HasBad;
HasFixOff = any(EC(AllIdx,:) == FIXOFF, 2);
HasSac    = any(EC(AllIdx,:) == SAC_OCCURRED, 2);
Comp      = Completed(AllIdx);

% Error stage from the task's outcome codes and their timing. Event
% presence alone is not enough: in the object task the fixation point is
% switched off (FIXOFF) as cleanup AFTER a fixation break, so "no FIXOFF"
% does not identify breaks there.
ET = AllData.EventTimeChannel;
Stage = repmat({'Completed'}, numel(AllIdx), 1);
for r = 1:numel(AllIdx)
    if Comp(r)
        continue
    end
    Ev = EC(AllIdx(r),:); Tm = ET(AllIdx(r),:);
    tBreak  = Tm(find(Ev == BROKEFIX, 1));
    tFixOff = Tm(find(Ev == FIXOFF, 1));
    if ~isempty(tBreak) && (isempty(tFixOff) || tFixOff >= tBreak)
        Stage{r} = 'DelayBreak';                   % broke fixation before go
    elseif any(Ev == FAILED_SAC) || (HasFixOff(r) && ~HasSac(r))
        Stage{r} = 'NoSaccade';                    % no saccade after go
    else
        Stage{r} = 'GazeBreak';                    % FAILED_GAZ / broke after go
    end
end

%% ---- Value for target-shown trials without a value code ----
IsLocation = ismember(AllData.TaskCode, [ONEDR_TASK, OBJ_LOCATIONVALUE_TASK]);
[InfGood, Inferred, Perf.ValueInference] = InferMissingValue(EC(AllIdx,:), ...
    HasGood, Known, HasTarget, IsLocation, Params.MinInferenceAccuracy);
HasGood = (HasGood | (Inferred & InfGood))  & HasTarget;
HasBad  = (HasBad  | (Inferred & ~InfGood)) & HasTarget;
Known   = HasGood | HasBad;                   % target shown, value known

%% ---- Saccade behavior on Good trials (anticipation + reversals) ----
% Run without reversal exclusion so reversal points are visible; the
% exclusion is applied below to both rates.
BehParams = Params;
BehParams.NumSwitchExclude = 0;
Behavior = ComputeValueSaccadeBehavior(GoodData, GoodRange(1), GoodRange(2), BehParams);

% Trials since reversal in All Data (location tasks; Inf for the object
% task, see ComputeValueSaccadeBehavior): reversal placed right after the
% last completed pre-reversal trial.
SinceSwitch = inf(numel(AllIdx), 1);
GoodSince   = Behavior.TrialsSinceSwitch;
FirstAfter  = Behavior.TrialIdx(GoodSince == 0);       % Good trial numbers
for s = 1:numel(FirstAfter)
    g = FirstAfter(s);
    if g < 2
        continue
    end
    SwitchAll = CorrectIdx(g - 1) + 1;                 % first trial after last pre-reversal correct
    Rows = find(AllIdx >= SwitchAll);
    SinceSwitch(Rows) = (0:numel(Rows)-1)';
end
Included = Known & SinceSwitch >= Params.NumSwitchExclude;

Perf.Params = Params;
Perf.TrialTable = table(AllIdx(Known), HasGood(Known), Inferred(Known), Comp(Known), ...
    Stage(Known), SinceSwitch(Known), Included(Known), 'VariableNames', ...
    {'AllTrialIdx', 'IsGood', 'ValueInferred', 'Completed', 'ErrorStage', ...
    'TrialsSinceSwitch', 'Included'});
Perf.NumTargetNoValue = sum(HasTarget & ~Known);
Perf.NumValueInferred = sum(Inferred & Known);

%% ---- Error rate ----
Stages = {'DelayBreak', 'NoSaccade', 'GazeBreak'};
Labels = {'High', 'Low'};
Masks  = {Included & HasGood, Included & HasBad};
ByStage = zeros(2, numel(Stages));
for v = 1:2
    E.N    = sum(Masks{v});
    E.NErr = sum(Masks{v} & ~Comp);
    E.Rate = E.NErr / E.N;
    Perf.Error.(Labels{v}) = E;
    for s = 1:numel(Stages)
        ByStage(v, s) = sum(Masks{v} & strcmp(Stage, Stages{s}));
    end
end
Perf.Error.ByStage = array2table(ByStage, 'VariableNames', Stages, 'RowNames', Labels);
Perf.Error.P = FisherP(Perf.Error.High.NErr, Perf.Error.High.N, Perf.Error.Low.NErr, Perf.Error.Low.N);

%% ---- Break-fixation rate (fixation broken before the go signal) ----
for v = 1:2
    F.N      = sum(Masks{v});
    F.NBreak = sum(Masks{v} & strcmp(Stage, 'DelayBreak'));
    F.Rate   = F.NBreak / F.N;
    Perf.BreakFix.(Labels{v}) = F;
end
Perf.BreakFix.P = FisherP(Perf.BreakFix.High.NBreak, Perf.BreakFix.High.N, ...
    Perf.BreakFix.Low.NBreak, Perf.BreakFix.Low.N);

%% ---- Anticipation (completed trials) ----
RT = Behavior.RT;
Keep = ~isnan(RT) & GoodSince >= Params.NumSwitchExclude;
Gm = {Keep & Behavior.IsGood, Keep & ~Behavior.IsGood};
for v = 1:2
    A.N         = sum(Gm{v});
    A.NAnt      = sum(Gm{v} & RT < Params.AnticipationMs);
    A.Rate      = A.NAnt / A.N;
    A.NPreGo    = sum(Gm{v} & RT < 0);
    A.RatePreGo = A.NPreGo / A.N;
    Perf.Anticipation.(Labels{v}) = A;
end
H = Perf.Anticipation.High; L = Perf.Anticipation.Low;
Perf.Anticipation.P      = FisherP(H.NAnt, H.N, L.NAnt, L.N);
Perf.Anticipation.PPreGo = FisherP(H.NPreGo, H.N, L.NPreGo, L.N);
Perf.Anticipation.ThresholdMs = Params.AnticipationMs;

Perf.Behavior = Behavior;
end

%% ======================= Local functions =======================
function P = DefaultParams()
P.AnticipationMs   = 80;  % saccade onset earlier than this re: go = anticipatory
P.NumSwitchExclude = 0;   % 5 for OneDR
P.MinInferenceAccuracy = 0.95; % value-recovery rule must reproduce coded trials this well
end

function [InfGood, Inferred, Info] = InferMissingValue(EC, HasGood, Known, HasTarget, IsLocation, MinAcc)
% Recover the value of target-shown trials without a GOODOBJ/BADOBJ code.
% Inputs: EC (event codes, trials x codes), HasGood/Known/HasTarget
% (logical per trial), IsLocation (task), MinAcc (validation threshold).
% Outputs: InfGood (recovered value), Inferred (recovered trials), Info.
ECode;
n = size(EC, 1);
InfGood = false(n, 1); Inferred = false(n, 1);
Info = struct('Method', '', 'Accuracy', NaN, 'NInferred', 0);

Side = nan(n, 1);
Side(any(EC == LEFTHEMI, 2))  = 1;
Side(any(EC == RIGHTHEMI, 2)) = 2;
ID = nan(n, 1);
for i = 1:n
    c = EC(i, EC(i,:) >= 1000 & EC(i,:) < 2000);
    if ~isempty(c)
        ID(i) = c(1);
    end
end

% Fill a missing side code from the target-ID -> side map, if consistent
Both = ~isnan(Side) & ~isnan(ID);
UID = unique(ID(Both));
SideOfID = nan(size(UID)); Consistent = ~isempty(UID);
for u = 1:numel(UID)
    S = Side(Both & ID == UID(u));
    SideOfID(u) = mode(S);
    Consistent = Consistent && mean(S == SideOfID(u)) >= MinAcc;
end
if Consistent
    Fill = isnan(Side) & ismember(ID, UID);
    [~, Loc] = ismember(ID(Fill), UID);
    Side(Fill) = SideOfID(Loc);
end

Target = HasTarget & ~Known;
if ~any(Target)
    return
end

if IsLocation
    Info.Method = 'side+block';
    KS = find(Known & ~isnan(Side));                 % coded trials with side
    GS = Side(KS);
    GS(~HasGood(KS)) = 3 - GS(~HasGood(KS));         % high-value side
    % Validation on coded trials whose two neighbours agree
    Ok = 0; Nv = 0;
    for j = 2:numel(KS) - 1
        if GS(j-1) == GS(j+1)
            Nv = Nv + 1;
            Ok = Ok + ((Side(KS(j)) == GS(j-1)) == HasGood(KS(j)));
        end
    end
    Info.Accuracy = Ok / Nv;
    if Nv == 0 || Info.Accuracy < MinAcc
        return
    end
    for i = find(Target & ~isnan(Side))'
        Cand = [GS(find(KS < i, 1, 'last')); GS(find(KS > i, 1, 'first'))];
        if isempty(Cand) || numel(unique(Cand)) > 1
            continue                                   % reversal in between: unknown
        end
        Inferred(i) = true;
        InfGood(i)  = Side(i) == Cand(1);
    end
else
    Info.Method = 'objectID';
    KI = find(Known & ~isnan(ID));
    Ok = 0; Nv = 0;
    for j = KI'
        Same = Known & ID == ID(j); Same(j) = false;   % leave-one-out
        if any(Same)
            Nv = Nv + 1;
            Ok = Ok + ((mean(HasGood(Same)) > 0.5) == HasGood(j));
        end
    end
    Info.Accuracy = Ok / Nv;
    if Nv == 0 || Info.Accuracy < MinAcc
        return
    end
    for i = find(Target & ~isnan(ID))'
        Same = Known & ID == ID(i);
        if any(Same) && (all(HasGood(Same)) || ~any(HasGood(Same)))
            Inferred(i) = true;
            InfGood(i)  = HasGood(find(Same, 1));
        end
    end
end
Info.NInferred = sum(Inferred);
end

function p = FisherP(k1, n1, k2, n2)
% Two-sided Fisher exact test on [k1 n1-k1; k2 n2-k2]
if n1 == 0 || n2 == 0
    p = NaN;
    return
end
[~, p] = fishertest([k1, n1 - k1; k2, n2 - k2]);
end
