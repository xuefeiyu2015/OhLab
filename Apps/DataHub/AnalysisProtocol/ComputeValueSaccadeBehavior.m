function Behavior = ComputeValueSaccadeBehavior(Data, StartTrial, EndTrial, Params)
%COMPUTEVALUESACCADEBEHAVIOR Trial-wise saccade metrics for high vs low value
%targets in the object-value (stable value) and location-value (OneDR) tasks.
%
% Pure computation: no plotting, printing or file I/O.
%
% Trial flow: FP on -> target on (TGTCD) -> FP off (FIXOFF, go signal)
%   -> saccade to target -> gaze on target -> reward -> target off (TGTOFF)
%
% Target saccade: detected by DetectTargetSaccade (ported from JLab
% CalculateRT): Savitzky-Golay-smoothed 2-D speed, baseline-relative speed
% threshold, onset searched back from the SAC_OCCURRED marker (eye left the
% fixation window), offset from the speed trace. Trials are flagged as
% Artifact (all metrics NaN) if the amplitude does not match the target
% eccentricity, eye samples are missing around the saccade (blink), or
% RT > MaxRT; NoiseReject if tracker noise falls right after the go cue or
% inside the saccade.
%
% Metrics (per trial):
%   RT       : saccade onset - FIXOFF (ms); can be < 0 (anticipation)
%   PeakV    : peak 2-D eye velocity during the target saccade (deg/s)
%   Amp      : saccade amplitude (deg)
%   GazeDur  : target-saccade end -> time the eye leaves the target window
%              (ms). Window = circle of radius TargetWindowRadius (deg)
%              around the target centre. The eye has left when it is
%              outside the window (or missing, i.e. blink) for at least
%              MinOutsideMs consecutive ms. NaN if the eye does not enter
%              the window within EntryGraceMs of saccade end (LandedOutside)
%              or never leaves before the end of the analysis window
%              (target off + PostTgtOffBuffer; GazeCensored).
%              Leaving the target before reward onset (RWDCD) is
%              impossible on a rewarded trial -> GazeTooShort, NaN. (The
%              task's 400 ms hold starts when the eye enters the task
%              window, mid-saccade, so gaze from saccade END can be
%              slightly < 400 ms on valid trials; reward time is the
%              window-independent check.)
%              NOTE: target stays on until after reward, and reward is
%              longer for high value, so GazeDur is partly task-imposed.
%
% Inputs:
%   Data       : DataHub Data struct (single session)
%   StartTrial : first trial index
%   EndTrial   : last trial index
%   Params     : (optional) struct overriding any default in DefaultParams()
%                e.g. Params.NumSwitchExclude = 5 for OneDR
%
% Output:
%   Behavior struct with fields
%     .Params         parameters used
%     .TaskCode       task code of the session
%     .TrialIdx       [nTrial x 1] original trial numbers kept
%     .IsGood         [nTrial x 1] logical, high value target
%     .TargetAngle    [nTrial x 1] target polar angle (deg)
%     .TrialsSinceSwitch [nTrial x 1] trials since last value switch (Inf if none)
%     .RT, .PeakV, .Amp, .GazeDur  [nTrial x 1]
%     .GazeCensored   [nTrial x 1] logical
%     .LandedOutside  [nTrial x 1] logical, eye never entered target window
%     .TargetXY       [nTrial x 2] target centre (deg)
%     .Artifact       [nTrial x 1] logical, target saccade rejected
%     .NoiseReject    [nTrial x 1] logical, tracker noise (DetectTargetSaccade)
%     .GazeTooShort   [nTrial x 1] logical, eye left target before reward
%     .OnsetStatus, .OffsetStatus  {nTrial x 1} detection path per trial
%     .TargetEcc      [nTrial x 1] target eccentricity (deg)
%     .Summary        struct: Summary.(Metric) with Good/Bad mean, sem, n,
%                     ranksum p and AUC (Good as positive class)
%     .SummaryBySide  struct array, one entry per unique target angle
%
% Xuefei Yu

ECode;

if nargin < 4 || isempty(Params)
    Params = struct();
end
Params = MergeParams(DefaultParams(), Params);

%% ---------------- Trial selection ----------------
Selection = StartTrial:EndTrial;

% Remove photodiode/abnormal trials
AbnormalTrials = AbnormalTrialCheck(Data, [TGTCD, TGTOFF]);
Selection = Selection(~ismember(Selection, AbnormalTrials));

% Keep single-target trials only (if task type info exists)
if isfield(Data, 'ParaLib') && isa(Data.ParaLib, 'containers.Map') && isKey(Data.ParaLib, 'TaskType')
    TaskTypeLib = Data.ParaLib('TaskType');
    SingleTarget = find(TaskTypeLib(:,2) == 1);
    Selection = Selection(ismember(Selection, SingleTarget));
end
Selection = Selection(:);

EventChannel     = Data.EventChannel(Selection,:);
EventTimeChannel = Data.EventTimeChannel(Selection,:);
EventBin         = Data.EventBin(Selection,:);

% Value label from event codes (works for both object and location value)
IsGood  = false(numel(Selection),1);
HasCode = false(numel(Selection),1);
for i = 1:numel(Selection)
    Code = EventChannel(i, ismember(EventChannel(i,:), [GOODOBJ, BADOBJ]));
    if ~isempty(Code)
        HasCode(i) = true;
        IsGood(i)  = Code(1) == GOODOBJ;
    end
end

% Target angle (used for side split in the location value task)
if ndims(Data.TargetPolarLoc) == 3
    TargetAngle = Data.TargetPolarLoc(Selection,2,1);
else
    TargetAngle = Data.TargetPolarLoc(Selection,2);
end
TargetAngle = TargetAngle(:);
if ndims(Data.TargetPolarLoc) == 3
    TargetEcc = Data.TargetPolarLoc(Selection,2,2);
else
    TargetEcc = Data.TargetPolarLoc(Selection,3);
end
TargetEcc = TargetEcc(:);
% Target centre [x y] (deg). 2-D layout (Robin): [ID x y]. 3-D layout
% (Adams): trials x [ID target1 target2] x [x y].
if ndims(Data.TargetLoc) == 3
    TargetXY = [Data.TargetLoc(Selection,2,1), Data.TargetLoc(Selection,2,2)];
else
    TargetXY = Data.TargetLoc(Selection, 2:3);
end

% Trials since value reversal (location value tasks only). A reversal is
% where the high-value side changes: high-value side = target angle on high
% trials, the opposite angle on low trials. Independent of the TargetIndex
% layout (2-D for Robin, 3-D for Adams; TargetIndex(:,2) is all zeros in the
% 3-D layout, so it cannot be used). Object task: Inf (value belongs to the
% object, not the side).
TrialsSinceSwitch = inf(numel(Selection),1);
if ismember(Data.TaskCode, [ONEDR_TASK, OBJ_LOCATIONVALUE_TASK])
    GoodSide = round(TargetAngle);
    GoodSide(~IsGood) = mod(GoodSide(~IsGood) + 180, 360);
    GoodSide(~HasCode) = NaN;
    SwitchIdx = find(diff(GoodSide) ~= 0 & ~isnan(diff(GoodSide))) + 1;
    for k = 1:numel(SwitchIdx)
        Rows = SwitchIdx(k):numel(Selection);
        TrialsSinceSwitch(Rows) = (0:numel(Rows)-1)';
    end
end
KeepMask = HasCode & (TrialsSinceSwitch >= Params.NumSwitchExclude);

Selection         = Selection(KeepMask);
IsGood            = IsGood(KeepMask);
TargetAngle       = TargetAngle(KeepMask);
TargetEcc         = TargetEcc(KeepMask);
TargetXY          = TargetXY(KeepMask,:);
TrialsSinceSwitch = TrialsSinceSwitch(KeepMask);
EventChannel      = EventChannel(KeepMask,:);
EventTimeChannel  = EventTimeChannel(KeepMask,:);
EventBin          = EventBin(KeepMask,:);
nTrial            = numel(Selection);

%% ---------------- Eye data: target on -> target off + buffer ----------------
EyeBinWidth = Data.EyeBinWidth;
TgtOnBin    = FindOutTime(EventChannel, EventBin, TGTCD);
TgtOffBin   = FindOutTime(EventChannel, EventBin, TGTOFF);
TgtOnTime   = FindOutTime(EventChannel, EventTimeChannel, TGTCD);
FixOffTime  = FindOutTime(EventChannel, EventTimeChannel, FIXOFF);
TgtOnBin = TgtOnBin(:,1); TgtOffBin = TgtOffBin(:,1);
FixOffRel   = FixOffTime(:,1) - TgtOnTime(:,1); % go signal re: target on (ms)

WinEndBin = TgtOffBin + round(Params.PostTgtOffBuffer / EyeBinWidth);
MaxLen    = max(WinEndBin - TgtOnBin + 1);
EyeX = nan(nTrial, MaxLen);
EyeY = nan(nTrial, MaxLen);
ValidLen = zeros(nTrial,1); % samples actually recorded per trial
for i = 1:nTrial
    x = Data.EyeDataRaw{Selection(i),1};
    y = Data.EyeDataRaw{Selection(i),2};
    if isnan(TgtOnBin(i)) || isnan(TgtOffBin(i)) || isempty(x)
        continue
    end
    Idx = TgtOnBin(i):min(WinEndBin(i), numel(x));
    EyeX(i,1:numel(Idx)) = x(Idx);
    EyeY(i,1:numel(Idx)) = y(Idx);
    ValidLen(i) = numel(Idx);
end

%% ---------------- Saccade detection (JLab-style, DetectTargetSaccade) ----------------
% Online marker of the eye leaving the fixation window (anchors onset search)
SacOccTime = FindOutTime(EventChannel, EventTimeChannel, SAC_OCCURRED);
ApproxRT   = SacOccTime(:,1) - FixOffTime(:,1);          % ms re: go cue

% Reward onset re: target onset (ms): the eye cannot leave before it
RwdTime = FindOutTime(EventChannel, EventTimeChannel, RWDCD);
RwdRel  = RwdTime(:,1) - TgtOnTime(:,1);

% Required gaze hold on the target (task parameter; reported only)
GazeHoldMs = Params.GazeHoldMs;
if isfield(Data, 'ParaLib') && isa(Data.ParaLib, 'containers.Map') && ...
        isKey(Data.ParaLib, 'GazeHoldForTarget,GazeHoldForTargetVar')
    Hold = Data.ParaLib('GazeHoldForTarget,GazeHoldForTargetVar');
    GazeHoldMs = min(Hold(:,2));
end

RT            = nan(nTrial,1);
PeakV         = nan(nTrial,1);
Amp           = nan(nTrial,1);
GazeDur       = nan(nTrial,1);
GazeCensored  = false(nTrial,1);
GazeTooShort  = false(nTrial,1);
LandedOutside = false(nTrial,1);
Artifact      = false(nTrial,1);
NoiseReject   = false(nTrial,1);
OnsetStatus   = repmat({''}, nTrial, 1);
OffsetStatus  = repmat({''}, nTrial, 1);
TimeAll = (0:MaxLen-1) * EyeBinWidth;                    % ms re: target onset

for i = 1:nTrial
    if isnan(FixOffRel(i)) || ValidLen(i) == 0
        continue
    end
    % Detection window re: go cue
    WinHi = ApproxRT(i) + Params.PostMarkerMs;
    if isnan(WinHi)
        WinHi = Params.MaxRT + Params.MaxSacDurMs;
    end
    Sel = TimeAll - FixOffRel(i) >= Params.PreGoMs & TimeAll - FixOffRel(i) <= WinHi;
    Sel(ValidLen(i)+1:end) = false;
    det = DetectTargetSaccade(TimeAll(Sel) - FixOffRel(i), EyeX(i,Sel), EyeY(i,Sel), ...
        Params, ApproxRT(i));
    NoiseReject(i)  = det.noiseReject;
    OnsetStatus{i}  = det.onsetStatus;
    OffsetStatus{i} = det.offsetStatus;
    if ~det.detected
        continue
    end

    % Artifact checks: amplitude must match target eccentricity, no missing
    % samples within +/- ArtifactPadMs of the saccade, RT within limit
    b1 = round((det.Onset  + FixOffRel(i)) / EyeBinWidth) + 1;   % index in EyeX
    b2 = round((det.Offset + FixOffRel(i)) / EyeBinWidth) + 1;
    Pad = round(Params.ArtifactPadMs / EyeBinWidth);
    Seg = max(1, b1-Pad):min(ValidLen(i), b2+Pad);
    AmpOK = abs(det.Amplitude - TargetEcc(i)) <= Params.AmpTolerance * TargetEcc(i);
    if ~AmpOK || any(isnan(EyeX(i,Seg))) || any(isnan(EyeY(i,Seg))) || det.Onset > Params.MaxRT
        Artifact(i) = true;
        continue
    end

    RT(i)    = det.Onset;
    Amp(i)   = det.Amplitude;
    PeakV(i) = det.PeakVelocity;

    % Gaze duration (target window rule): saccade end -> first sustained
    % exit from the target window. Missing samples (blink) count as outside.
    Post   = b2:ValidLen(i);
    Dist   = hypot(EyeX(i,Post) - TargetXY(i,1), EyeY(i,Post) - TargetXY(i,2));
    Inside = Dist <= Params.TargetWindowRadius; % NaN -> false (outside)
    Entry  = find(Inside(1:min(end, round(Params.EntryGraceMs/EyeBinWidth)+1)), 1, 'first');
    if isempty(Entry)
        LandedOutside(i) = true;
        continue
    end
    MinOut = max(1, round(Params.MinOutsideMs / EyeBinWidth));
    Outside = ~Inside;
    Outside(1:Entry) = false;
    RunLen = conv(double(Outside), ones(1, MinOut), 'valid');
    Leave  = find(RunLen == MinOut, 1, 'first');
    if isempty(Leave)
        GazeCensored(i) = true; % still on target at end of window
    elseif (b2 - 1 + Leave - 1) * EyeBinWidth < RwdRel(i)
        % Left the target before reward onset: impossible on a rewarded
        % trial -> measurement error (window too small, blink, offset)
        GazeTooShort(i) = true;
    else
        GazeDur(i) = (Leave - 1) * EyeBinWidth; % Post(1) = saccade end
    end
end

%% ---------------- Output ----------------
Behavior.Params            = Params;
Behavior.TaskCode          = Data.TaskCode;
Behavior.TrialIdx          = Selection;
Behavior.IsGood            = IsGood;
Behavior.TargetAngle       = TargetAngle;
Behavior.TargetEcc         = TargetEcc;
Behavior.TrialsSinceSwitch = TrialsSinceSwitch;
Behavior.RT                = RT;
Behavior.PeakV             = PeakV;
Behavior.Amp               = Amp;
Behavior.GazeDur           = GazeDur;
Behavior.GazeCensored      = GazeCensored;
Behavior.LandedOutside     = LandedOutside;
Behavior.TargetXY          = TargetXY;
Behavior.Artifact          = Artifact;
Behavior.NoiseReject       = NoiseReject;
Behavior.GazeTooShort      = GazeTooShort;
Behavior.GazeHoldMs        = GazeHoldMs;
Behavior.OnsetStatus       = OnsetStatus;
Behavior.OffsetStatus      = OffsetStatus;
Behavior.MetricNames       = {'RT','PeakV','GazeDur'};

Behavior.Summary = SummarizeGoodBad(Behavior, true(nTrial,1));

Angles = unique(TargetAngle(~isnan(TargetAngle)));
Behavior.SummaryBySide = struct('TargetAngle', {}, 'Summary', {});
for a = 1:numel(Angles)
    Behavior.SummaryBySide(a).TargetAngle = Angles(a);
    Behavior.SummaryBySide(a).Summary = SummarizeGoodBad(Behavior, TargetAngle == Angles(a));
end

end

%% ======================= Local functions =======================
function P = DefaultParams()
% Saccade detection (DetectTargetSaccade, ported from JLab CalculateRT)
P.PreGoMs            = -300;       % detection window start re: go cue (ms)
P.PostMarkerMs       = 300;        % window end = SAC_OCCURRED + this (ms)
P.SgolayWin          = 7;          % Savitzky-Golay window (samples, ~7 ms)
P.SpeedNoise         = 3000;       % deg/s, faster = tracker noise
P.NoiseWinMs         = [0 200];    % noise here (re: go) rejects the trial
P.BaseWinMs          = [-250 -100];% baseline speed window re: go (ms); ends
                                   % before the earliest anticipatory onsets
P.SpeedThr           = 30;         % deg/s above baseline mean (trigger)
P.NContig            = 5;          % samples a criterion must hold
P.OnsetBackMs        = -100;       % onset may be this early re: go (ms)
P.MaxSacDurMs        = 200;        % offset searched within this of onset (ms)
P.EndThrStep         = 10;         % deg/s, raise end threshold if no offset
P.MaxRT              = 600;        % ms re: go; slower -> excluded
% Artifact rejection
P.AmpTolerance       = 0.4;  % |amp - target ecc| <= AmpTolerance*ecc
P.ArtifactPadMs      = 10;   % no missing samples within +/- this of saccade
% Required gaze hold for reward (reported; read from ParaLib if available)
P.GazeHoldMs         = 400;  % ms
% Gaze duration: target window rule
P.TargetWindowRadius = 6;    % deg, window radius around target centre
                             % (Adams landing error up to ~5 deg at 20 deg ecc)
P.EntryGraceMs       = 100;  % eye must enter window within this of saccade end
P.MinOutsideMs       = 20;   % eye must stay outside this long to count as leaving
P.PostTgtOffBuffer   = 1000; % ms, analysis window past target offset
% Trials excluded after each value switch (0 for object value; 5 for OneDR)
P.NumSwitchExclude   = 0;
end

function P = MergeParams(P, Override)
f = fieldnames(Override);
for i = 1:numel(f)
    P.(f{i}) = Override.(f{i});
end
end

function Summary = SummarizeGoodBad(Behavior, Mask)
% Good vs bad mean/sem/n, rank-sum p, and AUC (Good = positive class)
for m = 1:numel(Behavior.MetricNames)
    Name = Behavior.MetricNames{m};
    Val  = Behavior.(Name);
    G = Val(Mask &  Behavior.IsGood); G = G(~isnan(G));
    B = Val(Mask & ~Behavior.IsGood); B = B(~isnan(B));

    S.GoodMean = mean(G);  S.BadMean = mean(B);
    S.GoodSem  = std(G)/sqrt(numel(G));
    S.BadSem   = std(B)/sqrt(numel(B));
    S.GoodN    = numel(G); S.BadN    = numel(B);
    if numel(G) > 1 && numel(B) > 1
        S.P   = ranksum(G, B);
        S.AUC = rocN(G, B, 100);
    else
        S.P = NaN; S.AUC = NaN;
    end
    Summary.(Name) = S;
end
end
