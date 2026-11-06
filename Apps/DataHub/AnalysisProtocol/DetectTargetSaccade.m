function det = DetectTargetSaccade(t, x, y, cfg, approxRT)
%DETECTTARGETSACCADE Detect the go-cue-triggered saccade in one trial and
%measure onset, offset, amplitude and peak velocity.
%
% Ported from JLab CalculateRT/detectSaccade (Xuefei Yu 2026), adapted to
% OHLab: times in ms, calibrated traces in deg, and no target-acquired
% (Choicetime) marker, so the offset is searched from the speed trace only.
%
% Algorithm (pure, no plotting):
%   1. Smooth x/y with a Savitzky-Golay filter (cfg.SgolayWin samples) and
%      take 2-D speed from the smoothed position.
%   2. Speed > cfg.SpeedNoise is tracker noise (blink/glitch): speed and the
%      two position samples around it are set NaN. Noise within cfg.NoiseWinMs
%      after the go cue rejects the trial.
%   3. Baseline speed mean/SD over cfg.BaseWinMs (re: go cue).
%      Trigger level = baseline mean + cfg.SpeedThr; back-trace level =
%      baseline mean + 3 SD.
%   4. Onset trigger: searching BACK from approxRT (SAC_OCCURRED - go cue,
%      i.e. the eye left the fixation window) to cfg.OnsetBackMs, the latest
%      sample i with speed above the trigger for NContig samples (i:i+N-1)
%      and below it before (i-N:i-1). If strict search fails, the pre-run
%      criterion is relaxed (JLab "release"). Without approxRT, the search
%      runs back from go + cfg.MaxRT.
%   5. Onset = trace back from the trigger while speed > trigger level, then
%      further back while speed > 3 SD level or eye eccentricity > half its
%      post-onset maximum.
%   6. Offset: first sample after onset (within cfg.MaxSacDurMs) where speed
%      stays below the end threshold for NContig samples; if none, raise the
%      threshold by cfg.EndThrStep and retry. Then trace forward while speed
%      > 3 SD level.
%   7. Noise between onset and offset rejects the trial.
%
% Inputs:
%   t        : 1 x n time (ms) relative to the go cue (FIXOFF)
%   x, y     : 1 x n eye position (deg)
%   cfg      : struct with SgolayWin, SpeedNoise, NoiseWinMs, BaseWinMs,
%              SpeedThr, NContig, OnsetBackMs, MaxRT, MaxSacDurMs, EndThrStep
%   approxRT : SAC_OCCURRED - go cue (ms); NaN if the marker is missing
%
% Output det struct:
%   .detected, .noiseReject, .onsetStatus ('regular'|'release'|'fail'),
%   .offsetStatus ('regular'|'release'|'raised'|'fail'),
%   .Onset, .Offset (ms re: go), .PeakVelocity (deg/s), .Amplitude (deg),
%   .StartPoint, .EndPoint ([x y] deg), .Duration (ms)
%
% Xuefei Yu

det = struct('detected', false, 'noiseReject', false, ...
    'onsetStatus', 'fail', 'offsetStatus', 'fail', ...
    'Onset', NaN, 'Offset', NaN, 'PeakVelocity', NaN, 'Amplitude', NaN, ...
    'StartPoint', [NaN NaN], 'EndPoint', [NaN NaN], 'Duration', NaN);

if numel(t) < 2 * cfg.NContig + cfg.SgolayWin || all(isnan(x))
    return
end

%% 1. Smoothed position and speed (speed on mid-sample axis tv)
dt = median(diff(t)) / 1000;                 % s
xs = smoothdata(x, 'sgolay', cfg.SgolayWin, 'omitnan');
ys = smoothdata(y, 'sgolay', cfg.SgolayWin, 'omitnan');
xs(isnan(x)) = NaN; ys(isnan(y)) = NaN;      % keep real gaps as gaps
speed = hypot(diff(xs), diff(ys)) / dt;      % deg/s
tv = (t(1:end-1) + t(2:end)) / 2;

%% 2. Tracker noise
noiseMask = speed > cfg.SpeedNoise | isnan(speed);
speed(noiseMask) = NaN;
if any(noiseMask & tv >= cfg.NoiseWinMs(1) & tv <= cfg.NoiseWinMs(2))
    det.noiseReject = true;
    return
end
dev = hypot(xs, ys);

%% 3. Baseline
inBase = tv >= cfg.BaseWinMs(1) & tv <= cfg.BaseWinMs(2);
baseMean = mean(speed(inBase), 'omitnan');
baseSd   = std(speed(inBase), 'omitnan');
if isnan(baseMean)
    return
end
trigLevel = baseMean + cfg.SpeedThr;
thr3sd    = baseMean + 3 * baseSd;
above     = speed > trigLevel;

%% 4. Onset trigger (search back from approxRT)
hi = approxRT;
if isnan(hi)
    hi = cfg.MaxRT;
end
[iTrig, det.onsetStatus] = FindSaccadeRun(above, tv, cfg.NContig, cfg.OnsetBackMs, hi);
if isnan(iTrig)
    return
end

%% 5. Onset back-trace
iCross = iTrig;
while iCross > 1 && speed(iCross - 1) >= trigLevel
    iCross = iCross - 1;
end
maxDev = max(dev(iCross:min(iCross + 500, numel(dev))));
iBase = iCross;
while iBase > 1 && (speed(iBase - 1) > thr3sd || dev(iBase - 1) > maxDev / 2)
    iBase = iBase - 1;
end
onsetT = tv(iBase);

%% 6. Offset (speed only; no target-acquired marker in OHLab)
iDownMax = find(tv <= onsetT + cfg.MaxSacDurMs, 1, 'last');
if isempty(iDownMax) || iDownMax <= iCross
    iDownMax = min(numel(speed), iCross + 1);
end
endThr = trigLevel;
winMax = max(speed(iCross:iDownMax), [], 'omitnan');
Raised = false;
while true
    [iDown, offOk, det.offsetStatus] = FindSaccadeEnd(speed, iCross, iDownMax, endThr, cfg.NContig);
    if offOk || isnan(winMax) || endThr >= winMax
        break
    end
    endThr = endThr + cfg.EndThrStep;
    Raised = true;
end
if offOk
    iOff = iDown;
    if Raised
        det.offsetStatus = 'raised';
    end
else
    iOff = iCross;
end
while iOff < numel(speed) && speed(iOff + 1) > thr3sd
    iOff = iOff + 1;
end
offsetT = tv(iOff);

%% 7. Noise inside the saccade
if any(noiseMask & tv >= onsetT & tv <= offsetT)
    det.noiseReject = true;
    return
end

%% Metrics (from the smoothed trace)
sp0 = [interp1(t, xs, onsetT),  interp1(t, ys, onsetT)];
sp1 = [interp1(t, xs, offsetT), interp1(t, ys, offsetT)];
det.detected     = true;
det.Onset        = onsetT;
det.Offset       = offsetT;
det.StartPoint   = sp0;
det.EndPoint     = sp1;
% Peak over the full onset->offset interval. (JLab takes it from the
% trigger, which misses the peak when the trigger lands on a post-saccadic
% bump after an anticipatory saccade.)
det.PeakVelocity = max(speed(iBase:iOff));
det.Amplitude    = hypot(sp1(1) - sp0(1), sp1(2) - sp0(2));
det.Duration     = offsetT - onsetT;
end

%% ======================= Local functions =======================
function [iStart, Status] = FindSaccadeRun(above, tv, nContig, lo, hi)
% Search BACK from hi to lo for the latest clean rising edge: above(i:i+N-1)
% all true and enough of above(i-N:i-1) false. Strict first (more than
% floor(N/2) of the pre-run below threshold), then released (more than 0).
iStart = NaN; Status = 'fail';
i0  = find(tv >= lo, 1);
iHi = find(tv <= hi, 1, 'last');
if isempty(i0) || isempty(iHi) || iHi < i0
    return
end
idxLast  = find(above, 1, 'last');
idxFirst = find(above, 1, 'first');
if isempty(idxLast)
    return
end
iTop = min([iHi, numel(above) - nContig + 1, idxLast]);
iBot = max([i0, nContig + 1, idxFirst]);
Turning = [floor(nContig/2), 0];
Labels  = {'regular', 'release'};
for pass = 1:2
    for i = iTop:-1:iBot
        if all(above(i:i + nContig - 1)) && sum(~above(i - nContig:i - 1)) > Turning(pass)
            iStart = i; Status = Labels{pass};
            return
        end
    end
end
end

function [iOff, ok, Status] = FindSaccadeEnd(speed, iCross, iMax, thr, nContig)
% First index in (iCross, iMax] where speed was above thr for the N samples
% before and is below thr for enough of the N samples after. Strict first,
% then released.
ok = false; iOff = NaN; Status = 'fail';
below = speed < thr;
if below(iCross)
    iOff = iCross; ok = true; Status = 'regular';
    return
end
Turning = [ceil(nContig/2), 0];
Labels  = {'regular', 'release'};
for pass = 1:2
    for i = (iCross + 1):min(iMax, numel(below) - nContig + 1)
        if i - nContig >= 1 && all(~below(i - nContig:i - 1)) && sum(below(i:i + nContig - 1)) > Turning(pass)
            iOff = i; ok = true; Status = Labels{pass};
            return
        end
    end
end
end
