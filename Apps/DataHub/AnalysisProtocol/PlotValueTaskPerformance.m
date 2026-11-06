function Fig = PlotValueTaskPerformance(Perf, FigName)
%PLOTVALUETASKPERFORMANCE Plot break-fixation rate (fixation broken before
%the go signal) and anticipatory saccade rate for high vs low value in one
%session.
%Rendering only; takes the output of ComputeValueTaskPerformance.
%
% Inputs:
%   Perf    : struct from ComputeValueTaskPerformance
%   FigName : figure name string
% Output:
%   Fig     : figure handle
%
% Xuefei Yu

% Lab colour convention: object red/blue, location purple/green
if isfield(Perf, 'Behavior') && isfield(Perf.Behavior, 'TaskCode')
    [ColorGood, ColorBad] = ValueColors(Perf.Behavior.TaskCode);
else
    [ColorGood, ColorBad] = ValueColors('Object');
end

Fig = PrepareFigure(112, 'w', [80, 120, 900, 420], 'Name', FigName);

%% Break-fixation rate (fixation broken before the go signal)
ax = subplot(1, 2, 1); hold(ax, 'on');
F = Perf.BreakFix;
R = [F.High.Rate, F.Low.Rate] * 100;
b = bar(ax, [1 2], R, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = [ColorGood; ColorBad];
Cnt = {[F.High.NBreak, F.High.N], [F.Low.NBreak, F.Low.N]};
for v = 1:2
    text(ax, v, R(v), sprintf('%.1f%%\n(%d/%d)', R(v), Cnt{v}(1), Cnt{v}(2)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
end
set(ax, 'XTick', [1 2], 'XTickLabel', {'High', 'Low'}, 'XLim', [0.4 2.6], 'TickDir', 'out');
ylim(ax, [0, max([R, 1]) * 1.35]);
ylabel(ax, 'Break-fixation rate (% of trials)');
title(ax, sprintf('Break fixation before go, p=%.2g (Fisher)', F.P));
box(ax, 'off');

%% Anticipation rate
ax = subplot(1, 2, 2); hold(ax, 'on');
A = Perf.Anticipation;
R = [A.High.Rate, A.Low.Rate] * 100;
b = bar(ax, [1 2], R, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = [ColorGood; ColorBad];
Pre = [A.High.RatePreGo, A.Low.RatePreGo] * 100;
plot(ax, [1 2], Pre, 'kd', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
Cnt = {[A.High.NAnt, A.High.N], [A.Low.NAnt, A.Low.N]};
for v = 1:2
    text(ax, v, R(v), sprintf('%.0f%%\n(%d/%d)', R(v), Cnt{v}(1), Cnt{v}(2)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
end
set(ax, 'XTick', [1 2], 'XTickLabel', {'High', 'Low'}, 'XLim', [0.4 2.6], 'TickDir', 'out');
ylim(ax, [0, max([R, 1]) * 1.35]);
ylabel(ax, '% of completed trials');
legend(ax, {sprintf('RT < %d ms', A.ThresholdMs), 'RT < 0 (before go)'}, 'Box', 'off', 'Location', 'northwest');
title(ax, sprintf('Anticipatory saccades, p=%.2g (Fisher)', A.P));
box(ax, 'off');

Note = sprintf('Trials with target shown; value recovered for %d uncoded trials', Perf.NumValueInferred);
if Perf.NumTargetNoValue > 0
    Note = sprintf('%s; %d with unknown value (not counted)', Note, Perf.NumTargetNoValue);
end
annotation(Fig, 'textbox', [0.05 0.00 0.94 0.05], 'String', Note, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'FontSize', 9);
end
