function fig = PlotGroupAUCTimecourse(TCObj, tObj, TCLoc, tLoc, groupNames, tLim, titleStr, figNum)
% PlotGroupAUCTimecourse  One panel per group: sign-aligned mean +/- SEM AUC
% for object (red) and location (purple) tasks, with per-bin significance bars.
%
% Inputs
%   TCObj/TCLoc — struct arrays from ComputeGroupAUCTimecourse (same groups)
%   tObj/tLoc   — 1-by-T time vectors (ms)
%   groupNames  — cellstr{G}
%   tLim        — [tmin tmax] (ms); titleStr; figNum
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-05

colObj = [0.85 0.15 0.15];
colLoc = [0.55 0.10 0.80];
G = numel(groupNames);
fig = PrepareFigure(figNum, 'w', [50 100 350*G 350], 'Name', titleStr);
yL = [0.4 1];
for g = 1:G
    subplot(1, G, g); hold on;
    DrawTrace(tObj, TCObj(g), colObj);
    DrawTrace(tLoc, TCLoc(g), colLoc);
    DrawSigBars(tObj, TCObj(g).sig, yL(2) - 0.03, colObj);
    DrawSigBars(tLoc, TCLoc(g).sig, yL(2) - 0.06, colLoc);
    yline(0.5, '--k'); xline(0, '--k');
    xlim(tLim); ylim(yL);
    xlabel('Time from target onset (ms)');
    if g == 1, ylabel('AUC (sign-aligned)'); end
    title(sprintf('%s\nn=%d', groupNames{g}, TCObj(g).n));
    box off;
end
legend({'Object (stable)', 'Location (oneDR)'}, 'Location', 'southeast', 'Box', 'off');
sgtitle(titleStr, 'Interpreter', 'none');
end

function DrawTrace(t, tc, col)
if tc.n == 0, return; end
fill([t, fliplr(t)], [tc.mean + tc.sem, fliplr(tc.mean - tc.sem)], col, ...
     'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t, tc.mean, '-', 'Color', col, 'LineWidth', 1.5);
end

function DrawSigBars(t, sig, y, col)
if ~any(sig), return; end
dt = median(diff(t));
d  = diff([false, sig(:)', false]);
s  = find(d == 1); e = find(d == -1) - 1;
for k = 1:numel(s)
    plot([t(s(k)) - dt/2, t(e(k)) + dt/2], [y y], '-', 'Color', col, ...
         'LineWidth', 3, 'HandleVisibility', 'off');
end
end
