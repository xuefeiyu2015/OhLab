function fig = PlotAUCHeatmapPair(aucObj, tObj, aucLoc, tLoc, tLim, cLim, titleStr, figNum, hLines, sortLabel)
% PlotAUCHeatmapPair  Side-by-side AUC-vs-time heatmaps (object | location),
% rows already sorted by the caller (same row order in both panels).
%
% Inputs
%   aucObj/aucLoc — N-by-T AUC matrices (pre-sorted, same row order)
%   tObj/tLoc     — 1-by-T time vectors (ms from target onset)
%   tLim          — [tmin tmax] display range (ms)
%   cLim          — color limits, centered on 0.5
%   titleStr      — figure title; figNum — figure number
%   hLines        — (optional) row positions for horizontal separators,
%                   e.g. where the sort key crosses 0.5
%   sortLabel     — (optional) name of the sort key for the y label (default 'object AUC')
% Output: fig — figure handle
%
% Xuefei Yu, 2026-11-05

if nargin < 9,  hLines = []; end
if nargin < 10, sortLabel = 'object AUC'; end
fig = PrepareFigure(figNum, 'w', [100 100 900 600], 'Name', titleStr);
nN  = size(aucObj, 1);

% Diverging blue-white-red map centered on 0.5
nC   = 128;
up   = linspace(0, 1, nC)';
cmap = [[up; ones(nC,1)], [up; flipud(up)], [ones(nC,1); flipud(up)]];

panels = {aucObj, tObj, 'Object value (stable)'; ...
          aucLoc, tLoc, 'Location value (oneDR)'};
for p = 1:2
    ax = subplot(1, 2, p);
    imagesc(panels{p,2}, 1:nN, panels{p,1}, 'AlphaData', ~isnan(panels{p,1}));
    set(ax, 'Color', [0.85 0.85 0.85], 'YDir', 'reverse');
    colormap(ax, cmap); clim(ax, cLim);
    xline(0, '--k', 'LineWidth', 1);
    for h = hLines(:)'
        yline(h, '-k', 'LineWidth', 1);
    end
    xlim(tLim);
    xlabel('Time from target onset (ms)');
    if p == 1, ylabel(sprintf('Neuron (sorted by %s), n=%d', sortLabel, nN)); end
    title(panels{p,3});
    box off;
end
cb = colorbar(ax); cb.Label.String = 'AUC (high vs low value)';
sgtitle(titleStr, 'Interpreter', 'none');
end
