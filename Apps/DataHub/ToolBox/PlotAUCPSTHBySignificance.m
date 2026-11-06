function PlotAUCPSTHBySignificance(AUCPSTH_cells, p_AUC, TgtAUC, tAxis, titleStr)
% PlotAUCPSTHBySignificance  Plot average AUC PSTH split by significance & direction.
%
% Inputs:
%   AUCPSTH_cells — cell array {1 x N}, each cell is a 1-by-T AUC time series
%   p_AUC         — 1-by-N vector of p-values for target AUC
%   TgtAUC        — 1-by-N vector of scalar AUC values
%   tAxis         — 1-by-T time vector (ms)
%   titleStr      — (optional) string for figure title
%
% Groups:
%   Red:  p < 0.05 & AUC > 0.5  (prefer good)
%   Blue: p < 0.05 & AUC < 0.5  (prefer bad)
%   Grey: p >= 0.05             (non-significant)

if nargin < 5
    titleStr = 'AUC PSTH by Significance';
end

% Ensure column logical vectors
sel_sig_high = (p_AUC < 0.05 & TgtAUC > 0.5)';
sel_sig_low  = (p_AUC < 0.05 & TgtAUC < 0.5)';
sel_ns       = (p_AUC >= 0.05)';

% Compute mean +/- SEM for each group
[AUC_ns_avg, AUC_ns_sem]     = groupAUCStats(AUCPSTH_cells, sel_ns);
[AUC_high_avg, AUC_high_sem] = groupAUCStats(AUCPSTH_cells, sel_sig_high);
[AUC_low_avg, AUC_low_sem]   = groupAUCStats(AUCPSTH_cells, sel_sig_low);

% Plot
figure('Color','w');
hold on;
h_lines = [];
h_labels = {};

if sum(sel_ns) > 0
    h1 = errorbar(tAxis, AUC_ns_avg, AUC_ns_sem, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
    h_lines(end+1) = h1;
    h_labels{end+1} = sprintf('n.s. (n=%d)', sum(sel_ns));
end
if sum(sel_sig_high) > 0
    h2 = errorbar(tAxis, AUC_high_avg, AUC_high_sem, '-r', 'LineWidth', 1.5);
    h_lines(end+1) = h2;
    h_labels{end+1} = sprintf('AUC>0.5, p<0.05 (n=%d)', sum(sel_sig_high));
end
if sum(sel_sig_low) > 0
    h3 = errorbar(tAxis, AUC_low_avg, AUC_low_sem, '-b', 'LineWidth', 1.5);
    h_lines(end+1) = h3;
    h_labels{end+1} = sprintf('AUC<0.5, p<0.05 (n=%d)', sum(sel_sig_low));
end

yline(0.5,'LineStyle', '--','Color','k', 'LineWidth', 0.8);
xline(0, 'LineStyle', '--','Color','k', 'LineWidth', 0.8);

% Significant time bins: ranksum test of AUC vs 0.5 at each time bin per group
yLims = get(gca, 'YLim');
sigY_high = yLims(2) + 0.01 * diff(yLims);
sigY_low  = yLims(2) + 0.03 * diff(yLims);
sigY_ns   = yLims(2) + 0.05 * diff(yLims);

barH = 0.005 * diff(yLims);  % bar thickness

if sum(sel_sig_high) > 1
    mat_high = groupAUCMat(AUCPSTH_cells, sel_sig_high);
    sig_high = timeBinRanksum(mat_high, 0.5);
    plotSigBar(tAxis, sig_high, sigY_high, barH, [1 0 0]);
end
if sum(sel_sig_low) > 1
    mat_low = groupAUCMat(AUCPSTH_cells, sel_sig_low);
    sig_low = timeBinRanksum(mat_low, 0.5);
    plotSigBar(tAxis, sig_low, sigY_low, barH, [0 0 1]);
end
if sum(sel_ns) > 1
    mat_ns = groupAUCMat(AUCPSTH_cells, sel_ns);
    sig_ns_bins = timeBinRanksum(mat_ns, 0.5);
    plotSigBar(tAxis, sig_ns_bins, sigY_ns, barH, [0.6 0.6 0.6]);
end

xlabel('Time from Target Onset (ms)');
ylabel('AUC');
title(sprintf('%s  [red=%d, blue=%d, grey=%d]', ...
    titleStr, sum(sel_sig_high), sum(sel_sig_low), sum(sel_ns)));
if ~isempty(h_lines)
    legend(h_lines, h_labels, 'Location', 'best');
    legend('boxoff');
end
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'off');

end


function mat = groupAUCMat(AUCPSTH_cells, sel)
% Extract AUC PSTH matrix [N x T] for selected neurons
    mat = cell2mat(AUCPSTH_cells(sel)');
    mat = (mat')';
end


function sigBins = timeBinRanksum(mat, chance)
% Ranksum test at each time bin comparing AUC values against chance (0.5)
% Returns logical vector [1 x T]
    nT = size(mat, 2);
    sigBins = false(1, nT);
    for t = 1:nT
        vals = mat(:, t);
        vals = vals(~isnan(vals));
        if length(vals) >= 3
            p = signrank(vals - chance);
            sigBins(t) = p < 0.05;
        end
    end
end


function [avg, sem] = groupAUCStats(AUCPSTH_cells, sel)
% Compute mean and SEM of AUC PSTH for selected neurons
    if sum(sel) == 0
        avg = [];
        sem = [];
        return;
    end
    mat = cell2mat(AUCPSTH_cells(sel)');
    mat = (mat')';
    avg = nanmean(mat, 1);
    sem = nanstd(mat, 0, 1) / sqrt(sum(sel));
end


function plotSigBar(tAxis, sigBins, yPos, barH, clr)
% Plot significance as horizontal filled bars for contiguous significant bins
%   tAxis   — 1-by-T time vector
%   sigBins — 1-by-T logical, true where p < 0.05
%   yPos    — y position of bar center
%   barH    — bar half-height
%   clr     — [r g b] color
    if ~any(sigBins), return; end
    dt = median(diff(tAxis));
    sigIdx = find(sigBins);
    breaks = find(diff(sigIdx) > 1);
    starts = [sigIdx(1), sigIdx(breaks + 1)];
    ends   = [sigIdx(breaks), sigIdx(end)];
    for k = 1:length(starts)
        x1 = tAxis(starts(k)) - dt/2;
        x2 = tAxis(ends(k))   + dt/2;
        fill([x1 x2 x2 x1], [yPos-barH yPos-barH yPos+barH yPos+barH], ...
            clr, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
