% PlotDecoderResults_ObjLoc_FEF.m
% Load saved decoder outputs from Batch_ObjLoc_PopPSTH_FEF and produce six figures:
%   Figure 1 — Within-context accuracy (Stable | OneDR)
%   Figure 2 — Cross-context accuracy (Train Stable→Test OneDR | Train OneDR→Test Stable)
%   Figure 3 — Axis angle curve (same time bin, Stable vs OneDR weight vectors)
%   Figure 4 — Axis angle 2D matrix (Stable t_train × OneDR t_test)
%   Figure 5 — Axis angle at max-beta time point
%   Figure 6 — Axis angle of average beta 100–200 ms
%
% Inputs (loaded from DataPath):
%   Decoder_Stable_RF_FEF_Others.mat   — acc_all, weights_all, acc_perm, time_out
%   Decoder_OneDR_RF_FEF_Others.mat    — acc_all_oneDR, weights_all_oneDR, acc_perm_oneDR, time_out_oneDR
%   Decoder_Cross_Stable2OneDR_RF_FEF.mat — acc_s2d, acc_perm_s2d, time_out
%   Decoder_Cross_OneDR2Stable_RF_FEF.mat — acc_d2s, acc_perm_d2s, time_out

clear all;

%% ── Paths ─────────────────────────────────────────────────────────────────────
DataPath = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/ObjectLocationPSTH/FEF';

%% ── Step 1: Load decoder files ────────────────────────────────────────────────
fprintf('Loading Stable decoder...\n');
S = load(fullfile(DataPath, 'Decoder_Stable_RF_FEF_Others.mat'), ...
    'acc_all', 'weights_all', 'acc_perm', 'time_out');

fprintf('Loading OneDR decoder...\n');
D = load(fullfile(DataPath, 'Decoder_OneDR_RF_FEF_Others.mat'), ...
    'acc_all_oneDR', 'weights_all_oneDR', 'acc_perm_oneDR', 'time_out_oneDR');

fprintf('Loading cross S->D decoder...\n');
S2D = load(fullfile(DataPath, 'Decoder_Cross_Stable2OneDR_RF_FEF.mat'));

fprintf('Loading cross D->S decoder...\n');
D2S = load(fullfile(DataPath, 'Decoder_Cross_OneDR2Stable_RF_FEF.mat'));

% Unpack for convenience
acc_stable    = S.acc_all;          % [n_rep × T]
w_stable_raw  = S.weights_all;      % [n_rep × T × N]
perm_stable   = S.acc_perm;         % [n_perm × T]
time_out      = S.time_out;         % [1 × T]

acc_oneDR     = D.acc_all_oneDR;
w_oneDR_raw   = D.weights_all_oneDR;
perm_oneDR    = D.acc_perm_oneDR;
time_out_d    = D.time_out_oneDR;

acc_s2d       = S2D.acc_s2d;
acc_d2s       = D2S.acc_d2s;

T        = numel(time_out);
n_rep    = size(acc_stable, 1);
n_perm   = size(perm_stable, 1);

% Cross-context permutation baselines — use stable perm as fallback if absent
if isfield(S2D, 'acc_perm_s2d')
    perm_s2d = S2D.acc_perm_s2d;
else
    warning('acc_perm_s2d not found in file — using stable permutation as baseline.');
    perm_s2d = perm_stable;
end

if isfield(D2S, 'acc_perm_d2s')
    perm_d2s = D2S.acc_perm_d2s;
else
    warning('acc_perm_d2s not found in file — using stable permutation as baseline.');
    perm_d2s = perm_stable;
end

%% ── Step 2: Compute significance (cluster-based permutation) ──────────────────
fprintf('Computing cluster-based significance...\n');

sig_stable = ComputeClusterSig(acc_stable, perm_stable, 0.05, n_rep, n_perm, T);
sig_oneDR  = ComputeClusterSig(acc_oneDR,  perm_oneDR,  0.05, n_rep, n_perm, T);
sig_s2d    = ComputeClusterSig(acc_s2d,    perm_s2d,    0.05, n_rep, n_perm, T);
sig_d2s    = ComputeClusterSig(acc_d2s,    perm_d2s,    0.05, n_rep, n_perm, T);

fprintf('Sig bins — Stable:%d  OneDR:%d  S->D:%d  D->S:%d\n', ...
    sum(sig_stable), sum(sig_oneDR), sum(sig_s2d), sum(sig_d2s));

%% ── Color palette ─────────────────────────────────────────────────────────────
col_stable = [0.20 0.50 0.90];
col_oneDR  = [0.90 0.40 0.10];
col_s2d    = [0.10 0.60 0.30];
col_d2s    = [0.70 0.30 0.10];
col_perm   = [0.60 0.60 0.60];

%% ════════════════════════════════════════════════════════════════════════════
%  Figure 1 — Within-context decoding accuracy
%% ════════════════════════════════════════════════════════════════════════════
figure('Name', 'Within-context accuracy', 'NumberTitle', 'off', ...
    'Position', [100 400 900 380]);

% --- Stable ---
subplot(1, 2, 1); hold on;
PlotAccuracy(time_out, acc_stable, perm_stable, sig_stable, col_stable, col_perm);
title('Value decoder — Stable task, RF, FEF');

% --- OneDR ---
subplot(1, 2, 2); hold on;
PlotAccuracy(time_out_d, acc_oneDR, perm_oneDR, sig_oneDR, col_oneDR, col_perm);
title('Value decoder — OneDR task, RF, FEF');

%% ════════════════════════════════════════════════════════════════════════════
%  Figure 2 — Cross-context decoding accuracy
%% ════════════════════════════════════════════════════════════════════════════
figure('Name', 'Cross-context accuracy', 'NumberTitle', 'off', ...
    'Position', [100 400 900 380]);

% --- Train Stable → Test OneDR ---
subplot(1, 2, 1); hold on;
PlotAccuracy(time_out, acc_s2d, perm_s2d, sig_s2d, col_s2d, col_perm);
title('Cross-context: Train Stable \rightarrow Test OneDR, RF, FEF');

% --- Train OneDR → Test Stable ---
subplot(1, 2, 2); hold on;
PlotAccuracy(time_out, acc_d2s, perm_d2s, sig_d2s, col_d2s, col_perm);
title('Cross-context: Train OneDR \rightarrow Test Stable, RF, FEF');

%% ── Step 5: Mean weight matrices ──────────────────────────────────────────────
% Average over repeats: [n_rep × T × N] → [T × N]
w_stable = squeeze(mean(w_stable_raw, 1));   % [T × N]
w_oneDR  = squeeze(mean(w_oneDR_raw,  1));   % [T × N]

%% ════════════════════════════════════════════════════════════════════════════
%  Figure 3 — Same-time axis angle curve
%% ════════════════════════════════════════════════════════════════════════════
% angle(t) = acosd( dot(w_s(t,:), w_d(t,:)) / (||w_s(t,:)|| * ||w_d(t,:)||) )
angle_curve = zeros(1, T);
for t = 1:T
    ws = w_stable(t, :);
    wd = w_oneDR(t, :);
    angle_curve(t) = acosd(dot(ws, wd) / (norm(ws) * norm(wd) + 1e-12));
end

figure('Name', 'Axis angle — same time', 'NumberTitle', 'off', ...
    'Position', [100 400 520 340]);
hold on;
plot(time_out, angle_curve, 'k-', 'LineWidth', 1.8);
yline(90, 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
text(time_out(end), 90, '  90° (orthogonal)', 'Color', [0.5 0.5 0.5], ...
    'VerticalAlignment', 'bottom', 'FontSize', 8);
xline(0, 'LineStyle', '--', 'Color', 'k');
xlabel('Time (ms)');
ylabel('Axis angle (°)');
title('Coding-axis angle: Stable vs OneDR (same time bin)');
ylim([0 180]);
xlim([time_out(1) time_out(end)]);
box off;

%% ════════════════════════════════════════════════════════════════════════════
%  Figure 4 — 2D angle matrix (Stable t_train × OneDR t_test)
%% ════════════════════════════════════════════════════════════════════════════
% Pre-normalise rows for speed
ws_norm = w_stable ./ (vecnorm(w_stable, 2, 2) + 1e-12);   % [T × N]
wd_norm = w_oneDR  ./ (vecnorm(w_oneDR,  2, 2) + 1e-12);   % [T × N]

cos_mat     = ws_norm * wd_norm';          % [T × T], dot products
cos_mat     = max(min(cos_mat, 1), -1);    % clamp to [-1,1]
angle_mat   = acosd(cos_mat);              % [T(stable) × T(oneDR)]

figure('Name', 'Axis angle — 2D matrix', 'NumberTitle', 'off', ...
    'Position', [100 400 480 400]);
imagesc(time_out, time_out, angle_mat);
axis xy;
hold on;
xline(0, 'LineStyle', '--', 'Color', 'w', 'LineWidth', 1);
yline(0, 'LineStyle', '--', 'Color', 'w', 'LineWidth', 1);
cb = colorbar;
cb.Label.String = 'Axis angle (°)';
colormap('jet');
clim([0 180]);
xlabel('OneDR time (ms)');
ylabel('Stable time (ms)');
title('Coding-axis angle: Stable t_{train} \times OneDR t_{test}');
axis square;

%% ════════════════════════════════════════════════════════════════════════════
%  Figure 5 — Axis angle at max-beta time point
%% ════════════════════════════════════════════════════════════════════════════
% Max-beta: time bin where L2 norm of mean weight vector is largest
[~, t_max_s] = max(vecnorm(w_stable, 2, 2));
[~, t_max_d] = max(vecnorm(w_oneDR,  2, 2));

ws_max = w_stable(t_max_s, :);
wd_max = w_oneDR(t_max_d, :);
angle_maxbeta = acosd(dot(ws_max, wd_max) / (norm(ws_max) * norm(wd_max) + 1e-12));

fprintf('Max-beta time: Stable = %.0f ms (bin %d)  |  OneDR = %.0f ms (bin %d)\n', ...
    time_out(t_max_s), t_max_s, time_out(t_max_d), t_max_d);
fprintf('Axis angle at max-beta: %.1f°\n', angle_maxbeta);

figure('Name', 'Axis angle — max beta', 'NumberTitle', 'off', ...
    'Position', [100 400 380 340]);
bar(1, angle_maxbeta, 0.5, 'FaceColor', [0.35 0.35 0.35], 'EdgeColor', 'none');
hold on;
yline(90, 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
ylim([0 180]);
yticks(0:30:180);
xticks(1);
xticklabels({sprintf('Stable t=%.0fms\nvs OneDR t=%.0fms', ...
    time_out(t_max_s), time_out(t_max_d))});
ylabel('Axis angle (°)');
title(sprintf('Coding-axis angle at max-\\beta: %.1f°', angle_maxbeta));
box off;

%% ════════════════════════════════════════════════════════════════════════════
%  Figure 6 — Axis angle of average beta 100–200 ms
%% ════════════════════════════════════════════════════════════════════════════
t_win = time_out >= 100 & time_out <= 200;
if sum(t_win) == 0
    warning('No time bins in 100–200 ms window. Check time_out range.');
    ws_avg = w_stable(1, :);
    wd_avg = w_oneDR(1, :);
else
    ws_avg = mean(w_stable(t_win, :), 1);   % [1 × N]
    wd_avg = mean(w_oneDR(t_win, :),  1);
end

angle_avgbeta = acosd(dot(ws_avg, wd_avg) / (norm(ws_avg) * norm(wd_avg) + 1e-12));

fprintf('Axis angle (avg beta 100–200 ms): %.1f°\n', angle_avgbeta);

figure('Name', 'Axis angle — avg beta 100-200ms', 'NumberTitle', 'off', ...
    'Position', [100 400 380 340]);
bar(1, angle_avgbeta, 0.5, 'FaceColor', [0.35 0.35 0.35], 'EdgeColor', 'none');
hold on;
yline(90, 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
ylim([0 180]);
yticks(0:30:180);
xticks(1);
xticklabels({'Stable vs OneDR (100–200 ms)'});
ylabel('Axis angle (°)');
title(sprintf('Coding-axis angle, avg \\beta 100–200 ms: %.1f°', angle_avgbeta));
box off;

fprintf('\nAll figures generated.\n');

%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
%% ════════════════════════════════════════════════════════════════════════════

function sig_corrected = ComputeClusterSig(acc_all, acc_perm, alpha, n_rep, n_perm, T)
% ComputeClusterSig  Cluster-based permutation test for decoder accuracy.
%
% Inputs:
%   acc_all   [n_rep × T]   — real decoder accuracies
%   acc_perm  [n_perm × T]  — permutation accuracies (shuffled labels)
%   alpha     scalar        — uncorrected threshold (e.g. 0.05)
%   n_rep     scalar
%   n_perm    scalar
%   T         scalar        — number of time bins
% Output:
%   sig_corrected [1 × T] logical — cluster-corrected significant bins

mean_acc  = mean(acc_all,  1);

% Per-bin p-value: fraction of permutation means >= real mean
p_val = zeros(1, T);
for t = 1:T
    p_val(t) = mean(acc_perm(:, t) >= mean_acc(t));
end

% Find contiguous clusters
sig_bins  = p_val < alpha;
clusters  = bwlabel(sig_bins);
n_cl      = max(clusters);

if n_cl == 0
    sig_corrected = false(1, T);
    return;
end

cluster_mass = zeros(1, n_cl);
for c = 1:n_cl
    cluster_mass(c) = sum(1 - p_val(clusters == c));
end

% Null cluster-mass distribution
null_mass = zeros(n_perm, 1);
for perm = 1:n_perm
    p_perm_t = zeros(1, T);
    for t = 1:T
        p_perm_t(t) = mean(acc_perm(:, t) >= acc_perm(perm, t));
    end
    sp = bwlabel(p_perm_t < alpha);
    nc = max(sp);
    if nc > 0
        masses = zeros(1, nc);
        for c = 1:nc
            masses(c) = sum(1 - p_perm_t(sp == c));
        end
        null_mass(perm) = max(masses);
    end
end

thresh = prctile(null_mass, 95);
sig_corrected = false(1, T);
for c = 1:n_cl
    if cluster_mass(c) > thresh
        sig_corrected(clusters == c) = true;
    end
end
end  % ComputeClusterSig


function PlotAccuracy(time_ax, acc_all, acc_perm, sig_corrected, col_line, col_perm)
% PlotAccuracy  Plot mean ± SEM accuracy with permutation baseline and sig bar.
%
% Inputs:
%   time_ax       [1 × T]        — time axis (ms)
%   acc_all       [n_rep × T]    — real accuracies
%   acc_perm      [n_perm × T]   — permutation accuracies
%   sig_corrected [1 × T] logical
%   col_line      [1×3]          — RGB for real accuracy line/fill
%   col_perm      [1×3]          — RGB for permutation line/fill

n_rep  = size(acc_all,  1);
n_perm = size(acc_perm, 1);

mean_acc  = mean(acc_all,  1);
sem_acc   = std(acc_all,  0, 1) / sqrt(n_rep);
mean_perm = mean(acc_perm, 1);
sem_perm  = std(acc_perm, 0, 1) / sqrt(n_perm);

% Permutation shading
fill([time_ax, fliplr(time_ax)], ...
     [mean_perm + sem_perm, fliplr(mean_perm - sem_perm)], ...
     col_perm, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(time_ax, mean_perm, '-', 'Color', col_perm * 0.7, 'LineWidth', 1);

% Real accuracy shading + line
fill([time_ax, fliplr(time_ax)], ...
     [mean_acc + sem_acc, fliplr(mean_acc - sem_acc)], ...
     col_line, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(time_ax, mean_acc, '-', 'Color', col_line, 'LineWidth', 1.8);

% Significance bar
dt = time_ax(2) - time_ax(1);
y_sig = 0.955;
sig_idx = find(sig_corrected);
for ii = 1:numel(sig_idx)
    tc = time_ax(sig_idx(ii));
    fill([tc - dt/2, tc + dt/2, tc + dt/2, tc - dt/2], ...
         [y_sig - 0.01, y_sig - 0.01, y_sig + 0.01, y_sig + 0.01], ...
         col_line, 'EdgeColor', 'none');
end

yline(0.5,'LineStyle', '--','Color','k');
xline(0, 'LineStyle',  '--','Color','k');
xlabel('Time (ms)');
ylabel('Decoding accuracy');
xlim([time_ax(1) time_ax(end)]);
ylim([0.3 1.0]);
box off;
end  % PlotAccuracy
