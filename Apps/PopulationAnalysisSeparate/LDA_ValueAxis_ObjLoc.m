% LDA_ValueAxis_ObjLoc.m
%
% Purpose : Load per-trial PSTH data (PopPSTH_RF_*.mat) for FEF and SC,
%           apply regularised LDA at the population peak window to find the
%           optimal value discriminant axis in each context, compute the
%           angle between axes, and evaluate within- and cross-context
%           decoding accuracy over time.
%
% Data expected (one file per area):
%   DataDir/{Area}/PopPSTH_RF_{Area}.mat  containing:
%     psth_targeton_raw_stable_RF_high   {1 x N}  each [n_trials x T_s]
%     psth_targeton_raw_stable_RF_low    {1 x N}
%     psth_targeton_raw_oneDR_RF_high    {1 x N}  each [n_trials x T_d]
%     psth_targeton_raw_oneDR_RF_low     {1 x N}
%     PSTH_Time_TgtOn_stable             [1 x T_s]  (ms, t=0 = target on)
%     PSTH_Time_TgtOn_OneDR              [1 x T_d]
%
% Xuefei Yu

%% ── Paths ────────────────────────────────────────────────────────────────────
DataDir = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/ObjectLocationPSTH';
Areas   = {'FEF', 'SC'};

% LDA regularisation (Gamma ∈ [0,1]; 0 = standard LDA, >0 shrinks towards
% diagonal.  Use ~0.5 when N_neurons > n_trials to avoid singular covariance.)
Gamma   = 0.5;
K_fold  = 5;      % cross-validation folds for within-context accuracy

%% ── Loop over areas ──────────────────────────────────────────────────────────
for ai = 1:numel(Areas)
    area = Areas{ai};
    fprintf('\n========== %s ==========\n', area);

    %% Load
    fpath = fullfile(DataDir, area, sprintf('PopPSTH_RF_%s.mat', area));
    D = load(fpath);

    psth_hi_s = D.psth_targeton_raw_stable_RF_high;   % {1 x N}
    psth_lo_s = D.psth_targeton_raw_stable_RF_low;
    psth_hi_d = D.psth_targeton_raw_oneDR_RF_high;
    psth_lo_d = D.psth_targeton_raw_oneDR_RF_low;
    time_s    = D.PSTH_Time_TgtOn_stable(:)';          % [1 x T_s]
    time_d    = D.PSTH_Time_TgtOn_OneDR(:)';           % [1 x T_d]

    N_neurons = numel(psth_hi_s);
    T_s = numel(time_s);
    T_d = numel(time_d);
    fprintf('  Neurons: %d\n', N_neurons);

    %% ── Find population peak window ─────────────────────────────────────────
    % Mean PSTH per neuron per condition: [N x T]
    mn = @(c) cell2mat(cellfun(@(x) mean(x, 1), c(:), 'UniformOutput', false));
    mean_hi_s = mn(psth_hi_s);   % [N x T_s]
    mean_lo_s = mn(psth_lo_s);
    mean_hi_d = mn(psth_hi_d);   % [N x T_d]
    mean_lo_d = mn(psth_lo_d);

    diff_s = mean_hi_s - mean_lo_s;   % [N x T_s] value difference timecourse
    diff_d = mean_hi_d - mean_lo_d;

    % Fixed window: 100–200 ms after target onset
    win_s = find(time_s >= 100 & time_s <= 200);
    win_d = find(time_d >= 100 & time_d <= 200);

    if isempty(win_s) || isempty(win_d)
        error('Fixed window 100–200 ms not found in time axis. Check PSTH_Time_TgtOn.');
    end

    fprintf('  Fixed window — Stable: %.0f–%.0f ms  |  OneDR: %.0f–%.0f ms\n', ...
        time_s(win_s(1)), time_s(win_s(end)), time_d(win_d(1)), time_d(win_d(end)));

    %% ── Build pseudo-population trial matrices at peak window ────────────────
    % For each neuron: average firing rate across peak window bins → scalar
    % per trial.  Stack across neurons → [n_trials x N_neurons].
    % Pseudo-population trials are independent across neurons (not simultaneous).
    win_mean = @(psth_cell, win) ...
        cellfun(@(x) mean(x(:, win), 2), psth_cell(:), 'UniformOutput', false);

    fr_hi_s = win_mean(psth_hi_s, win_s);   % {N x 1} each [n_trials_n x 1]
    fr_lo_s = win_mean(psth_lo_s, win_s);
    fr_hi_d = win_mean(psth_hi_d, win_d);
    fr_lo_d = win_mean(psth_lo_d, win_d);

    % Minimum trials per class across neurons (balance classes)
    n_min_s = min(min(cellfun(@numel, fr_hi_s)), min(cellfun(@numel, fr_lo_s)));
    n_min_d = min(min(cellfun(@numel, fr_hi_d)), min(cellfun(@numel, fr_lo_d)));
    fprintf('  Min trials/class — Stable: %d  |  OneDR: %d\n', n_min_s, n_min_d);

    % Build [n_min x N_neurons] matrices by sampling without replacement
    sample = @(v, n) v(randperm(numel(v), n));
    build_X = @(fr_cell, n_min) ...
        cell2mat(cellfun(@(v) sample(v, n_min)', fr_cell(:), 'UniformOutput', false))';
    % result: [n_min x N_neurons]

    X_hi_s = build_X(fr_hi_s, n_min_s);   % [n_min_s x N_neurons]
    X_lo_s = build_X(fr_lo_s, n_min_s);
    X_hi_d = build_X(fr_hi_d, n_min_d);
    X_lo_d = build_X(fr_lo_d, n_min_d);

    %% ── LDA: find optimal value axis per context ─────────────────────────────
    % Labels: 1 = high value, 0 = low value
    X_s = [X_hi_s; X_lo_s];
    y_s = [ones(n_min_s, 1); zeros(n_min_s, 1)];

    X_d = [X_hi_d; X_lo_d];
    y_d = [ones(n_min_d, 1); zeros(n_min_d, 1)];

    % Remove neurons with zero within-class variance in either context —
    % they are uninformative and cause fitcdiscr to fail even with Gamma > 0.
    var_s = var(X_hi_s) + var(X_lo_s);   % [1 x N_neurons] sum of within-class vars
    var_d = var(X_hi_d) + var(X_lo_d);
    keep  = (var_s > 0) & (var_d > 0);
    if any(~keep)
        fprintf('  Dropping %d zero-variance neuron(s) before LDA.\n', sum(~keep));
    end

    mdl_s = fitcdiscr(X_s(:, keep), y_s, 'DiscrimType', 'linear', 'Gamma', Gamma);
    mdl_d = fitcdiscr(X_d(:, keep), y_d, 'DiscrimType', 'linear', 'Gamma', Gamma);

    % Embed discriminant vector back into full N_neurons space (zeros for dropped)
    w_lda_s = zeros(N_neurons, 1);
    w_lda_d = zeros(N_neurons, 1);
    w_lda_s(keep) = mdl_s.Coeffs(1,2).Linear;
    w_lda_d(keep) = mdl_d.Coeffs(1,2).Linear;
    w_lda_s = w_lda_s / norm(w_lda_s);
    w_lda_d = w_lda_d / norm(w_lda_d);

    %% ── Angle between LDA axes ───────────────────────────────────────────────
    cos_theta = dot(w_lda_s, w_lda_d);
    theta_deg = acosd(cos_theta);
    fprintf('  Angle between LDA value axes: %.2f deg\n', theta_deg);

    %% ── Bootstrap CI for angle ───────────────────────────────────────────────
    N_boot = 500;
    boot_theta = zeros(1, N_boot);
    for b = 1:N_boot
        idx_s = randi(n_min_s * 2, n_min_s * 2, 1);
        idx_d = randi(n_min_d * 2, n_min_d * 2, 1);

        Xbs = X_s(idx_s, keep);   Xbd = X_d(idx_d, keep);
        % Drop any neuron that became zero-variance in this resample
        kbs = (var(Xbs(y_s(idx_s)==1,:)) + var(Xbs(y_s(idx_s)==0,:))) > 0;
        kbd = (var(Xbd(y_d(idx_d)==1,:)) + var(Xbd(y_d(idx_d)==0,:))) > 0;

        mdl_bs = fitcdiscr(Xbs(:,kbs), y_s(idx_s), 'DiscrimType','linear','Gamma',Gamma);
        mdl_bd = fitcdiscr(Xbd(:,kbd), y_d(idx_d), 'DiscrimType','linear','Gamma',Gamma);

        w_bs = zeros(N_neurons, 1);  w_bs(keep) = 0;  w_bs_k = zeros(sum(keep),1);
        w_bd = zeros(N_neurons, 1);  w_bd_k = zeros(sum(keep),1);
        w_bs_k(kbs) = mdl_bs.Coeffs(1,2).Linear;
        w_bd_k(kbd) = mdl_bd.Coeffs(1,2).Linear;
        w_bs(keep) = w_bs_k;
        w_bd(keep) = w_bd_k;
        w_bs = w_bs / norm(w_bs);
        w_bd = w_bd / norm(w_bd);

        % Sign-align to originals
        if dot(w_bs, w_lda_s) < 0, w_bs = -w_bs; end
        if dot(w_bd, w_lda_d) < 0, w_bd = -w_bd; end

        boot_theta(b) = acosd(dot(w_bs, w_bd));
    end
    ci_t = prctile(boot_theta, [5 95]);
    fprintf('  Bootstrap 90%% CI: [%.1f°, %.1f°]\n', ci_t(1), ci_t(2));

    %% ── Decoding accuracy over time ──────────────────────────────────────────
    % At each time bin t, project the per-trial firing rates onto the LDA axis
    % and classify. Threshold = midpoint between projected class means from the
    % TRAINING context (no leakage from test data).
    %
    % Within-context  (k-fold CV): axis trained on context A, CV within A
    % Cross-context               : axis trained on context A, tested on B

    % Build time-resolved trial matrices: [n_min x N_neurons x T]
    % Xt(trial, neuron, time) — explicit loop avoids dimension-order ambiguity.
    make_Xt = @(psth_cell, n_min, T) deal_Xt(psth_cell, n_min, T, N_neurons);

    Xt_hi_s = make_Xt(psth_hi_s, n_min_s, T_s);   % [n_min_s x N x T_s]
    Xt_lo_s = make_Xt(psth_lo_s, n_min_s, T_s);
    Xt_hi_d = make_Xt(psth_hi_d, n_min_d, T_d);
    Xt_lo_d = make_Xt(psth_lo_d, n_min_d, T_d);

    % Thresholds from true training-context mean PSTHs projected onto LDA axes
    thr_s_t = (w_lda_s' * mean_hi_s + w_lda_s' * mean_lo_s) / 2;   % [1 x T_s]
    thr_d_t = (w_lda_d' * mean_hi_d + w_lda_d' * mean_lo_d) / 2;   % [1 x T_d]

    acc_within_s = zeros(1, T_s);
    acc_within_d = zeros(1, T_d);
    acc_cross_s_on_d = zeros(1, T_d);   % stable axis → oneDR trials
    acc_cross_d_on_s = zeros(1, T_s);   % oneDR  axis → stable trials

    for t = 1:T_s
        X_hi = reshape(Xt_hi_s(:,:,t), n_min_s, N_neurons);   % [n_min_s x N]
        X_lo = reshape(Xt_lo_s(:,:,t), n_min_s, N_neurons);
        thr  = thr_s_t(t);
        acc_within_s(t)    = mean([(X_hi * w_lda_s) > thr; (X_lo * w_lda_s) < thr]);
        acc_cross_d_on_s(t)= mean([(X_hi * w_lda_d) > thr_d_t(min(t,T_d)); ...
                                    (X_lo * w_lda_d) < thr_d_t(min(t,T_d))]);
    end

    for t = 1:T_d
        X_hi = reshape(Xt_hi_d(:,:,t), n_min_d, N_neurons);   % [n_min_d x N]
        X_lo = reshape(Xt_lo_d(:,:,t), n_min_d, N_neurons);
        thr  = thr_d_t(t);
        acc_within_d(t)    = mean([(X_hi * w_lda_d) > thr; (X_lo * w_lda_d) < thr]);
        acc_cross_s_on_d(t)= mean([(X_hi * w_lda_s) > thr_s_t(min(t,T_s)); ...
                                    (X_lo * w_lda_s) < thr_s_t(min(t,T_s))]);
    end

    %% ── Plots ────────────────────────────────────────────────────────────────
    % Decoding accuracy
    figure('Name', sprintf('%s — LDA decoding accuracy', area), ...
           'Position', [100 100 950 420]);

    subplot(1,2,1); hold on;
    plot(time_s, acc_within_s*100, 'b-', 'LineWidth', 2, 'DisplayName', 'Stable context');
    plot(time_d, acc_within_d*100, 'r-', 'LineWidth', 2, 'DisplayName', 'OneDR context');
    yline(50,'color', 'k','LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xline(0, 'color', 'k','LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Time (ms)'); ylabel('Accuracy (%)');
    title(sprintf('%s — Within-context LDA decoding', area));
    legend('Location', 'best'); box off; ylim([40 100]);

    subplot(1,2,2); hold on;
    plot(time_d, acc_cross_s_on_d*100, 'b-', 'LineWidth', 2, 'DisplayName', 'Stable axis → OneDR');
    plot(time_s, acc_cross_d_on_s*100, 'r-', 'LineWidth', 2, 'DisplayName', 'OneDR axis → Stable');
    yline(50, 'color', 'k','LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xline(0,  'color', 'k','LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Time (ms)'); ylabel('Accuracy (%)');
    title(sprintf('%s — Cross-context LDA decoding', area));
    legend('Location', 'best'); box off; ylim([40 100]);

    % Bootstrap CI histogram for angle
    figure('Name', sprintf('%s — LDA axis angle bootstrap CI', area), ...
           'Position', [100 100 500 360]);
    hold on;
    histogram(boot_theta, 30, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');
    xline(theta_deg, 'color', 'r','LineStyle','-', 'LineWidth', 2,   'DisplayName', sprintf('Observed: %.1f°', theta_deg));
    xline(ci_t(1),   'color', 'b','LineStyle','--', 'LineWidth', 1.5, 'DisplayName', sprintf('90%% CI: [%.1f°, %.1f°]', ci_t(1), ci_t(2)));
    xline(ci_t(2),   'color', 'b','LineStyle','--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('Angle between LDA value axes (deg)');
    ylabel('Bootstrap count');
    title(sprintf('%s — LDA axis angle  \\theta=%.1f° [%.1f°, %.1f°]', area, theta_deg, ci_t(1), ci_t(2)));
    legend('Location', 'best'); box off;
end

%% ── Helper ───────────────────────────────────────────────────────────────────
function Xt = deal_Xt(psth_cell, n_min, T, N_neurons)
% Build [n_min x N_neurons x T] array from cell array of [n_trials x T] matrices.
% Randomly samples n_min trials per neuron without replacement.
Xt = zeros(n_min, N_neurons, T);
for n = 1:N_neurons
    trials = psth_cell{n};                          % [n_trials x T]
    idx    = randperm(size(trials, 1), n_min);
    Xt(:, n, :) = trials(idx, :);                  % [n_min x T] → slot n
end
end
