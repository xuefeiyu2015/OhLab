% Batch_ObjLoc_PopPSTH_FEF.m
% Load population PSTH data from FEF recordings (all modulation types).
%
% Data source: ObjectLocationPSTH/FEF/PopPSTH_RF_FEF.mat
% Variables loaded (per-trial cell arrays, one cell per neuron):
%   psth_targeton_raw_stable_RF_high    {N x 1} — Stable task, high value, RF trials
%   psth_targeton_raw_stable_RF_low     {N x 1} — Stable task, low value, RF trials
%   psth_targeton_raw_stable_nonRF_high {N x 1} — Stable task, high value, nonRF trials
%   psth_targeton_raw_stable_nonRF_low  {N x 1} — Stable task, low value, nonRF trials
%   psth_targeton_raw_oneDR_RF_high     {N x 1} — OneDR task, high value, RF trials
%   psth_targeton_raw_oneDR_RF_low      {N x 1} — OneDR task, low value, RF trials
%   psth_targeton_raw_oneDR_nonRF_high  {N x 1} — OneDR task, high value, nonRF trials
%   psth_targeton_raw_oneDR_nonRF_low   {N x 1} — OneDR task, low value, nonRF trials
%   ModuType_Sel                        [N x 1] — modulation type per neuron
%     (1=Fast Excitation, 2=Slow Excitation, 3=Inhibition, 4=No Effect)
%   PSTH_Time_TgtOn_OneDR              [1 x T] — time axis for OneDR
%   PSTH_Time_TgtOn_stable             [1 x T] — time axis for Stable
clear all

%% Decoder handle
StableDecoder = 0;
OneDRDecoder = 1;
StableToOneDR = 1;
OneDRToStable = 1;



%% Paths
DataPath   = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/ObjectLocationPSTH/FEF';
DataFile   = fullfile(DataPath, 'PopPSTH_RF_FEF.mat');
brain_area = 'FEF';   % tag saved into every result file

%% Load data (all neurons, with ModuType)
load(DataFile, ...
    'psth_targeton_raw_stable_RF_high', ...
    'psth_targeton_raw_stable_RF_low', ...
    'psth_targeton_raw_stable_nonRF_high', ...
    'psth_targeton_raw_stable_nonRF_low', ...
    'psth_targeton_raw_oneDR_RF_high', ...
    'psth_targeton_raw_oneDR_RF_low', ...
    'psth_targeton_raw_oneDR_nonRF_high', ...
    'psth_targeton_raw_oneDR_nonRF_low', ...
    'ModuType_Sel', ...
    'PSTH_Time_TgtOn_OneDR', ...
    'PSTH_Time_TgtOn_stable');

fprintf('Loaded: %s\n', DataFile);
fprintf('  Total neurons: %d\n', numel(psth_targeton_raw_stable_RF_high));
fprintf('  ModuType counts: Fast=%d, Slow=%d, Inhib=%d, NoEffect=%d\n', ...
    sum(ModuType_Sel==1), sum(ModuType_Sel==2), sum(ModuType_Sel==3), sum(ModuType_Sel==4));

%% Select neuron group (change here to filter by ModuType)
% e.g. Others only: sel = ModuType_Sel == 4;
%      Fast only:   sel = ModuType_Sel == 1;
%      All neurons: sel = true(size(ModuType_Sel));
sel =true(size(ModuType_Sel));
%sel = true(size(ModuType_Sel));

psth_high_stable_sel = psth_targeton_raw_stable_RF_high(sel);
psth_low_stable_sel  = psth_targeton_raw_stable_RF_low(sel);

psth_high_oneDR_sel = psth_targeton_raw_oneDR_RF_high(sel);
psth_low_oneDR_sel  = psth_targeton_raw_oneDR_RF_low(sel);




%% Preprocess (smooth)
[psth_high_stable_smoothed, time_full] = PreprocessPSTH(psth_high_stable_sel, PSTH_Time_TgtOn_stable, 'kernel', 'gaussian', 'win_ms', 30);
[psth_low_stable_smoothed, ~]          = PreprocessPSTH(psth_low_stable_sel,  PSTH_Time_TgtOn_stable,'kernel', 'gaussian', 'win_ms', 30);

%% Select analysis time window [-500, 300] ms
t_sel = time_full >= -500 & time_full <= 300;
time_out = time_full(t_sel);
for n = 1:numel(psth_high_stable_smoothed)
    psth_high_stable_smoothed{n} = psth_high_stable_smoothed{n}(:, t_sel);
    psth_low_stable_smoothed{n}  = psth_low_stable_smoothed{n}(:, t_sel);
end

N     = numel(psth_high_stable_smoothed);
T_out = numel(time_out);

%% Preprocess oneDR (smooth + select same time window [-500, 300] ms)
[psth_high_oneDR_smoothed, time_full_oneDR] = PreprocessPSTH(psth_high_oneDR_sel, PSTH_Time_TgtOn_OneDR, 'kernel', 'gaussian', 'win_ms', 30);
[psth_low_oneDR_smoothed, ~]               = PreprocessPSTH(psth_low_oneDR_sel,   PSTH_Time_TgtOn_OneDR, 'kernel', 'gaussian', 'win_ms', 30);

t_sel_oneDR    = time_full_oneDR >= -500 & time_full_oneDR <= 300;
time_out_oneDR = time_full_oneDR(t_sel_oneDR);
for n = 1:numel(psth_high_oneDR_smoothed)
    psth_high_oneDR_smoothed{n} = psth_high_oneDR_smoothed{n}(:, t_sel_oneDR);
    psth_low_oneDR_smoothed{n}  = psth_low_oneDR_smoothed{n}(:, t_sel_oneDR);
end
T_out_oneDR = numel(time_out_oneDR);

% n_min per neuron for oneDR: equal trial count for high and low
n_min_oneDR = zeros(N, 1);
for n = 1:N
    n_min_oneDR(n) = min(size(psth_high_oneDR_smoothed{n}, 1), ...
                         size(psth_low_oneDR_smoothed{n},  1));
end










%% Pseudo-population decoder — stable task, RF trials
% Repeated stratified k-fold CV per neuron.
% Each neuron's high/low trials are independently split into k folds.
% Within each training fold, trials are resampled to equalise trial counts
% across neurons. The test fold uses separate held-out trials (also resampled
% to keep matrix dimensions fixed). Accuracy is averaged across k folds to
% give one estimate per outer repeat, reducing Monte Carlo variance.
k_folds         = 5;   % number of CV folds per neuron
n_resample      = 10;  % resampled trials per condition per neuron (train, with replacement)
n_test_resample = 10;  % resampled trials per condition per neuron (test, with replacement)
n_repeats       = 50;  % outer repetitions of the full k-fold
n_perm          = 50;  % permutation iterations

acc_all     = zeros(n_repeats, T_out);         % [repeats x T]
weights_all = zeros(n_repeats, T_out, N);      % [repeats x T x N]
acc_cross   = zeros(n_repeats, T_out, T_out);  % [repeats x T_train x T_test]


if StableDecoder

% Trial counts per neuron (high and low kept separate for per-class fold assignment)
n_hi_stable = zeros(N, 1);
n_lo_stable = zeros(N, 1);
for n = 1:N
    n_hi_stable(n) = size(psth_high_stable_smoothed{n}, 1);
    n_lo_stable(n) = size(psth_low_stable_smoothed{n},  1);
end

% Pre-compute fold assignments for all repeats in serial (local script
% functions are not visible to parfor workers, so we generate them here).
pf_hi = cell(n_repeats, N);   % pf_hi{rep, n}: fold labels for high trials of neuron n in repeat rep
pf_lo = cell(n_repeats, N);
for rep_pre = 1:n_repeats
    for n = 1:N
        pf_hi{rep_pre, n} = AssignFolds(n_hi_stable(n), k_folds);
        pf_lo{rep_pre, n} = AssignFolds(n_lo_stable(n), k_folds);
    end
end

parfor rep = 1:n_repeats
    % --- Step 1: Retrieve pre-computed fold assignments for this repeat ---
    fold_hi = pf_hi(rep, :);   % {1 x N} cell — fold label per high trial
    fold_lo = pf_lo(rep, :);   % {1 x N} cell — fold label per low trial

    % Accumulators across folds (averaged at the end)
    rep_acc_sum   = zeros(1, T_out);
    rep_w_sum     = zeros(T_out, N);
    rep_cross_sum = zeros(T_out, T_out);

    for f = 1:k_folds
        % --- Step 2: Assemble train and test pseudo-populations for fold f ---
        X_tr_all = zeros(2*n_resample,      N, T_out, 'single');
        X_te_all = zeros(2*n_test_resample, N, T_out, 'single');

        for n = 1:N
            hi_data = psth_high_stable_smoothed{n};
            lo_data = psth_low_stable_smoothed{n};

            hi_tr_idx = find(fold_hi{n} ~= f);
            hi_te_idx = find(fold_hi{n} == f);
            lo_tr_idx = find(fold_lo{n} ~= f);
            lo_te_idx = find(fold_lo{n} == f);

            % Resample train trials with replacement
            rs_hi_tr = hi_tr_idx(randi(numel(hi_tr_idx), n_resample,      1));
            rs_lo_tr = lo_tr_idx(randi(numel(lo_tr_idx), n_resample,      1));
            X_tr_all(:, n, :) = single([hi_data(rs_hi_tr, :); lo_data(rs_lo_tr, :)]);

            % Resample test trials with replacement (held-out fold only)
            rs_hi_te = hi_te_idx(randi(numel(hi_te_idx), n_test_resample, 1));
            rs_lo_te = lo_te_idx(randi(numel(lo_te_idx), n_test_resample, 1));
            X_te_all(:, n, :) = single([hi_data(rs_hi_te, :); lo_data(rs_lo_te, :)]);
        end

        y_tr = [ones(n_resample,      1); zeros(n_resample,      1)];
        y_te = [ones(n_test_resample, 1); zeros(n_test_resample, 1)];

        % --- Step 3: Z-score using training statistics only ---
        mu = mean(X_tr_all, 1);
        sd = std(X_tr_all, 0, 1) + 1e-8;
        X_tr_z = (X_tr_all - mu) ./ sd;
        X_te_z = (X_te_all - mu) ./ sd;

        % --- Step 4: Train at each time bin, test on all time bins ---
        fold_acc   = zeros(1, T_out);
        fold_w     = zeros(T_out, N);
        fold_cross = zeros(T_out, T_out);

        for t_train = 1:T_out
            X_tr = squeeze(X_tr_z(:, :, t_train));
            mdl  = fitclinear(X_tr, y_tr, 'Learner', 'logistic', ...
                'Regularization', 'ridge', 'Lambda', 'auto', 'Solver', 'bfgs');
            fold_w(t_train, :)    = mdl.Beta';
            fold_acc(t_train)     = mean(predict(mdl, squeeze(X_te_z(:, :, t_train))) == y_te);

            % Cross-temporal generalization
            for t_test = 1:T_out
                fold_cross(t_train, t_test) = mean(predict(mdl, squeeze(X_te_z(:, :, t_test))) == y_te);
            end
        end

        rep_acc_sum   = rep_acc_sum   + fold_acc;
        rep_w_sum     = rep_w_sum     + fold_w;
        rep_cross_sum = rep_cross_sum + fold_cross;
    end

    % --- Step 5: Average across folds ---
    acc_all(rep, :)        = rep_acc_sum   / k_folds;
    weights_all(rep, :, :) = rep_w_sum     / k_folds;
    acc_cross(rep, :, :)   = rep_cross_sum / k_folds;
    fprintf('Stable rep %d/%d done\n', rep, n_repeats);
end

%% Permutation test — shuffle training labels, same k-fold structure
% Each permutation iteration draws one random 80/20 split per neuron
% (consistent with the k=5 fold fractions) to build train and test sets,
% then shuffles the training labels before fitting.
acc_perm = zeros(n_perm, T_out);

parfor perm = 1:n_perm
    X_tr_all = zeros(2*n_resample,      N, T_out, 'single');
    X_te_all = zeros(2*n_test_resample, N, T_out, 'single');

    for n = 1:N
        hi_data = psth_high_stable_smoothed{n};
        lo_data = psth_low_stable_smoothed{n};

        % Random 80/20 split per class (mirrors one k=5 fold held-out)
        perm_hi   = randperm(n_hi_stable(n));
        n_tr_hi   = round(n_hi_stable(n) * 0.8);
        hi_tr_idx = perm_hi(1:n_tr_hi);
        hi_te_idx = perm_hi(n_tr_hi+1:end);

        perm_lo   = randperm(n_lo_stable(n));
        n_tr_lo   = round(n_lo_stable(n) * 0.8);
        lo_tr_idx = perm_lo(1:n_tr_lo);
        lo_te_idx = perm_lo(n_tr_lo+1:end);

        rs_hi_tr = hi_tr_idx(randi(numel(hi_tr_idx), n_resample,      1));
        rs_lo_tr = lo_tr_idx(randi(numel(lo_tr_idx), n_resample,      1));
        X_tr_all(:, n, :) = single([hi_data(rs_hi_tr, :); lo_data(rs_lo_tr, :)]);

        rs_hi_te = hi_te_idx(randi(numel(hi_te_idx), n_test_resample, 1));
        rs_lo_te = lo_te_idx(randi(numel(lo_te_idx), n_test_resample, 1));
        X_te_all(:, n, :) = single([hi_data(rs_hi_te, :); lo_data(rs_lo_te, :)]);
    end

    % Shuffle training labels only — test labels stay intact
    y_shuf = [ones(n_resample,     1); zeros(n_resample,     1)];
    y_shuf = y_shuf(randperm(2*n_resample));
    y_true = [ones(n_test_resample,1); zeros(n_test_resample,1)];

    mu = mean(X_tr_all, 1);
    sd = std(X_tr_all, 0, 1) + 1e-8;
    X_tr_z = (X_tr_all - mu) ./ sd;
    X_te_z = (X_te_all - mu) ./ sd;

    perm_acc = zeros(1, T_out);
    for t = 1:T_out
        X_tr = squeeze(X_tr_z(:, :, t));
        X_te = squeeze(X_te_z(:, :, t));
        mdl = fitclinear(X_tr, y_shuf, 'Learner', 'logistic', ...
            'Regularization', 'ridge', 'Lambda', 'auto', 'Solver', 'bfgs');
        perm_acc(t) = mean(predict(mdl, X_te) == y_true);
    end
    acc_perm(perm, :) = perm_acc;
    if mod(perm, 10) == 0
        fprintf('Perm %d/%d done\n', perm, n_perm);
    end
end

%% Save
SaveFile_Dec = fullfile(DataPath, 'Decoder_Stable_RF_FEF_Others.mat');
save(SaveFile_Dec, 'acc_all', 'weights_all', 'acc_perm', 'acc_cross', 'time_out', ...
    'n_resample', 'n_test_resample', 'n_repeats', 'k_folds', 'ModuType_Sel', 'n_perm', 'brain_area', '-v7.3');
fprintf('Saved -> %s\n', SaveFile_Dec);

%% Significance test (cluster-based permutation)
mean_acc  = mean(acc_all,  1);
sem_acc   = std(acc_all,  0, 1) / sqrt(n_repeats);
mean_perm = mean(acc_perm, 1);
sem_perm  = std(acc_perm, 0, 1) / sqrt(n_perm);

% Per-bin p-value: proportion of permutation means >= real mean
p_val = zeros(1, T_out);
for t = 1:T_out
    p_val(t) = mean(acc_perm(:, t) >= mean_acc(t));
end

% Cluster-based correction for multiple comparisons
alpha = 0.05;
sig_bins = p_val < alpha;   % uncorrected significant bins

% Find clusters of contiguous significant bins
clusters = bwlabel(sig_bins);
n_clusters = max(clusters);

% Cluster mass: sum of (1 - p_val) within each cluster
cluster_mass = zeros(1, n_clusters);
for c = 1:n_clusters
    cluster_mass(c) = sum(1 - p_val(clusters == c));
end

% Null cluster mass distribution: for each permutation, compute max cluster mass
null_cluster_mass = zeros(n_perm, 1);
for perm = 1:n_perm
    % p-value of this permutation vs all others
    p_perm_t = zeros(1, T_out);
    for t = 1:T_out
        p_perm_t(t) = mean(acc_perm(:, t) >= acc_perm(perm, t));
    end
    sig_perm = p_perm_t < alpha;
    cl_perm  = bwlabel(sig_perm);
    n_cl     = max(cl_perm);
    if n_cl > 0
        masses = zeros(1, n_cl);
        for c = 1:n_cl
            masses(c) = sum(1 - p_perm_t(cl_perm == c));
        end
        null_cluster_mass(perm) = max(masses);
    end
end

% Keep only clusters with mass exceeding 95th percentile of null
mass_thresh = prctile(null_cluster_mass, 95);
sig_corrected = false(1, T_out);
for c = 1:n_clusters
    if cluster_mass(c) > mass_thresh
        sig_corrected(clusters == c) = true;
    end
end

fprintf('Significant time bins (cluster-corrected): %d / %d\n', sum(sig_corrected), T_out);

%% Plot
figure; hold on;
fill([time_out, fliplr(time_out)], ...
     [mean_perm + sem_perm, fliplr(mean_perm - sem_perm)], ...
     [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(time_out, mean_perm, 'k-', 'LineWidth', 1);
fill([time_out, fliplr(time_out)], ...
     [mean_acc + sem_acc, fliplr(mean_acc - sem_acc)], ...
     [0.3 0.6 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(time_out, mean_acc, 'b-', 'LineWidth', 1.5);

% Mark significant time bins (cluster-corrected) as thick horizontal bar
y_sig = 0.95;  % fixed y-position within ylim
dt = time_out(2) - time_out(1);
sig_idx = find(sig_corrected);
for ii = 1:numel(sig_idx)
    t_c = time_out(sig_idx(ii));
    fill([t_c - dt/2, t_c + dt/2, t_c + dt/2, t_c - dt/2], ...
         [y_sig - 0.01, y_sig - 0.01, y_sig + 0.01, y_sig + 0.01], ...
         'b', 'EdgeColor', 'none');
end

yline(0.5, 'LineStyle','--','Color','k');
xline(0,  'LineStyle','--', 'Color','k');
xlabel('Time (ms)');
ylabel('Decoding accuracy');
title('Value decoder — Stable task, RF, FEF');
xlim([time_out(1), time_out(end)]);
ylim([0.3 1]);
box off;


%
%% Cross-temporal generalization heatmap
mean_cross = squeeze(mean(acc_cross, 1));   % [T_train x T_test]

figure;
imagesc(time_out, time_out, mean_cross);
axis xy;
hold on;
xline(0, 'LineStyle','--','Color','w', 'LineWidth', 1);
yline(0, 'LineStyle','--','Color','w', 'LineWidth', 1);
colorbar;
%caxis([0.4 0.8]);
colormap('jet');
xlabel('Test time (ms)');
ylabel('Train time (ms)');
title('Cross-temporal generalization — Stable, RF, FEF');
axis square;

%keyboard;
end

if OneDRDecoder
%% ============================================================
%% OneDR task decoder — RF trials
%% Repeated stratified k-fold CV (same structure as stable decoder).
%% ============================================================

% Per-neuron trial counts for OneDR (high and low kept separate)
n_hi_oneDR = zeros(N, 1);
n_lo_oneDR = zeros(N, 1);
for n = 1:N
    n_hi_oneDR(n) = size(psth_high_oneDR_smoothed{n}, 1);
    n_lo_oneDR(n) = size(psth_low_oneDR_smoothed{n},  1);
end

% Pre-compute fold assignments for all repeats in serial
pf_hi_oneDR = cell(n_repeats, N);
pf_lo_oneDR = cell(n_repeats, N);
for rep_pre = 1:n_repeats
    for n = 1:N
        pf_hi_oneDR{rep_pre, n} = AssignFolds(n_hi_oneDR(n), k_folds);
        pf_lo_oneDR{rep_pre, n} = AssignFolds(n_lo_oneDR(n), k_folds);
    end
end

acc_all_oneDR     = zeros(n_repeats, T_out_oneDR);
weights_all_oneDR = zeros(n_repeats, T_out_oneDR, N);
acc_cross_oneDR   = zeros(n_repeats, T_out_oneDR, T_out_oneDR);

parfor rep = 1:n_repeats
    fold_hi = pf_hi_oneDR(rep, :);
    fold_lo = pf_lo_oneDR(rep, :);

    rep_acc_sum   = zeros(1, T_out_oneDR);
    rep_w_sum     = zeros(T_out_oneDR, N);
    rep_cross_sum = zeros(T_out_oneDR, T_out_oneDR);

    for f = 1:k_folds
        X_tr_all = zeros(2*n_resample,      N, T_out_oneDR, 'single');
        X_te_all = zeros(2*n_test_resample, N, T_out_oneDR, 'single');

        for n = 1:N
            hi_data = psth_high_oneDR_smoothed{n};
            lo_data = psth_low_oneDR_smoothed{n};

            hi_tr_idx = find(fold_hi{n} ~= f);
            hi_te_idx = find(fold_hi{n} == f);
            lo_tr_idx = find(fold_lo{n} ~= f);
            lo_te_idx = find(fold_lo{n} == f);

            rs_hi_tr = hi_tr_idx(randi(numel(hi_tr_idx), n_resample,      1));
            rs_lo_tr = lo_tr_idx(randi(numel(lo_tr_idx), n_resample,      1));
            X_tr_all(:, n, :) = single([hi_data(rs_hi_tr, :); lo_data(rs_lo_tr, :)]);

            rs_hi_te = hi_te_idx(randi(numel(hi_te_idx), n_test_resample, 1));
            rs_lo_te = lo_te_idx(randi(numel(lo_te_idx), n_test_resample, 1));
            X_te_all(:, n, :) = single([hi_data(rs_hi_te, :); lo_data(rs_lo_te, :)]);
        end

        y_tr = [ones(n_resample,      1); zeros(n_resample,      1)];
        y_te = [ones(n_test_resample, 1); zeros(n_test_resample, 1)];

        mu = mean(X_tr_all, 1);
        sd = std(X_tr_all, 0, 1) + 1e-8;
        X_tr_z = (X_tr_all - mu) ./ sd;
        X_te_z = (X_te_all - mu) ./ sd;

        fold_acc   = zeros(1, T_out_oneDR);
        fold_w     = zeros(T_out_oneDR, N);
        fold_cross = zeros(T_out_oneDR, T_out_oneDR);

        for t_train = 1:T_out_oneDR
            X_tr = squeeze(X_tr_z(:, :, t_train));
            mdl  = fitclinear(X_tr, y_tr, 'Learner', 'logistic', ...
                'Regularization', 'ridge', 'Lambda', 'auto', 'Solver', 'bfgs');
            fold_w(t_train, :) = mdl.Beta';
            fold_acc(t_train)  = mean(predict(mdl, squeeze(X_te_z(:, :, t_train))) == y_te);
            for t_test = 1:T_out_oneDR
                fold_cross(t_train, t_test) = mean(predict(mdl, squeeze(X_te_z(:, :, t_test))) == y_te);
            end
        end

        rep_acc_sum   = rep_acc_sum   + fold_acc;
        rep_w_sum     = rep_w_sum     + fold_w;
        rep_cross_sum = rep_cross_sum + fold_cross;
    end

    acc_all_oneDR(rep, :)        = rep_acc_sum   / k_folds;
    weights_all_oneDR(rep, :, :) = rep_w_sum     / k_folds;
    acc_cross_oneDR(rep, :, :)   = rep_cross_sum / k_folds;
    fprintf('OneDR rep %d/%d done\n', rep, n_repeats);
end

%% OneDR permutation test — shuffle training labels, same k-fold structure
acc_perm_oneDR = zeros(n_perm, T_out_oneDR);

parfor perm = 1:n_perm
    X_tr_all = zeros(2*n_resample,      N, T_out_oneDR, 'single');
    X_te_all = zeros(2*n_test_resample, N, T_out_oneDR, 'single');

    for n = 1:N
        hi_data = psth_high_oneDR_smoothed{n};
        lo_data = psth_low_oneDR_smoothed{n};

        perm_hi   = randperm(n_hi_oneDR(n));
        n_tr_hi   = round(n_hi_oneDR(n) * 0.8);
        hi_tr_idx = perm_hi(1:n_tr_hi);
        hi_te_idx = perm_hi(n_tr_hi+1:end);

        perm_lo   = randperm(n_lo_oneDR(n));
        n_tr_lo   = round(n_lo_oneDR(n) * 0.8);
        lo_tr_idx = perm_lo(1:n_tr_lo);
        lo_te_idx = perm_lo(n_tr_lo+1:end);

        rs_hi_tr = hi_tr_idx(randi(numel(hi_tr_idx), n_resample,      1));
        rs_lo_tr = lo_tr_idx(randi(numel(lo_tr_idx), n_resample,      1));
        X_tr_all(:, n, :) = single([hi_data(rs_hi_tr, :); lo_data(rs_lo_tr, :)]);

        rs_hi_te = hi_te_idx(randi(numel(hi_te_idx), n_test_resample, 1));
        rs_lo_te = lo_te_idx(randi(numel(lo_te_idx), n_test_resample, 1));
        X_te_all(:, n, :) = single([hi_data(rs_hi_te, :); lo_data(rs_lo_te, :)]);
    end

    y_shuf = [ones(n_resample,     1); zeros(n_resample,     1)];
    y_shuf = y_shuf(randperm(2*n_resample));
    y_true = [ones(n_test_resample,1); zeros(n_test_resample,1)];

    mu = mean(X_tr_all, 1);
    sd = std(X_tr_all, 0, 1) + 1e-8;
    X_tr_z = (X_tr_all - mu) ./ sd;
    X_te_z = (X_te_all - mu) ./ sd;

    perm_acc = zeros(1, T_out_oneDR);
    for t = 1:T_out_oneDR
        X_tr = squeeze(X_tr_z(:, :, t));
        X_te = squeeze(X_te_z(:, :, t));
        mdl  = fitclinear(X_tr, y_shuf, 'Learner', 'logistic', ...
            'Regularization', 'ridge', 'Lambda', 'auto', 'Solver', 'bfgs');
        perm_acc(t) = mean(predict(mdl, X_te) == y_true);
    end
    acc_perm_oneDR(perm, :) = perm_acc;
    if mod(perm, 10) == 0
        fprintf('OneDR Perm %d/%d done\n', perm, n_perm);
    end
end

%% Save OneDR decoder
SaveFile_OneDR = fullfile(DataPath, 'Decoder_OneDR_RF_FEF_Others.mat');
save(SaveFile_OneDR, 'acc_all_oneDR', 'weights_all_oneDR', 'acc_perm_oneDR', ...
    'acc_cross_oneDR', 'time_out_oneDR', ...
    'n_resample', 'n_repeats', 'n_perm', 'brain_area', '-v7.3');
fprintf('Saved -> %s\n', SaveFile_OneDR);

%% Significance test + plot: OneDR decoder
mean_acc_oneDR  = mean(acc_all_oneDR,  1);
sem_acc_oneDR   = std(acc_all_oneDR,  0, 1) / sqrt(n_repeats);
mean_perm_oneDR = mean(acc_perm_oneDR, 1);
sem_perm_oneDR  = std(acc_perm_oneDR, 0, 1) / sqrt(n_perm);

p_val_oneDR = zeros(1, T_out_oneDR);
for t = 1:T_out_oneDR
    p_val_oneDR(t) = mean(acc_perm_oneDR(:, t) >= mean_acc_oneDR(t));
end

alpha_oneDR      = 0.05;
sig_bins_oneDR   = p_val_oneDR < alpha_oneDR;
clusters_oneDR   = bwlabel(sig_bins_oneDR);
n_clusters_oneDR = max(clusters_oneDR);

cluster_mass_oneDR = zeros(1, n_clusters_oneDR);
for c = 1:n_clusters_oneDR
    cluster_mass_oneDR(c) = sum(1 - p_val_oneDR(clusters_oneDR == c));
end

null_cluster_mass_oneDR = zeros(n_perm, 1);
for perm = 1:n_perm
    p_perm_t = zeros(1, T_out_oneDR);
    for t = 1:T_out_oneDR
        p_perm_t(t) = mean(acc_perm_oneDR(:, t) >= acc_perm_oneDR(perm, t));
    end
    sig_perm = p_perm_t < alpha_oneDR;
    cl_perm  = bwlabel(sig_perm);
    n_cl     = max(cl_perm);
    if n_cl > 0
        masses = zeros(1, n_cl);
        for c = 1:n_cl
            masses(c) = sum(1 - p_perm_t(cl_perm == c));
        end
        null_cluster_mass_oneDR(perm) = max(masses);
    end
end

mass_thresh_oneDR   = prctile(null_cluster_mass_oneDR, 95);
sig_corrected_oneDR = false(1, T_out_oneDR);
for c = 1:n_clusters_oneDR
    if cluster_mass_oneDR(c) > mass_thresh_oneDR
        sig_corrected_oneDR(clusters_oneDR == c) = true;
    end
end

fprintf('OneDR significant bins (cluster-corrected): %d / %d\n', sum(sig_corrected_oneDR), T_out_oneDR);

figure; hold on;
fill([time_out_oneDR, fliplr(time_out_oneDR)], ...
     [mean_perm_oneDR + sem_perm_oneDR, fliplr(mean_perm_oneDR - sem_perm_oneDR)], ...
     [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(time_out_oneDR, mean_perm_oneDR, 'k-', 'LineWidth', 1);
fill([time_out_oneDR, fliplr(time_out_oneDR)], ...
     [mean_acc_oneDR + sem_acc_oneDR, fliplr(mean_acc_oneDR - sem_acc_oneDR)], ...
     [1 0.5 0.2], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(time_out_oneDR, mean_acc_oneDR, 'Color', [0.9 0.4 0.1], 'LineWidth', 1.5);

dt_oneDR      = time_out_oneDR(2) - time_out_oneDR(1);
sig_idx_oneDR = find(sig_corrected_oneDR);
for ii = 1:numel(sig_idx_oneDR)
    t_c = time_out_oneDR(sig_idx_oneDR(ii));
    fill([t_c - dt_oneDR/2, t_c + dt_oneDR/2, t_c + dt_oneDR/2, t_c - dt_oneDR/2], ...
         [0.94, 0.94, 0.96, 0.96], [0.9 0.4 0.1], 'EdgeColor', 'none');
end

yline(0.5, 'LineStyle','--','Color','k');
xline(0,   'LineStyle','--','Color','k');
xlabel('Time (ms)');
ylabel('Decoding accuracy');
title('Value decoder — OneDR task, RF, FEF');
xlim([time_out_oneDR(1), time_out_oneDR(end)]);
ylim([0.3 1]);
box off;

mean_cross_oneDR = squeeze(mean(acc_cross_oneDR, 1));
figure;
imagesc(time_out_oneDR, time_out_oneDR, mean_cross_oneDR);
axis xy; hold on;
xline(0, 'LineStyle','--','Color','w','LineWidth',1);
yline(0, 'LineStyle','--','Color','w','LineWidth',1);
colorbar; colormap('jet');
xlabel('Test time (ms)');
ylabel('Train time (ms)');
title('Cross-temporal generalization — OneDR, RF, FEF');
axis square;


end
%% ============================================================
%% Cross-context decoder: Train Stable -> Test OneDR
%% Repeated stratified k-fold CV on stable (training) context.
%% Z-score from stable training fold stats applied to oneDR test set.
%% ============================================================

if StableToOneDR

assert(T_out == T_out_oneDR, ...
    'Stable (%d bins) and OneDR (%d bins) must match for cross-decoding.', ...
    T_out, T_out_oneDR);

% Trial counts for stable training context (compute if not already available)
if ~exist('n_hi_stable', 'var')
    n_hi_stable = zeros(N, 1);
    n_lo_stable = zeros(N, 1);
    for n = 1:N
        n_hi_stable(n) = size(psth_high_stable_smoothed{n}, 1);
        n_lo_stable(n) = size(psth_low_stable_smoothed{n},  1);
    end
end

% Pre-compute fold assignments for stable training context in serial.
% Local functions are not visible to parfor workers, so folds must be
% assigned here before the parfor loop.
if ~exist('pf_hi', 'var')
    pf_hi = cell(n_repeats, N);
    pf_lo = cell(n_repeats, N);
    for rep_pre = 1:n_repeats
        for n = 1:N
            pf_hi{rep_pre, n} = AssignFolds(n_hi_stable(n), k_folds);
            pf_lo{rep_pre, n} = AssignFolds(n_lo_stable(n), k_folds);
        end
    end
end

acc_s2d    = zeros(n_repeats, T_out);
acc_ct_s2d = zeros(n_repeats, T_out, T_out);

parfor rep = 1:n_repeats
    fold_hi = pf_hi(rep, :);   % {1 x N} fold labels for stable high trials
    fold_lo = pf_lo(rep, :);   % {1 x N} fold labels for stable low trials

    rep_acc_sum = zeros(1, T_out);
    rep_ct_sum  = zeros(T_out, T_out);

    for f = 1:k_folds
        % --- Training: Stable, folds ~= f (resampled) ---
        X_tr_stable = zeros(2*n_resample, N, T_out, 'single');
        for n = 1:N
            hi_data   = psth_high_stable_smoothed{n};
            lo_data   = psth_low_stable_smoothed{n};
            hi_tr_idx = find(fold_hi{n} ~= f);
            lo_tr_idx = find(fold_lo{n} ~= f);
            rs_hi = hi_tr_idx(randi(numel(hi_tr_idx), n_resample, 1));
            rs_lo = lo_tr_idx(randi(numel(lo_tr_idx), n_resample, 1));
            X_tr_stable(:, n, :) = single([hi_data(rs_hi, :); lo_data(rs_lo, :)]);
        end

        % --- Test: OneDR, all trials (resampled) ---
        X_te_oneDR = zeros(2*n_test_resample, N, T_out, 'single');
        for n = 1:N
            hi_data = psth_high_oneDR_smoothed{n};
            lo_data = psth_low_oneDR_smoothed{n};
            idx_hi  = randi(size(hi_data, 1), n_test_resample, 1);
            idx_lo  = randi(size(lo_data, 1), n_test_resample, 1);
            X_te_oneDR(:, n, :) = single([hi_data(idx_hi, :); lo_data(idx_lo, :)]);
        end

        y_tr = [ones(n_resample,      1); zeros(n_resample,      1)];
        y_te = [ones(n_test_resample, 1); zeros(n_test_resample, 1)];

        % Z-score from stable training fold stats; apply same transform to oneDR
        mu = mean(X_tr_stable, 1);
        sd = std(X_tr_stable, 0, 1) + 1e-8;
        X_tr_z = (X_tr_stable - mu) ./ sd;
        X_te_z = (X_te_oneDR  - mu) ./ sd;

        fold_acc = zeros(1, T_out);
        fold_ct  = zeros(T_out, T_out);

        for t_train = 1:T_out
            X_tr = squeeze(X_tr_z(:, :, t_train));
            mdl  = fitclinear(X_tr, y_tr, 'Learner', 'logistic', ...
                'Regularization', 'ridge', 'Lambda', 'auto', 'Solver', 'bfgs');
            fold_acc(t_train) = mean(predict(mdl, squeeze(X_te_z(:, :, t_train))) == y_te);
            for t_test = 1:T_out
                fold_ct(t_train, t_test) = mean(predict(mdl, squeeze(X_te_z(:, :, t_test))) == y_te);
            end
        end

        rep_acc_sum = rep_acc_sum + fold_acc;
        rep_ct_sum  = rep_ct_sum  + fold_ct;
    end

    acc_s2d(rep, :)       = rep_acc_sum / k_folds;
    acc_ct_s2d(rep, :, :) = rep_ct_sum  / k_folds;
    fprintf('Cross S->D Rep %d/%d done\n', rep, n_repeats);
end

%% Permutation test: Train Stable -> Test OneDR (shuffle training labels)
% One random 80/20 split per permutation (mirrors one k=5 fold held-out).
acc_perm_s2d = zeros(n_perm, T_out);

parfor perm = 1:n_perm
    X_tr_stable = zeros(2*n_resample,      N, T_out, 'single');
    X_te_oneDR  = zeros(2*n_test_resample, N, T_out, 'single');

    for n = 1:N
        hi_data   = psth_high_stable_smoothed{n};
        lo_data   = psth_low_stable_smoothed{n};
        perm_hi   = randperm(n_hi_stable(n));
        n_tr_hi   = round(n_hi_stable(n) * 0.8);
        hi_tr_idx = perm_hi(1:n_tr_hi);
        perm_lo   = randperm(n_lo_stable(n));
        n_tr_lo   = round(n_lo_stable(n) * 0.8);
        lo_tr_idx = perm_lo(1:n_tr_lo);
        rs_hi = hi_tr_idx(randi(n_tr_hi, n_resample, 1));
        rs_lo = lo_tr_idx(randi(n_tr_lo, n_resample, 1));
        X_tr_stable(:, n, :) = single([hi_data(rs_hi, :); lo_data(rs_lo, :)]);

        hi_d   = psth_high_oneDR_smoothed{n};
        lo_d   = psth_low_oneDR_smoothed{n};
        idx_hi = randi(size(hi_d, 1), n_test_resample, 1);
        idx_lo = randi(size(lo_d, 1), n_test_resample, 1);
        X_te_oneDR(:, n, :) = single([hi_d(idx_hi, :); lo_d(idx_lo, :)]);
    end

    y_shuf = [ones(n_resample,      1); zeros(n_resample,      1)];
    y_shuf = y_shuf(randperm(2*n_resample));
    y_true = [ones(n_test_resample, 1); zeros(n_test_resample, 1)];

    mu = mean(X_tr_stable, 1);
    sd = std(X_tr_stable, 0, 1) + 1e-8;
    X_tr_z = (X_tr_stable - mu) ./ sd;
    X_te_z = (X_te_oneDR  - mu) ./ sd;

    perm_acc = zeros(1, T_out);
    for t = 1:T_out
        mdl = fitclinear(squeeze(X_tr_z(:,:,t)), y_shuf, 'Learner', 'logistic', ...
            'Regularization', 'ridge', 'Lambda', 'auto', 'Solver', 'bfgs');
        perm_acc(t) = mean(predict(mdl, squeeze(X_te_z(:,:,t))) == y_true);
    end
    acc_perm_s2d(perm, :) = perm_acc;
    if mod(perm, 10) == 0
        fprintf('Cross S->D Perm %d/%d done\n', perm, n_perm);
    end
end

%% Save cross stable->oneDR
SaveFile_S2D = fullfile(DataPath, 'Decoder_Cross_Stable2OneDR_RF_FEF.mat');
save(SaveFile_S2D, 'acc_s2d', 'acc_ct_s2d', 'acc_perm_s2d', 'time_out', ...
    'n_resample', 'n_test_resample', 'n_repeats', 'k_folds', 'n_perm', 'brain_area', '-v7.3');
fprintf('Saved -> %s\n', SaveFile_S2D);

%% Significance test: cross stable->oneDR (cluster correction)
mean_s2d      = mean(acc_s2d,      1);
sem_s2d       = std(acc_s2d,      0, 1) / sqrt(n_repeats);
mean_perm_s2d = mean(acc_perm_s2d, 1);
sem_perm_s2d  = std(acc_perm_s2d, 0, 1) / sqrt(n_perm);

p_val_s2d = zeros(1, T_out);
for t = 1:T_out
    p_val_s2d(t) = mean(acc_perm_s2d(:, t) >= mean_s2d(t));
end

alpha_s2d      = 0.05;
sig_bins_s2d   = p_val_s2d < alpha_s2d;
clusters_s2d   = bwlabel(sig_bins_s2d);
n_clusters_s2d = max(clusters_s2d);

cluster_mass_s2d = zeros(1, n_clusters_s2d);
for c = 1:n_clusters_s2d
    cluster_mass_s2d(c) = sum(1 - p_val_s2d(clusters_s2d == c));
end

null_cluster_mass_s2d = zeros(n_perm, 1);
for perm = 1:n_perm
    p_perm_t = zeros(1, T_out);
    for t = 1:T_out
        p_perm_t(t) = mean(acc_perm_s2d(:, t) >= acc_perm_s2d(perm, t));
    end
    sig_perm   = p_perm_t < alpha_s2d;
    cl_perm    = bwlabel(sig_perm);
    n_cl       = max(cl_perm);
    if n_cl > 0
        masses = zeros(1, n_cl);
        for c = 1:n_cl
            masses(c) = sum(1 - p_perm_t(cl_perm == c));
        end
        null_cluster_mass_s2d(perm) = max(masses);
    end
end

mass_thresh_s2d   = prctile(null_cluster_mass_s2d, 95);
sig_corrected_s2d = false(1, T_out);
for c = 1:n_clusters_s2d
    if cluster_mass_s2d(c) > mass_thresh_s2d
        sig_corrected_s2d(clusters_s2d == c) = true;
    end
end

fprintf('Cross S->D significant bins (cluster-corrected): %d / %d\n', sum(sig_corrected_s2d), T_out);

%% Plot cross stable->oneDR
dt_s2d      = time_out(2) - time_out(1);
sig_idx_s2d = find(sig_corrected_s2d);

figure; hold on;
fill([time_out, fliplr(time_out)], ...
     [mean_perm_s2d + sem_perm_s2d, fliplr(mean_perm_s2d - sem_perm_s2d)], ...
     [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(time_out, mean_perm_s2d, 'k-', 'LineWidth', 1);
fill([time_out, fliplr(time_out)], ...
     [mean_s2d + sem_s2d, fliplr(mean_s2d - sem_s2d)], ...
     [0.3 0.8 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(time_out, mean_s2d, 'Color', [0.1 0.6 0.3], 'LineWidth', 1.5);
for ii = 1:numel(sig_idx_s2d)
    t_c = time_out(sig_idx_s2d(ii));
    fill([t_c - dt_s2d/2, t_c + dt_s2d/2, t_c + dt_s2d/2, t_c - dt_s2d/2], ...
         [0.94, 0.94, 0.96, 0.96], [0.1 0.6 0.3], 'EdgeColor', 'none');
end
yline(0.5, 'LineStyle','--','Color','k');
xline(0,   'LineStyle','--','Color','k');
xlabel('Time (ms)');
ylabel('Decoding accuracy');
title('Cross-context: Train Stable \rightarrow Test OneDR, RF, FEF');
xlim([time_out(1), time_out(end)]);
ylim([0.3 1]);
box off;

mean_ct_s2d = squeeze(mean(acc_ct_s2d, 1));
figure;
imagesc(time_out, time_out, mean_ct_s2d);
axis xy; hold on;
xline(0, 'LineStyle','--','Color','w','LineWidth',1);
yline(0, 'LineStyle','--','Color','w','LineWidth',1);
colorbar; colormap('jet');
xlabel('Test time — OneDR (ms)');
ylabel('Train time — Stable (ms)');
title('Cross-temporal: Train Stable \rightarrow Test OneDR, RF, FEF');
axis square;

end  % if StableToOneDR

%% ============================================================
%% Cross-context decoder: Train OneDR -> Test Stable
%% Repeated stratified k-fold CV on oneDR (training) context.
%% Z-score from oneDR training fold stats applied to stable test set.
%% ============================================================

if OneDRToStable

% Trial counts for oneDR training context (compute if not already available)
if ~exist('n_hi_oneDR', 'var')
    n_hi_oneDR = zeros(N, 1);
    n_lo_oneDR = zeros(N, 1);
    for n = 1:N
        n_hi_oneDR(n) = size(psth_high_oneDR_smoothed{n}, 1);
        n_lo_oneDR(n) = size(psth_low_oneDR_smoothed{n},  1);
    end
end

% Trial counts for stable test context (compute if not already available)
if ~exist('n_hi_stable', 'var')
    n_hi_stable = zeros(N, 1);
    n_lo_stable = zeros(N, 1);
    for n = 1:N
        n_hi_stable(n) = size(psth_high_stable_smoothed{n}, 1);
        n_lo_stable(n) = size(psth_low_stable_smoothed{n},  1);
    end
end

% Pre-compute fold assignments for oneDR training context in serial.
if ~exist('pf_hi_oneDR', 'var')
    pf_hi_oneDR = cell(n_repeats, N);
    pf_lo_oneDR = cell(n_repeats, N);
    for rep_pre = 1:n_repeats
        for n = 1:N
            pf_hi_oneDR{rep_pre, n} = AssignFolds(n_hi_oneDR(n), k_folds);
            pf_lo_oneDR{rep_pre, n} = AssignFolds(n_lo_oneDR(n), k_folds);
        end
    end
end

acc_d2s    = zeros(n_repeats, T_out);
acc_ct_d2s = zeros(n_repeats, T_out, T_out);

parfor rep = 1:n_repeats
    fold_hi = pf_hi_oneDR(rep, :);   % {1 x N} fold labels for oneDR high trials
    fold_lo = pf_lo_oneDR(rep, :);   % {1 x N} fold labels for oneDR low trials

    rep_acc_sum = zeros(1, T_out);
    rep_ct_sum  = zeros(T_out, T_out);

    for f = 1:k_folds
        % --- Training: OneDR, folds ~= f (resampled) ---
        X_tr_oneDR = zeros(2*n_resample, N, T_out, 'single');
        for n = 1:N
            hi_data   = psth_high_oneDR_smoothed{n};
            lo_data   = psth_low_oneDR_smoothed{n};
            hi_tr_idx = find(fold_hi{n} ~= f);
            lo_tr_idx = find(fold_lo{n} ~= f);
            rs_hi = hi_tr_idx(randi(numel(hi_tr_idx), n_resample, 1));
            rs_lo = lo_tr_idx(randi(numel(lo_tr_idx), n_resample, 1));
            X_tr_oneDR(:, n, :) = single([hi_data(rs_hi, :); lo_data(rs_lo, :)]);
        end

        % --- Test: Stable, all trials (resampled) ---
        X_te_stable = zeros(2*n_test_resample, N, T_out, 'single');
        for n = 1:N
            hi_data = psth_high_stable_smoothed{n};
            lo_data = psth_low_stable_smoothed{n};
            idx_hi  = randi(size(hi_data, 1), n_test_resample, 1);
            idx_lo  = randi(size(lo_data, 1), n_test_resample, 1);
            X_te_stable(:, n, :) = single([hi_data(idx_hi, :); lo_data(idx_lo, :)]);
        end

        y_tr = [ones(n_resample,      1); zeros(n_resample,      1)];
        y_te = [ones(n_test_resample, 1); zeros(n_test_resample, 1)];

        % Z-score from oneDR training fold stats; apply same transform to stable
        mu = mean(X_tr_oneDR, 1);
        sd = std(X_tr_oneDR, 0, 1) + 1e-8;
        X_tr_z = (X_tr_oneDR  - mu) ./ sd;
        X_te_z = (X_te_stable - mu) ./ sd;

        fold_acc = zeros(1, T_out);
        fold_ct  = zeros(T_out, T_out);

        for t_train = 1:T_out
            X_tr = squeeze(X_tr_z(:, :, t_train));
            mdl  = fitclinear(X_tr, y_tr, 'Learner', 'logistic', ...
                'Regularization', 'ridge', 'Lambda', 'auto', 'Solver', 'bfgs');
            fold_acc(t_train) = mean(predict(mdl, squeeze(X_te_z(:, :, t_train))) == y_te);
            for t_test = 1:T_out
                fold_ct(t_train, t_test) = mean(predict(mdl, squeeze(X_te_z(:, :, t_test))) == y_te);
            end
        end

        rep_acc_sum = rep_acc_sum + fold_acc;
        rep_ct_sum  = rep_ct_sum  + fold_ct;
    end

    acc_d2s(rep, :)       = rep_acc_sum / k_folds;
    acc_ct_d2s(rep, :, :) = rep_ct_sum  / k_folds;
    fprintf('Cross D->S Rep %d/%d done\n', rep, n_repeats);
end

%% Permutation test: Train OneDR -> Test Stable (shuffle training labels)
% One random 80/20 split per permutation (mirrors one k=5 fold held-out).
acc_perm_d2s = zeros(n_perm, T_out);

parfor perm = 1:n_perm
    X_tr_oneDR  = zeros(2*n_resample,      N, T_out, 'single');
    X_te_stable = zeros(2*n_test_resample, N, T_out, 'single');

    for n = 1:N
        hi_data   = psth_high_oneDR_smoothed{n};
        lo_data   = psth_low_oneDR_smoothed{n};
        perm_hi   = randperm(n_hi_oneDR(n));
        n_tr_hi   = round(n_hi_oneDR(n) * 0.8);
        hi_tr_idx = perm_hi(1:n_tr_hi);
        perm_lo   = randperm(n_lo_oneDR(n));
        n_tr_lo   = round(n_lo_oneDR(n) * 0.8);
        lo_tr_idx = perm_lo(1:n_tr_lo);
        rs_hi = hi_tr_idx(randi(n_tr_hi, n_resample, 1));
        rs_lo = lo_tr_idx(randi(n_tr_lo, n_resample, 1));
        X_tr_oneDR(:, n, :) = single([hi_data(rs_hi, :); lo_data(rs_lo, :)]);

        hi_s   = psth_high_stable_smoothed{n};
        lo_s   = psth_low_stable_smoothed{n};
        idx_hi = randi(size(hi_s, 1), n_test_resample, 1);
        idx_lo = randi(size(lo_s, 1), n_test_resample, 1);
        X_te_stable(:, n, :) = single([hi_s(idx_hi, :); lo_s(idx_lo, :)]);
    end

    y_shuf = [ones(n_resample,      1); zeros(n_resample,      1)];
    y_shuf = y_shuf(randperm(2*n_resample));
    y_true = [ones(n_test_resample, 1); zeros(n_test_resample, 1)];

    mu = mean(X_tr_oneDR, 1);
    sd = std(X_tr_oneDR, 0, 1) + 1e-8;
    X_tr_z = (X_tr_oneDR  - mu) ./ sd;
    X_te_z = (X_te_stable - mu) ./ sd;

    perm_acc = zeros(1, T_out);
    for t = 1:T_out
        mdl = fitclinear(squeeze(X_tr_z(:,:,t)), y_shuf, 'Learner', 'logistic', ...
            'Regularization', 'ridge', 'Lambda', 'auto', 'Solver', 'bfgs');
        perm_acc(t) = mean(predict(mdl, squeeze(X_te_z(:,:,t))) == y_true);
    end
    acc_perm_d2s(perm, :) = perm_acc;
    if mod(perm, 10) == 0
        fprintf('Cross D->S Perm %d/%d done\n', perm, n_perm);
    end
end

%% Save cross oneDR->stable
SaveFile_D2S = fullfile(DataPath, 'Decoder_Cross_OneDR2Stable_RF_FEF.mat');
save(SaveFile_D2S, 'acc_d2s', 'acc_ct_d2s', 'acc_perm_d2s', 'time_out', ...
    'n_resample', 'n_test_resample', 'n_repeats', 'k_folds', 'n_perm', 'brain_area', '-v7.3');
fprintf('Saved -> %s\n', SaveFile_D2S);

%% Significance test: cross oneDR->stable (cluster correction)
mean_d2s      = mean(acc_d2s,      1);
sem_d2s       = std(acc_d2s,      0, 1) / sqrt(n_repeats);
mean_perm_d2s = mean(acc_perm_d2s, 1);
sem_perm_d2s  = std(acc_perm_d2s, 0, 1) / sqrt(n_perm);

p_val_d2s = zeros(1, T_out);
for t = 1:T_out
    p_val_d2s(t) = mean(acc_perm_d2s(:, t) >= mean_d2s(t));
end

alpha_d2s      = 0.05;
sig_bins_d2s   = p_val_d2s < alpha_d2s;
clusters_d2s   = bwlabel(sig_bins_d2s);
n_clusters_d2s = max(clusters_d2s);

cluster_mass_d2s = zeros(1, n_clusters_d2s);
for c = 1:n_clusters_d2s
    cluster_mass_d2s(c) = sum(1 - p_val_d2s(clusters_d2s == c));
end

null_cluster_mass_d2s = zeros(n_perm, 1);
for perm = 1:n_perm
    p_perm_t = zeros(1, T_out);
    for t = 1:T_out
        p_perm_t(t) = mean(acc_perm_d2s(:, t) >= acc_perm_d2s(perm, t));
    end
    sig_perm   = p_perm_t < alpha_d2s;
    cl_perm    = bwlabel(sig_perm);
    n_cl       = max(cl_perm);
    if n_cl > 0
        masses = zeros(1, n_cl);
        for c = 1:n_cl
            masses(c) = sum(1 - p_perm_t(cl_perm == c));
        end
        null_cluster_mass_d2s(perm) = max(masses);
    end
end

mass_thresh_d2s   = prctile(null_cluster_mass_d2s, 95);
sig_corrected_d2s = false(1, T_out);
for c = 1:n_clusters_d2s
    if cluster_mass_d2s(c) > mass_thresh_d2s
        sig_corrected_d2s(clusters_d2s == c) = true;
    end
end

fprintf('Cross D->S significant bins (cluster-corrected): %d / %d\n', sum(sig_corrected_d2s), T_out);

%% Plot cross oneDR->stable
dt_d2s      = time_out(2) - time_out(1);
sig_idx_d2s = find(sig_corrected_d2s);

figure; hold on;
fill([time_out, fliplr(time_out)], ...
     [mean_perm_d2s + sem_perm_d2s, fliplr(mean_perm_d2s - sem_perm_d2s)], ...
     [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(time_out, mean_perm_d2s, 'k-', 'LineWidth', 1);
fill([time_out, fliplr(time_out)], ...
     [mean_d2s + sem_d2s, fliplr(mean_d2s - sem_d2s)], ...
     [0.8 0.5 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(time_out, mean_d2s, 'Color', [0.7 0.3 0.1], 'LineWidth', 1.5);
for ii = 1:numel(sig_idx_d2s)
    t_c = time_out(sig_idx_d2s(ii));
    fill([t_c - dt_d2s/2, t_c + dt_d2s/2, t_c + dt_d2s/2, t_c - dt_d2s/2], ...
         [0.94, 0.94, 0.96, 0.96], [0.7 0.3 0.1], 'EdgeColor', 'none');
end
yline(0.5, 'LineStyle','--','Color','k');
xline(0,   'LineStyle','--','Color','k');
xlabel('Time (ms)');
ylabel('Decoding accuracy');
title('Cross-context: Train OneDR \rightarrow Test Stable, RF, FEF');
xlim([time_out(1), time_out(end)]);
ylim([0.3 1]);
box off;

mean_ct_d2s = squeeze(mean(acc_ct_d2s, 1));
figure;
imagesc(time_out, time_out, mean_ct_d2s);
axis xy; hold on;
xline(0, 'LineStyle','--','Color','w','LineWidth',1);
yline(0, 'LineStyle','--','Color','w','LineWidth',1);
colorbar; colormap('jet');
xlabel('Test time — Stable (ms)');
ylabel('Train time — OneDR (ms)');
title('Cross-temporal: Train OneDR \rightarrow Test Stable, RF, FEF');
axis square;

end  % if OneDRToStable

