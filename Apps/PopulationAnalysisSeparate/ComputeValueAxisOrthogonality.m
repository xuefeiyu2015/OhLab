% ComputeValueAxisOrthogonality.m
% Standalone script — quantify the angle between value-encoding axes in the
% object (stable) vs. location (oneDR) contexts.
%
% Data source: ObjectLocationPSTH/[AreaLabel]/PopPSTH_RF_[AreaLabel].mat
% Prerequisites: run Batch_ObjLocValue_Organize with OutputFlag=true first.
%
% Value axis at each time bin is defined as:
%   w_obj(t) = mean(high_obj) - mean(low_obj)
%   w_loc(t) = mean(high_loc) - mean(low_loc)
%
% Reports unsigned angle in [0, 90]:
%   theta(t) = acosd( abs( dot(w_obj_norm, w_loc_norm) ) )

%% ── Configuration ────────────────────────────────────────────────────────────
AreaLabel       = 'FEF';    % 'FEF' or 'SC'
n_repeats       = 100;
n_boot          = 5;        % trials drawn with replacement per neuron per condition
window_mode     = 'peak';    % 'peak': auto-detect from grand mean | 'fixed': use window_ms_fixed
window_ms_fixed = [100 200]; % ms — used only when window_mode = 'fixed'
peak_search_range = [0 300]; % ms post-stimulus to search for peak angle (peak mode only)
peak_half_win   = 50;        % ± ms around peak for summary window (peak mode only)
seed            = 'shuffle'; % random seed each run
do_shuffle_null = true;
n_null          = 100;
plot_result     = true;
save_figs       = false;     % set true to export all figures as EPS to fig_save_dir
fig_save_dir    = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/DataHub_YXF/Apps/PopulationAnalysisSeparate/ExportFig';        % destination folder (created if needed); e.g. 'Figures/FEF'
save_path       = '';        % set to a full .mat path to save Results, or leave ''
do_sanity_check = false;
do_regression_analysis = true;   % set false to skip (~10 s for N≈100 neurons)
do_decoder             = false;    % plot within- and cross-context decoding accuracy
do_proj_auc            = true;    % compute value discriminability AUC along each axis

%% ── Paths & load ─────────────────────────────────────────────────────────────
DataRoot = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/ObjectLocationPSTH';
DataFile = fullfile(DataRoot, AreaLabel, sprintf('PopPSTH_RF_%s.mat', AreaLabel));

load(DataFile, ...
    'psth_targeton_raw_stable_RF_high', ...
    'psth_targeton_raw_stable_RF_low', ...
    'psth_targeton_raw_oneDR_RF_high', ...
    'psth_targeton_raw_oneDR_RF_low', ...
    'ModuType_Sel', ...
    'PSTH_Time_TgtOn_stable', ...
    'PSTH_Time_TgtOn_OneDR');

fprintf('Loaded: %s\n', DataFile);
fprintf('  Total neurons: %d\n', numel(psth_targeton_raw_stable_RF_high));
fprintf('  ModuType counts: Fast=%d, Slow=%d, Inhib=%d, NoEffect=%d\n', ...
    sum(ModuType_Sel==1), sum(ModuType_Sel==2), ...
    sum(ModuType_Sel==3), sum(ModuType_Sel==4));

%% ── Neuron group selection ───────────────────────────────────────────────────
% All neurons: sel = true(size(ModuType_Sel));
% Fast only:   sel = ModuType_Sel == 1;
% Others only: sel = ModuType_Sel == 4;
sel = true(size(ModuType_Sel));
%sel = ModuType_Sel == 1;
%sel = ModuType_Sel == 4;
psth_high_value_object_RF   = psth_targeton_raw_stable_RF_high(sel);
psth_low_value_object_RF    = psth_targeton_raw_stable_RF_low(sel);
psth_high_value_location_RF = psth_targeton_raw_oneDR_RF_high(sel);
psth_low_value_location_RF  = psth_targeton_raw_oneDR_RF_low(sel);

%% ── Preprocess: Gaussian smooth (30 ms) + crop to [-500, 300] ms ─────────────
[psth_high_value_object_RF,   time_stable] = PreprocessPSTH( ...
    psth_high_value_object_RF,   PSTH_Time_TgtOn_stable, 'kernel', 'gaussian', 'win_ms', 30);
[psth_low_value_object_RF,    ~          ] = PreprocessPSTH( ...
    psth_low_value_object_RF,    PSTH_Time_TgtOn_stable, 'kernel', 'gaussian', 'win_ms', 30);

[psth_high_value_location_RF, time_oneDR ] = PreprocessPSTH( ...
    psth_high_value_location_RF, PSTH_Time_TgtOn_OneDR,  'kernel', 'gaussian', 'win_ms', 30);
[psth_low_value_location_RF,  ~          ] = PreprocessPSTH( ...
    psth_low_value_location_RF,  PSTH_Time_TgtOn_OneDR,  'kernel', 'gaussian', 'win_ms', 30);

% Crop stable to [-500, 300] ms
t_sel_stable = time_stable >= -500 & time_stable <= 300;
time_axis    = time_stable(t_sel_stable);
for n = 1:numel(psth_high_value_object_RF)
    psth_high_value_object_RF{n} = psth_high_value_object_RF{n}(:, t_sel_stable);
    psth_low_value_object_RF{n}  = psth_low_value_object_RF{n}(:,  t_sel_stable);
end

% Crop oneDR to the same [-500, 300] ms window
t_sel_oneDR = time_oneDR >= -500 & time_oneDR <= 300;
for n = 1:numel(psth_high_value_location_RF)
    psth_high_value_location_RF{n} = psth_high_value_location_RF{n}(:, t_sel_oneDR);
    psth_low_value_location_RF{n}  = psth_low_value_location_RF{n}(:,  t_sel_oneDR);
end

%% ── Validate inputs ──────────────────────────────────────────────────────────
N     = numel(psth_high_value_object_RF);
nTime = numel(time_axis);

assert(numel(psth_low_value_object_RF)    == N, 'Neuron count mismatch: low object.');
assert(numel(psth_high_value_location_RF) == N, 'Neuron count mismatch: high location.');
assert(numel(psth_low_value_location_RF)  == N, 'Neuron count mismatch: low location.');

rng(seed);

fprintf('\n=== ComputeValueAxisOrthogonality ===\n');
fprintf('  Area            : %s\n', AreaLabel);
fprintf('  Neurons         : %d\n', N);
fprintf('  Time bins       : %d\n', nTime);
fprintf('  Repeats         : %d\n', n_repeats);
fprintf('  Bootstrap trials: %d (with replacement)\n', n_boot);
fprintf('  Window mode     : %s\n', window_mode);
%{
%% ── Grand mean axes (all trials, no bootstrap) ───────────────────────────────
% Used to determine the peak timepoint and as the displayed timecourse.
% Theta timecourse uses the average-first approach: at each time bin, average
% the value axes over a ±peak_half_win window, then compute theta — identical
% to the angle_at_window logic used in the bootstrap and null sections.
w_obj_grand = zeros(N, nTime);
w_loc_grand = zeros(N, nTime);
for n = 1:N
    w_obj_grand(n, :) = mean(psth_high_value_object_RF{n}, 1) - mean(psth_low_value_object_RF{n}, 1);
    w_loc_grand(n, :) = mean(psth_high_value_location_RF{n}, 1) - mean(psth_low_value_location_RF{n}, 1);
end



theta_grand = nan(1, nTime);
for ti = 1:nTime
    cwin   = time_axis(ti) + [-peak_half_win, peak_half_win];
    t_cwin = time_axis >= cwin(1) & time_axis <= cwin(2);
    if sum(t_cwin) == 0, continue; end
    theta_grand(ti) = angle_at_window(w_obj_grand, w_loc_grand, t_cwin);
end

%% ── Window selection ─────────────────────────────────────────────────────────
if strcmp(window_mode, 'peak')
    t_search = time_axis >= peak_search_range(1) & time_axis <= peak_search_range(2);
    [~, peak_idx_local] = max(theta_grand(t_search));
    search_times = time_axis(t_search);
    peak_time = search_times(peak_idx_local);
    window_ms = [peak_time - peak_half_win, peak_time + peak_half_win];
    fprintf('  Peak angle time : %.0f ms → summary window [%.0f, %.0f] ms\n', ...
        peak_time, window_ms(1), window_ms(2));
else  % 'fixed'
    window_ms = window_ms_fixed;
    peak_time = NaN;
    fprintf('  Fixed window    : [%.0f, %.0f] ms\n', window_ms(1), window_ms(2));
end
t_win = time_axis >= window_ms(1) & time_axis <= window_ms(2);
%}

%% Direct mean
avgTrial = @(cellArr) cell2mat( ...
    cellfun(@(m) mean(m(:,:), 1), cellArr(:), 'UniformOutput', false) );

meanpsth_high_value_object_RF   = avgTrial(psth_high_value_object_RF);
meanpsth_low_value_object_RF    = avgTrial(psth_low_value_object_RF);
meanpsth_high_value_location_RF = avgTrial(psth_high_value_location_RF);
meanpsth_low_value_location_RF  = avgTrial(psth_low_value_location_RF);


%Project the activities onto each axis



%% ── Bootstrap repeats (angle computed only at summary window) ────────────────
theta_win = nan(n_repeats, 1);
mean_w_loc = nan(n_repeats,N,nTime);
mean_w_obj = nan(n_repeats,N,nTime);
for rep = 1:n_repeats
    [w_obj, w_loc] = build_axes_once( ...
        psth_high_value_object_RF, psth_low_value_object_RF, ...
        psth_high_value_location_RF, psth_low_value_location_RF, ...
         false, n_boot);

    mean_w_loc(rep,:,:) = w_loc;
    mean_w_obj(rep,:,:) = w_obj;

   
end

avg_w_loc = squeeze(mean(mean_w_loc,1));
avg_w_obj = squeeze(mean(mean_w_obj,1));

%Search for the maximum time window
[~,peak_idx_loc] = max(abs(mean(avg_w_loc,1)));
[~,peak_idx_obj] = max(abs(mean(avg_w_obj,1)));

time_interval_loc = [time_axis(peak_idx_loc)-50 time_axis(peak_idx_loc)+50]
time_interval_obj = [time_axis(peak_idx_obj)-50 time_axis(peak_idx_obj)+50]

if contains(window_mode,'peak')
    twin_loc = round([peak_idx_loc-50/10, min(peak_idx_loc+50/10,80)]);
    twin_obj = round([peak_idx_obj-50/10, min(peak_idx_obj+50/10,80)]);
else
    twin_loc = [find(time_axis > 100, 1, 'first'), find(time_axis < 200, 1, 'last')];
    twin_obj = twin_loc;
end
%{
%Population Average
meanpsth_high_value_object_RF   = avgTrial(psth_high_value_object_RF);
meanpsth_low_value_object_RF    = avgTrial(psth_low_value_object_RF);
meanpsth_high_value_location_RF = avgTrial(psth_high_value_location_RF);
meanpsth_low_value_location_RF  = avgTrial(psth_low_value_location_RF);

value_obj = meanpsth_high_value_object_RF-meanpsth_low_value_object_RF ;
location_obj = meanpsth_high_value_location_RF-meanpsth_low_value_location_RF ;

mean_obj = mean(value_obj(:,twin_obj),2);
mean_loc = mean(location_obj(:,twin_loc),2);

norm_obj = norm(mean_obj);
norm_loc = norm(mean_loc);

c  = dot(mean_obj, mean_loc) / (norm_obj * norm_loc);

poptheta = acosd((c));
%}



 theta_grand = angle_at_window(mean(mean_w_loc,1), mean(mean_w_obj,1), twin_loc ,twin_obj);
 theta_bootstrap = angle_at_window(mean_w_loc, mean_w_obj, twin_loc ,twin_obj);
 theta_win_ci   = prctile(theta_bootstrap, [2.5 97.5]);
 theta_mean = theta_grand;   % point estimate: average axes first, then compute angle

%% ── Null model (trial-label shuffle, angle at summary window only) ───────────
null = struct();
if do_shuffle_null
    theta_null_win = nan(n_null, 1);
    w_obj_n = nan(n_null,N,nTime);
     w_loc_n = nan(n_null,N,nTime);

    for i = 1:n_null
        [w_obj_n(i,:,:), w_loc_n(i,:,:)] = build_axes_once( ...
            psth_high_value_object_RF, psth_low_value_object_RF, ...
            psth_high_value_location_RF, psth_low_value_location_RF, ...
            true,  n_boot);
    end

    theta_null_win = angle_at_window(w_loc_n, w_obj_n, twin_loc, twin_obj);
    null.theta_win       = theta_null_win;
    null.theta_win_mean  = nanmean(theta_null_win);
    null.theta_win_ci    = prctile(theta_null_win, [2.5 97.5]);
    null.p_window_orth      = mean(theta_null_win >= theta_grand);
    null.p_window_align     = mean(theta_null_win <= theta_grand);
    null.p_window_two_sided = 2 * min(null.p_window_orth, null.p_window_align);
end

%% ── Null for theta = 0°: same-axis resample ─────────────────────────────────
% Draw two independent bootstrap resamples from the SAME context (obj).
% The angle between them reflects sampling noise when the true axes are
% perfectly parallel — i.e., the smallest theta could be by chance alone.
% p_nonparallel = fraction of same-axis angles >= theta_grand.
fprintf('\n=== Computing theta=0 null (same-axis resampling, n=%d) ===\n', n_null);
%theta_null0 = nan(n_null, 1);
%theta_null_win = nan(n_null, 1);
     w_obj_n10 = nan(n_null,N,nTime);
     w_obj_n20 = nan(n_null,N,nTime);
     w_loc_n10 = nan(n_null,N,nTime);
     w_loc_n20 = nan(n_null,N,nTime);
for i = 1:n_null
    
     

    [w_obj_n10(i,:,:), w_obj_n20(i,:,:)] = build_axes_once( ...
            psth_high_value_object_RF, psth_low_value_object_RF, ...
            psth_high_value_object_RF, psth_low_value_object_RF, ...
            false,  n_boot);

     [w_loc_n10(i,:,:), w_loc_n20(i,:,:)] = build_axes_once( ...
            psth_high_value_location_RF, psth_low_value_location_RF, ...
            psth_high_value_location_RF, psth_low_value_location_RF, ...
            false,  n_boot);

    %{
    for n = 1:N
        [mu_hi_a, mu_lo_a] = bootstrap_and_mean( ...
            psth_high_value_object_RF{n}, psth_low_value_object_RF{n}, n_boot);
        [mu_hi_b, mu_lo_b] = bootstrap_and_mean( ...
            psth_high_value_object_RF{n}, psth_low_value_object_RF{n},n_boot);
        w_a(n, :) = mu_hi_a - mu_lo_a;
        w_b(n, :) = mu_hi_b - mu_lo_b;
    end
    
    a_win = mean(w_a(:, twin_obj(1):twin_obj(2)), 2);
    b_win = mean(w_b(:, twin_obj(1):twin_obj(2)), 2);
    c = dot(a_win, b_win) / (norm(a_win) * norm(b_win));
    theta_null0(i) = acosd(abs(max(-1, min(1, c))));
%}
end
theta_null0_win_obj = angle_at_window(w_obj_n10, w_obj_n20, twin_obj, twin_obj);
theta_null0_win_loc = angle_at_window(w_loc_n10, w_loc_n20, twin_loc, twin_loc);
theta_null0= [theta_null0_win_obj+theta_null0_win_loc]/2;
null.theta_null0    = theta_null0;
null.theta_null0_ci = prctile(theta_null0, [2.5 97.5]);
null.p_nonparallel  = mean(theta_null0 >= theta_grand);
fprintf('  theta=0 null range: %.1f – %.1f deg (95%% CI: %.1f – %.1f)\n', ...
    min(theta_null0), max(theta_null0), null.theta_null0_ci(1), null.theta_null0_ci(2));
fprintf('  p_nonparallel (theta >= observed): %.4f\n', null.p_nonparallel);

%% ── Pack Results ─────────────────────────────────────────────────────────────
Results = struct();
Results.theta_win       = theta_bootstrap;  % per-repeat angle in summary window
Results.theta_win_mean  = theta_grand;
Results.theta_win_ci    = theta_win_ci;
Results.time_axis       = time_axis;
Results.null            = null;
Results.n_neurons       = N;
Results.area            = AreaLabel;
Results.run_timestamp   = datestr(now, 'yyyy-mm-dd HH:MM:SS');

% Data for plots
Results.avg_w_obj       = avg_w_obj;        % [N × T] averaged obj-value axis
Results.avg_w_loc       = avg_w_loc;        % [N × T] averaged loc-value axis
Results.boot_w_obj      = mean_w_obj;       % [n_repeats × N × T] per-repeat obj-value axes
Results.boot_w_loc      = mean_w_loc;       % [n_repeats × N × T] per-repeat loc-value axes
Results.twin_obj        = twin_obj;         % index range for obj peak window
Results.twin_loc        = twin_loc;         % index range for loc peak window
Results.theta_mean      = theta_mean;       % scalar summary angle
Results.psth_high_obj   = meanpsth_high_value_object_RF;    % [N × T]
Results.psth_low_obj    = meanpsth_low_value_object_RF;     % [N × T]
Results.psth_high_loc   = meanpsth_high_value_location_RF;  % [N × T]
Results.psth_low_loc    = meanpsth_low_value_location_RF;   % [N × T]

%% ── Decoder: within & cross-context accuracy over time ───────────────────────
% Use the single peak-window-averaged, normalised axis vectors — the same
% vectors whose angle defines theta. Fixing one axis for all time bins
% guarantees that cos(theta) ≈ 0 → chance cross-decoding, by construction.
% (Using per-bin avg_w_obj(:,t) introduced a different axis at every t,
%  breaking the link to the window-averaged theta.)
decoder_K_fold = 5;
if do_decoder
    % Peak-window axes (mirrors the projection plots in plot_results)
    w_obj_win = mean(avg_w_obj(:, round(twin_obj(1)):round(twin_obj(2))), 2);  % [N x 1]
    w_loc_win = mean(avg_w_loc(:, round(twin_loc(1)):round(twin_loc(2))), 2);
    w_obj_dec = w_obj_win / norm(w_obj_win);   % unit vector — same as used for theta
    w_loc_dec = w_loc_win / norm(w_loc_win);

    fprintf('\n=== Computing decoder accuracy timecourses (K=%d-fold CV) ===\n', decoder_K_fold);
    fprintf('  dot(w_obj_dec, w_loc_dec) = %.4f  (cos theta = %.4f)\n', ...
        dot(w_obj_dec, w_loc_dec), cosd(theta_mean));
    [acc_wo, acc_wl, acc_o2l, acc_l2o, acc_wo_sem, acc_wl_sem, acc_o2l_sem, acc_l2o_sem] = ...
        compute_decoder_accuracy( ...
        psth_high_value_object_RF,  psth_low_value_object_RF, ...
        psth_high_value_location_RF, psth_low_value_location_RF, ...
        w_obj_dec, w_loc_dec, decoder_K_fold);

    Results.acc_within_obj        = acc_wo;        % [1 × T]
    Results.acc_within_loc        = acc_wl;
    Results.acc_cross_obj2loc     = acc_o2l;       % obj axis → loc trials
    Results.acc_cross_loc2obj     = acc_l2o;       % loc axis → obj trials
    Results.acc_within_obj_sem    = acc_wo_sem;    % [1 × T] SEM
    Results.acc_within_loc_sem    = acc_wl_sem;
    Results.acc_cross_obj2loc_sem = acc_o2l_sem;
    Results.acc_cross_loc2obj_sem = acc_l2o_sem;
    fprintf('  Peak within obj:  %.3f | within loc:  %.3f\n', max(acc_wo), max(acc_wl));
    fprintf('  Peak cross o2l:   %.3f | cross l2o:   %.3f\n', max(acc_o2l), max(acc_l2o));
end

%% ── Regression variance explained ────────────────────────────────────────────
if do_regression_analysis
    fprintf('\n=== Running regression analysis (partial R²) ===\n');
    [pR2_value, pR2_context, pR2_interact, pR2_all_neurons] = ...
        run_regression_analysis( ...
            psth_high_value_object_RF,  psth_low_value_object_RF, ...
            psth_high_value_location_RF, psth_low_value_location_RF);

    pR2_sem_value    = squeeze(nanstd(pR2_all_neurons(:,:,1), 0, 1)) / sqrt(N);
    pR2_sem_context  = squeeze(nanstd(pR2_all_neurons(:,:,2), 0, 1)) / sqrt(N);
    pR2_sem_interact = squeeze(nanstd(pR2_all_neurons(:,:,3), 0, 1)) / sqrt(N);

    Results.pR2_value        = pR2_value;           % [1 × T]
    Results.pR2_context      = pR2_context;         % [1 × T]
    Results.pR2_interact     = pR2_interact;         % [1 × T]
    Results.pR2_sem_value    = pR2_sem_value;        % [1 × T]
    Results.pR2_sem_context  = pR2_sem_context;      % [1 × T]
    Results.pR2_sem_interact = pR2_sem_interact;     % [1 × T]
    Results.pR2_all_neurons  = pR2_all_neurons;      % [N × T × 3]

    %% ── Shared and context-dependent value axes ──────────────────────────────
    fprintf('\n=== Computing shared and context-dependent value axes ===\n');
    [w_shared, w_ctxdep, w_obj_win_reg, w_loc_win_reg] = ...
        compute_shared_axes(avg_w_obj, avg_w_loc, twin_obj, twin_loc);

    proj_shared = struct();
    proj_shared.high_obj = meanpsth_high_value_object_RF'   * w_shared;   % [T × 1]
    proj_shared.low_obj  = meanpsth_low_value_object_RF'    * w_shared;
    proj_shared.high_loc = meanpsth_high_value_location_RF' * w_shared;
    proj_shared.low_loc  = meanpsth_low_value_location_RF'  * w_shared;

    proj_ctxdep = struct();
    proj_ctxdep.high_obj = meanpsth_high_value_object_RF'   * w_ctxdep;
    proj_ctxdep.low_obj  = meanpsth_low_value_object_RF'    * w_ctxdep;
    proj_ctxdep.high_loc = meanpsth_high_value_location_RF' * w_ctxdep;
    proj_ctxdep.low_loc  = meanpsth_low_value_location_RF'  * w_ctxdep;

    Results.w_shared    = w_shared;           % [N × 1] normalized
    Results.w_ctxdep    = w_ctxdep;           % [N × 1] normalized
    Results.w_obj_win_r = w_obj_win_reg;      % [N × 1] unnormalized obj window vector
    Results.w_loc_win_r = w_loc_win_reg;      % [N × 1] unnormalized loc window vector
    Results.proj_shared = proj_shared;        % struct of [T × 1]
    Results.proj_ctxdep = proj_ctxdep;        % struct of [T × 1]

    %% ── Orthogonalize shared and context-dep axes (Gram-Schmidt) ─────────────
    % Keep w_shared; project its component out of w_ctxdep, then renormalize.
    w_shared_orth = w_shared;   % already unit length
    w_ctxdep_orth = w_ctxdep - dot(w_ctxdep, w_shared_orth) * w_shared_orth;
    w_ctxdep_orth = w_ctxdep_orth / norm(w_ctxdep_orth);
    fprintf('  dot(shared, ctxdep) before orth: %.4f  after: %.4f\n', ...
        dot(w_shared, w_ctxdep), dot(w_shared_orth, w_ctxdep_orth));

    % Projections of mean PSTHs onto orthogonalized axes → [T × 1] each
    proj_orth_shared = struct();
    proj_orth_shared.high_obj = meanpsth_high_value_object_RF'   * w_shared_orth;
    proj_orth_shared.low_obj  = meanpsth_low_value_object_RF'    * w_shared_orth;
    proj_orth_shared.high_loc = meanpsth_high_value_location_RF' * w_shared_orth;
    proj_orth_shared.low_loc  = meanpsth_low_value_location_RF'  * w_shared_orth;

    proj_orth_ctxdep = struct();
    proj_orth_ctxdep.high_obj = meanpsth_high_value_object_RF'   * w_ctxdep_orth;
    proj_orth_ctxdep.low_obj  = meanpsth_low_value_object_RF'    * w_ctxdep_orth;
    proj_orth_ctxdep.high_loc = meanpsth_high_value_location_RF' * w_ctxdep_orth;
    proj_orth_ctxdep.low_loc  = meanpsth_low_value_location_RF'  * w_ctxdep_orth;

    Results.w_shared_orth      = w_shared_orth;      % [N × 1]
    Results.w_ctxdep_orth      = w_ctxdep_orth;      % [N × 1]
    Results.proj_orth_shared   = proj_orth_shared;   % struct of [T × 1]
    Results.proj_orth_ctxdep   = proj_orth_ctxdep;   % struct of [T × 1]

    %% ── Window-averaged regression ───────────────────────────────────────────
    % Average firing rate over the context-appropriate peak window per trial,
    % then fit the value × context model on those scalar responses.
    fprintf('\n=== Running window-averaged regression (partial R²) ===\n');
    [pR2_win_val, pR2_win_ctx, pR2_win_int, pR2_win_neurons, pval_win_neurons] = ...
        run_window_regression( ...
            psth_high_value_object_RF,  psth_low_value_object_RF, ...
            psth_high_value_location_RF, psth_low_value_location_RF, ...
            twin_obj, twin_loc);

    pR2_win_sem_val = nanstd(pR2_win_neurons(:,1)) / sqrt(sum(~isnan(pR2_win_neurons(:,1))));
    pR2_win_sem_ctx = nanstd(pR2_win_neurons(:,2)) / sqrt(sum(~isnan(pR2_win_neurons(:,2))));
    pR2_win_sem_int = nanstd(pR2_win_neurons(:,3)) / sqrt(sum(~isnan(pR2_win_neurons(:,3))));

    % Population-level: Wilcoxon signed-rank test (median partial R² > 0)
    [~, p_pop_val] = signrank(pR2_win_neurons(:,1), 0, 'tail', 'right');
    [~, p_pop_ctx] = signrank(pR2_win_neurons(:,2), 0, 'tail', 'right');
    [~, p_pop_int] = signrank(pR2_win_neurons(:,3), 0, 'tail', 'right');

    Results.pR2_win_value    = pR2_win_val;          % scalar — population mean
    Results.pR2_win_context  = pR2_win_ctx;
    Results.pR2_win_interact = pR2_win_int;
    Results.pR2_win_sem_value    = pR2_win_sem_val;  % scalar — SEM across neurons
    Results.pR2_win_sem_context  = pR2_win_sem_ctx;
    Results.pR2_win_sem_interact = pR2_win_sem_int;
    Results.pR2_win_neurons      = pR2_win_neurons;  % [N × 3]
    Results.pval_win_neurons     = pval_win_neurons; % [N × 3] per-neuron F-test p-values
    Results.pval_win_pop         = [p_pop_val, p_pop_ctx, p_pop_int]; % [1 × 3] Wilcoxon
    fprintf('  pR²  value=%.4f (p=%.4f)  context=%.4f (p=%.4f)  interaction=%.4f (p=%.4f)\n', ...
        pR2_win_val, p_pop_val, pR2_win_ctx, p_pop_ctx, pR2_win_int, p_pop_int);

    % ── Cell-type classification ─────────────────────────────────────────────
    % Each predictor tested independently at p<0.05; 8 overlap groups.
    [sig_flags, group_id, group_names, group_cols] = classify_cell_types(pval_win_neurons);
    Results.sig_flags  = sig_flags;    % [N × 3] boolean (val / ctx / int)
    Results.cell_type  = group_id;     % [N × 1] integer 1–8
    Results.type_names = group_names;  % {8 × 1} label strings
    Results.type_cols  = group_cols;   % [8 × 3] RGB colours
    fprintf('\n  Cell-type counts:\n');
    for tc = 1:8
        if sum(group_id==tc) > 0
            fprintf('    %s: %d (%.1f%%)\n', group_names{tc}, sum(group_id==tc), ...
                100*sum(group_id==tc)/N);
        end
    end
end

%% ── Value discriminability AUC along projection axes ────────────────────────
% For each time bin, pool w(n)*r_n(trial,t) across all N neurons, then
% compute AUC(high value vs low value). Captures global discriminability of
% value along each fixed axis without requiring simultaneously recorded neurons.
if do_proj_auc
    w_obj_win_auc = mean(avg_w_obj(:, round(twin_obj(1)):round(twin_obj(2))), 2);
    w_loc_win_auc = mean(avg_w_loc(:, round(twin_loc(1)):round(twin_loc(2))), 2);
    w_obj_n_auc   = w_obj_win_auc / norm(w_obj_win_auc);
    w_loc_n_auc   = w_loc_win_auc / norm(w_loc_win_auc);

    fprintf('\n=== Computing value discriminability AUC ===\n');
    Results.auc_obj_on_obj = auc_timecourse(psth_high_value_object_RF,   psth_low_value_object_RF,   w_obj_n_auc, nTime);
    Results.auc_loc_on_obj = auc_timecourse(psth_high_value_location_RF, psth_low_value_location_RF, w_obj_n_auc, nTime);
    Results.auc_obj_on_loc = auc_timecourse(psth_high_value_object_RF,   psth_low_value_object_RF,   w_loc_n_auc, nTime);
    Results.auc_loc_on_loc = auc_timecourse(psth_high_value_location_RF, psth_low_value_location_RF, w_loc_n_auc, nTime);
    fprintf('  Peak AUC  obj-on-obj: %.3f | loc-on-obj: %.3f\n', ...
        max(Results.auc_obj_on_obj), max(Results.auc_loc_on_obj));
    fprintf('  Peak AUC  obj-on-loc: %.3f | loc-on-loc: %.3f\n', ...
        max(Results.auc_obj_on_loc), max(Results.auc_loc_on_loc));

    % Null: shuffle high/low labels within each neuron, recompute AUC with the
    % same fixed axes.  95% CI shows what chance looks like given trial counts.
    fprintf('  Computing AUC null (n=%d shuffles)...\n', n_null);
    [null_oo, null_lo, null_ol, null_ll] = auc_null_timecourses( ...
        psth_high_value_object_RF,   psth_low_value_object_RF, ...
        psth_high_value_location_RF, psth_low_value_location_RF, ...
        w_obj_n_auc, w_loc_n_auc, nTime, n_null);
    Results.auc_null_ci_obj_on_obj = prctile(null_oo, [2.5 97.5], 1);  % [2 × nTime]
    Results.auc_null_ci_loc_on_obj = prctile(null_lo, [2.5 97.5], 1);
    Results.auc_null_ci_obj_on_loc = prctile(null_ol, [2.5 97.5], 1);
    Results.auc_null_ci_loc_on_loc = prctile(null_ll, [2.5 97.5], 1);
    fprintf('  Done.\n');
end

%% ── Sanity check (optional) ──────────────────────────────────────────────────
if do_sanity_check
    Results.sanity = run_sanity_check(nTime, n_repeats);
end

%% ── Plot ─────────────────────────────────────────────────────────────────────
if plot_result
    so = struct('enabled', save_figs && ~isempty(fig_save_dir), 'dir', fig_save_dir);
    plot_results(Results, so);
end

%% ── Save ─────────────────────────────────────────────────────────────────────
if ~isempty(save_path)
    save(save_path, 'Results', '-v7.3');
    fprintf('Saved Results -> %s\n', save_path);
end

fprintf('Window angle: %.2f deg (95%% CI: %.2f, %.2f)\n', ...
    theta_mean, theta_win_ci(1), theta_win_ci(2));
if isfield(null, 'p_window_orth')
    fprintf('Shuffle null p(orth): %.4f | p(two-sided): %.4f\n', ...
        null.p_window_orth, null.p_window_two_sided);
end
fprintf('Theta=0 null p(non-parallel): %.4f\n', null.p_nonparallel);

%% ── Local helper functions ───────────────────────────────────────────────────

function [w_obj, w_loc] = build_axes_once( ...
    obj_hi, obj_lo, loc_hi, loc_lo, do_shuffle, n_boot)
% Compute value-difference vectors [N x T] for one bootstrap repeat.
% Draws n_boot trials with replacement from each condition per neuron.
% do_shuffle=true: shuffle high/low labels within each neuron (null model).

if ~do_shuffle
    N     = numel(obj_hi);
    nTime = size(obj_hi{1}, 2);
    w_obj = nan(N, nTime);
    w_loc = nan(N, nTime);
%{
    r_obj_high = nan(N, n_boot,nTime);
    r_obj_low = nan(N, n_boot,nTime);

    r_loc_high = nan(N, n_boot,nTime);
    r_loc_low = nan(N, n_boot,nTime);
%}
    for n = 1:N
        [mu_hi_o, mu_lo_o] = bootstrap_and_mean(obj_hi{n}, obj_lo{n}, n_boot);
        [mu_hi_l, mu_lo_l] = bootstrap_and_mean(loc_hi{n}, loc_lo{n}, n_boot);

        w_obj(n, :) = mu_hi_o - mu_lo_o;
        w_loc(n, :) = mu_hi_l - mu_lo_l;
        %{

        %Save the resample

        r_obj_high(n,:,:) = resp_hi_o;
        r_obj_low(n,:,:) = resp_lo_o;

        r_loc_high(n,:,:) = resp_hi_l;
        r_loc_low(n,:,:) = resp_lo_l;
        %}




    end



else
    % Null model: resample within each context independently, then shuffle
    % value labels within that context. Obj trials never mix with loc trials.
    N     = numel(obj_hi);
    nTime = size(obj_hi{1}, 2);
    w_obj = nan(N, nTime);
    w_loc = nan(N, nTime);

    for n = 1:N
        [mu_hi_o, mu_lo_o] = bootstrap_and_mean(obj_hi{n}, obj_lo{n}, n_boot, true);
        [mu_hi_l, mu_lo_l] = bootstrap_and_mean(loc_hi{n}, loc_lo{n}, n_boot, true);

        w_obj(n, :) = mu_hi_o - mu_lo_o;
        w_loc(n, :) = mu_hi_l - mu_lo_l;




    end


    
end

end

function [mu_hi, mu_lo] = bootstrap_and_mean(x_hi, x_lo, n_boot, do_shuffle)
% Draw n_boot trials with replacement from each condition within this context.
% do_shuffle=true: permute the value labels after resampling (null model).
% do_shuffle=false (default): keep labels intact (bootstrap estimate).
if nargin < 4, do_shuffle = false; end

n_hi = size(x_hi, 1);
n_lo = size(x_lo, 1);

xh = x_hi(randi(n_hi, n_boot, 1), :);
xl = x_lo(randi(n_lo, n_boot, 1), :);

if do_shuffle
    all_x = [xh; xl];
    perm  = randperm(2 * n_boot);
    xh    = all_x(perm(1:n_boot),       :);
    xl    = all_x(perm(n_boot+1:end),   :);
end

mu_hi = mean(xh, 1);
mu_lo = mean(xl, 1);
end
%{
function [theta, valid] = angle_curve(w_a, w_b)
% Compute unsigned angle [0, 90] between w_a and w_b at every time bin.
% w_a/w_b: [N x T]
T     = size(w_a, 2);
theta = nan(1, T);
valid = false(1, T);

eps_n = 1e-12;
for t = 1:T
    a  = w_a(:, t);
    b  = w_b(:, t);
    na = norm(a);
    nb = norm(b);
    if na <= eps_n || nb <= eps_n
        continue;
    end
    c        = dot(a, b) / (na * nb);
    c        = max(-1, min(1, c));
    theta(t) = acosd(abs(c));
    valid(t) = true;
end
end
%}

  
function th = angle_at_window(w_a, w_b, t_a,t_b)

% Average w_a and w_b across the window time bins (collapsing to [N x 1]),
% then compute the angle between the two averaged axis vectors.

[NRep, N,nTime] = size(w_a);
for i = 1: NRep

    a = squeeze(mean(w_a(i,:, t_a(1):t_a(2)), 3));
    b = squeeze(mean(w_b(i,:, t_b(1):t_b(2)), 3));
    na = norm(a);
    nb = norm(b);
    eps_n = 1e-12;
    if na <= eps_n || nb <= eps_n
        th = NaN;
        return;
    end

    c  = dot(a, b) / (na * nb);
    c  = max(-1, min(1, c));
    th(i) = acosd(abs(c));



end



end

function plot_results(R, save_opts)
if nargin < 2, save_opts = struct('enabled', false, 'dir', ''); end
% ── Plot 1: Vector angle diagram ─────────────────────────────────────────────
% Symbolic 2D unit-vector plot showing theta_mean between the two value axes,
% with a grey arc marking the null distribution range (min → max).

theta_obs = R.theta_mean;

hist_scale = 0.22;   % max bar extent (in unit-circle units)
n_hist     = 18;     % number of histogram bins

fig1 = figure('Name', sprintf('Value Axis Angle — %s', R.area), ...
    'Color', 'w', 'Position', [100 200 500 500]);
ax1 = axes(fig1);
hold(ax1, 'on');
axis(ax1, 'equal');

% ── Parallel null histogram (theta_null0) — along horizontal axis (0°) ───────
% Bins in degree-space; x position = 1 - theta/90 to match unit-circle geometry:
% 0° → x=1 (obj-value axis tip), 90° → x=0 (vertical reference line).
% Bars extend downward from y = 0.
if isstruct(R.null) && isfield(R.null, 'theta_null0') && ~isempty(R.null.theta_null0)
    deg_edges  = linspace(0, 90, n_hist + 1);
    pos_edges  = 1 - deg_edges / 90;        % 0°→1, 90°→0 (matches unit-circle x-axis)
    [cnts, ~]  = histcounts(R.null.theta_null0, deg_edges);
    cnts_norm  = cnts / max(cnts + eps) * hist_scale;
    col_par    = [0.53 0.81 0.98];
    for b = 1:n_hist
        if cnts(b) == 0, continue; end
        xl = min(pos_edges(b), pos_edges(b+1));
        xr = max(pos_edges(b), pos_edges(b+1));
        patch(ax1, [xl xr xr xl], [0 0 -cnts_norm(b) -cnts_norm(b)], col_par, ...
            'EdgeColor', 'w', 'LineWidth', 0.5, 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
    end
    % Dotted line at observed angle
    xobs = 1 - theta_obs / 90;
    plot(ax1, [xobs xobs], [0 -(hist_scale + 0.03)], ...
        'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    text(ax1, 0.5, -(hist_scale + 0.07), 'parallel null (theta=0)', ...
        'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', [0.10 0.45 0.70]);
end

% ── Orthogonal null histogram (theta_null_win) — along vertical axis (90°) ───
% Bins in degree-space; y position = theta/90 (linear: 0°→0, 90°→1).
% Bars extend leftward from x = 0.
if isstruct(R.null) && isfield(R.null, 'theta_win') && ~isempty(R.null.theta_win)
    deg_edges  = linspace(0, 90, n_hist + 1);
    pos_edges  = deg_edges / 90;             % linear map to [0, 1]
    [cnts, ~]  = histcounts(R.null.theta_win, deg_edges);
    cnts_norm  = cnts / max(cnts + eps) * hist_scale;
    col_orth   = [0.75 0.75 0.75];
    for b = 1:n_hist
        if cnts(b) == 0, continue; end
        patch(ax1, [-cnts_norm(b) 0 0 -cnts_norm(b)], ...
            [pos_edges(b) pos_edges(b) pos_edges(b+1) pos_edges(b+1)], col_orth, ...
            'EdgeColor', 'w', 'LineWidth', 0.5, 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
    end
    % Dotted line at observed angle
    yobs = theta_obs / 90;
    plot(ax1, [-(hist_scale + 0.03) 0], [yobs yobs], ...
        'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    text(ax1, -(hist_scale + 0.07), 0.5, 'shuffle null', ...
        'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4], 'Rotation', 90);
end

% Reference line: 90°
plot(ax1, [0 cos(pi/2)], [0 sin(pi/2)], '--', 'Color', [0.3 0.3 0.3], ...
    'LineWidth', 1.2);
text(ax1, cos(pi/2)+0.05, sin(pi/2), '90°', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

% Reference line: 0°
plot(ax1, [0 1], [0 0], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
text(ax1, 1.05, 0, '0°', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

% Obj-value axis: unit vector at 0° (blue)
quiver(ax1, 0, 0, 1, 0, 0, 'Color', [0.12 0.47 0.71], 'LineWidth', 2.5, ...
    'MaxHeadSize', 0.4);
text(ax1, 1.08, -0.06, 'Obj-value axis', 'FontSize', 9, ...
    'Color', [0.12 0.47 0.71], 'HorizontalAlignment', 'left');

% Loc-value axis: unit vector at theta_mean° (red)
theta_rad = theta_obs * pi / 180;
quiver(ax1, 0, 0, cos(theta_rad), sin(theta_rad), 0, ...
    'Color', [0.84 0.15 0.16], 'LineWidth', 2.5, 'MaxHeadSize', 0.4);
text(ax1, cos(theta_rad)*1.08, sin(theta_rad)*1.08, 'Loc-value axis', ...
    'FontSize', 9, 'Color', [0.84 0.15 0.16], 'HorizontalAlignment', 'left');

% Angle annotation arc
arc_r2  = 0.3;
phi_obs = linspace(0, theta_obs, 60) * pi / 180;
plot(ax1, arc_r2*cos(phi_obs), arc_r2*sin(phi_obs), 'k-', 'LineWidth', 1);
phi_label = (theta_obs / 2) * pi / 180;
text(ax1, (arc_r2+0.08)*cos(phi_label), (arc_r2+0.08)*sin(phi_label), ...
    sprintf('%.1f°', theta_obs), 'FontSize', 9, 'HorizontalAlignment', 'center');

% Stats annotation
stats_txt = sprintf('θ = %.1f° (CI %.1f–%.1f°)', ...
    theta_obs, R.theta_win_ci(1), R.theta_win_ci(2));
if isstruct(R.null)
    if isfield(R.null, 'p_window_two_sided')
        stats_txt = [stats_txt, sprintf('\np_{shuffle,two-sided} = %.3f', ...
            R.null.p_window_two_sided)];
    end
    if isfield(R.null, 'p_nonparallel')
        stats_txt = [stats_txt, sprintf('\np_{non-parallel} = %.3f', ...
            R.null.p_nonparallel)];
    end
end
text(ax1, -0.05, -(hist_scale + 0.12), stats_txt, 'FontSize', 8, 'VerticalAlignment', 'top');

xlim(ax1, [-(hist_scale + 0.15) 1.4]);
ylim(ax1, [-(hist_scale + 0.25) 1.3]);
axis(ax1, 'off');
title(ax1, sprintf('Value-axis angle — %s', R.area), 'FontSize', 11);
save_figure(fig1, save_opts, sprintf('%s_angle_diagram', R.area));

% ── Plot 2: Population projection onto each value axis ────────────────────────
% Normalize avg_w_obj and avg_w_loc at their respective peak windows.
ta        = R.time_axis;
w_obj_win = mean(R.avg_w_obj(:, R.twin_obj(1):R.twin_obj(2)), 2);  % [N × 1]
w_loc_win = mean(R.avg_w_loc(:, R.twin_loc(1):R.twin_loc(2)), 2);  % [N × 1]

w_obj_norm = w_obj_win / norm(w_obj_win);
w_loc_norm = w_loc_win / norm(w_loc_win);

col_high_obj = [0.84 0.15 0.16];   % red    — obj high value
col_low_obj  = [0.12 0.47 0.71];   % blue   — obj low value
col_high_loc = [0.50 0.15 0.70];   % purple — loc high value
col_low_loc  = [0.17 0.63 0.17];   % green  — loc low value

% Bootstrap 95% CIs: for each repeat, re-derive normalized axis and project
% mean PSTHs → distribution across repeats captures axis-direction uncertainty
if isfield(R, 'boot_w_obj') && isfield(R, 'boot_w_loc')
    n_rep = size(R.boot_w_obj, 1);
    nT    = length(ta);
    psths = {R.psth_high_obj, R.psth_low_obj, R.psth_high_loc, R.psth_low_loc};
    boot_on_obj = nan(n_rep, nT, 4);
    boot_on_loc = nan(n_rep, nT, 4);
    for r = 1:n_rep
        w_o = squeeze(R.boot_w_obj(r,:,:));
        w_l = squeeze(R.boot_w_loc(r,:,:));
        wo_w = mean(w_o(:, round(R.twin_obj(1)):round(R.twin_obj(2))), 2);
        wl_w = mean(w_l(:, round(R.twin_loc(1)):round(R.twin_loc(2))), 2);
        wo_n = wo_w / norm(wo_w);
        wl_n = wl_w / norm(wl_w);
        for c = 1:4
            boot_on_obj(r,:,c) = (psths{c}' * wo_n)';
            boot_on_loc(r,:,c) = (psths{c}' * wl_n)';
        end
    end
    ci_obj = prctile(boot_on_obj, [2.5 97.5], 1);  % [2 × nT × 4]
    ci_loc = prctile(boot_on_loc, [2.5 97.5], 1);

    mean_on_obj = squeeze(mean(boot_on_obj, 1));    % [nT × 4]
    mean_on_loc = squeeze(mean(boot_on_loc, 1));
    proj_high_obj_on_obj = mean_on_obj(:,1);
    proj_low_obj_on_obj  = mean_on_obj(:,2);
    proj_high_loc_on_obj = mean_on_obj(:,3);
    proj_low_loc_on_obj  = mean_on_obj(:,4);
    proj_high_obj_on_loc = mean_on_loc(:,1);
    proj_low_obj_on_loc  = mean_on_loc(:,2);
    proj_high_loc_on_loc = mean_on_loc(:,3);
    proj_low_loc_on_loc  = mean_on_loc(:,4);
    has_ci = true;
else
    proj_high_obj_on_obj = R.psth_high_obj' * w_obj_norm;
    proj_low_obj_on_obj  = R.psth_low_obj'  * w_obj_norm;
    proj_high_loc_on_obj = R.psth_high_loc' * w_obj_norm;
    proj_low_loc_on_obj  = R.psth_low_loc'  * w_obj_norm;
    proj_high_obj_on_loc = R.psth_high_obj' * w_loc_norm;
    proj_low_obj_on_loc  = R.psth_low_obj'  * w_loc_norm;
    proj_high_loc_on_loc = R.psth_high_loc' * w_loc_norm;
    proj_low_loc_on_loc  = R.psth_low_loc'  * w_loc_norm;
    has_ci = false;
end

fig2 = figure('Name', sprintf('Population Projection — %s', R.area), ...
    'Color', 'w', 'Position', [550 200 820 380]);

% Panel A: projection onto obj-value axis
ax2a = subplot(1, 2, 1, 'Parent', fig2);
hold(ax2a, 'on');
if has_ci
    fill_ci_lo_hi(ax2a, ta, ci_obj(1,:,1), ci_obj(2,:,1), col_high_obj);
    fill_ci_lo_hi(ax2a, ta, ci_obj(1,:,2), ci_obj(2,:,2), col_low_obj);
    fill_ci_lo_hi(ax2a, ta, ci_obj(1,:,3), ci_obj(2,:,3), col_high_loc);
    fill_ci_lo_hi(ax2a, ta, ci_obj(1,:,4), ci_obj(2,:,4), col_low_loc);
end
plot(ax2a, ta, proj_high_obj_on_obj, '-', 'Color', col_high_obj, 'LineWidth', 2, 'DisplayName', 'High obj');
plot(ax2a, ta, proj_low_obj_on_obj,  '-', 'Color', col_low_obj,  'LineWidth', 2, 'DisplayName', 'Low obj');
plot(ax2a, ta, proj_high_loc_on_obj, '-', 'Color', col_high_loc, 'LineWidth', 2, 'DisplayName', 'High loc');
plot(ax2a, ta, proj_low_loc_on_obj,  '-', 'Color', col_low_loc,  'LineWidth', 2, 'DisplayName', 'Low loc');
xline(0, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);
xlabel(ax2a, 'Time (ms)');
ylabel(ax2a, 'Projection (a.u.)');
title(ax2a, 'Projection onto obj-value axis');
legend(ax2a, 'Location', 'best', 'FontSize', 8);
box(ax2a, 'off');
xlim([-500 300]);

% Panel B: projection onto loc-value axis
ax2b = subplot(1, 2, 2, 'Parent', fig2);
hold(ax2b, 'on');
if has_ci
    fill_ci_lo_hi(ax2b, ta, ci_loc(1,:,1), ci_loc(2,:,1), col_high_obj);
    fill_ci_lo_hi(ax2b, ta, ci_loc(1,:,2), ci_loc(2,:,2), col_low_obj);
    fill_ci_lo_hi(ax2b, ta, ci_loc(1,:,3), ci_loc(2,:,3), col_high_loc);
    fill_ci_lo_hi(ax2b, ta, ci_loc(1,:,4), ci_loc(2,:,4), col_low_loc);
end
plot(ax2b, ta, proj_high_obj_on_loc, '-', 'Color', col_high_obj, 'LineWidth', 2, 'DisplayName', 'High obj');
plot(ax2b, ta, proj_low_obj_on_loc,  '-', 'Color', col_low_obj,  'LineWidth', 2, 'DisplayName', 'Low obj');
plot(ax2b, ta, proj_high_loc_on_loc, '-', 'Color', col_high_loc, 'LineWidth', 2, 'DisplayName', 'High loc');
plot(ax2b, ta, proj_low_loc_on_loc,  '-', 'Color', col_low_loc,  'LineWidth', 2, 'DisplayName', 'Low loc');
xline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1);
xlabel(ax2b, 'Time (ms)');
ylabel(ax2b, 'Projection (a.u.)');
title(ax2b, 'Projection onto loc-value axis');
legend(ax2b, 'Location', 'best', 'FontSize', 8);
box(ax2b, 'off');
xlim([-500 300]);

sgtitle(fig2, sprintf('Population responses projected onto value axes — %s', R.area), ...
    'FontSize', 12);
save_figure(fig2, save_opts, sprintf('%s_population_projection', R.area));

% ── Value discriminability AUC ────────────────────────────────────────────────
if isfield(R, 'auc_obj_on_obj')
    plot_proj_auc(R, save_opts);
end

% ── Decoding accuracy (within & cross-context) ────────────────────────────────
if isfield(R, 'acc_within_obj')
    plot_decoding_accuracy(R, save_opts);
end

% ── New plots if regression analysis was run ──────────────────────────────────
if isfield(R, 'pR2_value')
    plot_regression_variance(R, save_opts);
end
if isfield(R, 'pR2_win_value')
    plot_window_regression(R, save_opts);
end
if isfield(R, 'cell_type')
    plot_cell_type_analysis(R, save_opts);
end
if isfield(R, 'proj_shared')
    plot_shared_ctxdep_axes(R, save_opts);
end
if isfield(R, 'proj_orth_shared')
    plot_2d_trajectory(R, save_opts);
end
end

function S = run_sanity_check(nTime, nRep)
% Quick synthetic checks: aligned vs orthogonal known axes.
N = 60;
rng(1234);
base = randn(N, nTime);

% Aligned
w1    = base;
w2    = base + 0.05*randn(N, nTime);
th_al = nan(nRep, nTime);
for i = 1:nRep
    idx = randi(N, N, 1);
    th_al(i, :) = angle_curve(w1(idx,:), w2(idx,:));
end

% Orthogonal by construction (Gram-Schmidt per time bin)
w3 = randn(N, nTime);
w4 = randn(N, nTime);
for t = 1:nTime
    a = w3(:, t);
    b = w4(:, t);
    b = b - (dot(a,b) / (dot(a,a)+1e-12)) * a;
    w4(:, t) = b;
end
th_or = nan(nRep, nTime);
for i = 1:nRep
    idx = randi(N, N, 1);
    th_or(i, :) = angle_curve(w3(idx,:), w4(idx,:));
end

S = struct();
S.aligned_mean    = nanmean(th_al, 'all');
S.orthogonal_mean = nanmean(th_or, 'all');
S.aligned_ci      = prctile(mean(th_al, 2, 'omitnan'), [2.5 97.5]);
S.orthogonal_ci   = prctile(mean(th_or, 2, 'omitnan'), [2.5 97.5]);

fprintf('Sanity aligned mean angle:    %.2f deg\n', S.aligned_mean);
fprintf('Sanity orthogonal mean angle: %.2f deg\n', S.orthogonal_mean);
end

% ═════════════════════════════════════════════════════════════════════════════
%  run_regression_analysis
% ═════════════════════════════════════════════════════════════════════════════
function [pR2_value, pR2_context, pR2_interact, pR2_all_neurons] = ...
    run_regression_analysis(psth_hi_obj, psth_lo_obj, psth_hi_loc, psth_lo_loc)
% run_regression_analysis  Partial R² for value, context, and value×context
%   interaction at every time bin, using trial-level data from all 4 conditions.
%
% Inputs (cell arrays [neurons × 1], each cell [n_trials × T]):
%   psth_hi_obj, psth_lo_obj  — object context, high/low value
%   psth_hi_loc, psth_lo_loc  — location context, high/low value
%
% Outputs:
%   pR2_value    [1 × T] — population-mean partial R² for value main effect
%   pR2_context  [1 × T] — population-mean partial R² for context main effect
%   pR2_interact [1 × T] — population-mean partial R² for value×context interaction
%   pR2_all_neurons [N × T × 3] — per-neuron partial R² (dim3: value/context/interaction)
%
% Model: y = β₀ + β_V·V + β_C·C + β_VC·(V×C) + ε
%   V = ±1 (high=+1, low=−1),  C = ±1 (object=+1, location=−1)
%
% Partial R² for predictor p = (SS_resid_reduced_p − SS_resid_full) / SS_total

N = numel(psth_hi_obj);
T = size(psth_hi_obj{1}, 2);
pR2_neurons = nan(N, T, 3);

for n = 1:N
    n_ho = size(psth_hi_obj{n}, 1);
    n_lo = size(psth_lo_obj{n}, 1);
    n_hl = size(psth_hi_loc{n}, 1);
    n_ll = size(psth_lo_loc{n}, 1);
    n_total = n_ho + n_lo + n_hl + n_ll;

    % Design matrix columns: [intercept, value, context, value×context]
    X = [ones(n_total, 1), ...
         [+ones(n_ho,1); -ones(n_lo,1); +ones(n_hl,1); -ones(n_ll,1)], ...
         [+ones(n_ho,1); +ones(n_lo,1); -ones(n_hl,1); -ones(n_ll,1)], ...
         [+ones(n_ho,1); -ones(n_lo,1); -ones(n_hl,1); +ones(n_ll,1)]];

    for t = 1:T
        y = [psth_hi_obj{n}(:,t); psth_lo_obj{n}(:,t); ...
             psth_hi_loc{n}(:,t); psth_lo_loc{n}(:,t)];
        SS_total = sum((y - mean(y)).^2);
        if SS_total < 1e-12, continue; end   % silent neuron at this bin

        SS_full = sum((y - X*(X\y)).^2);

        % Partial R² for value (drop col 2)
        Xr = X(:,[1 3 4]);
        pR2_neurons(n, t, 1) = (sum((y - Xr*(Xr\y)).^2) - SS_full) / SS_total;

        % Partial R² for context (drop col 3)
        Xr = X(:,[1 2 4]);
        pR2_neurons(n, t, 2) = (sum((y - Xr*(Xr\y)).^2) - SS_full) / SS_total;

        % Partial R² for interaction (drop col 4)
        Xr = X(:,[1 2 3]);
        pR2_neurons(n, t, 3) = (sum((y - Xr*(Xr\y)).^2) - SS_full) / SS_total;
    end

    if mod(n, 25) == 0
        fprintf('  Regression: neuron %d / %d\n', n, N);
    end
end

pR2_value    = squeeze(nanmean(pR2_neurons(:,:,1), 1));   % [1 × T]
pR2_context  = squeeze(nanmean(pR2_neurons(:,:,2), 1));
pR2_interact = squeeze(nanmean(pR2_neurons(:,:,3), 1));
pR2_all_neurons = pR2_neurons;
end

% ═════════════════════════════════════════════════════════════════════════════
%  run_window_regression
% ═════════════════════════════════════════════════════════════════════════════
function [pR2_val, pR2_ctx, pR2_int, pR2_neurons, pval_neurons] = run_window_regression( ...
    psth_hi_obj, psth_lo_obj, psth_hi_loc, psth_lo_loc, twin_obj, twin_loc)
% run_window_regression  Partial R² + per-neuron F-test p-values from
% window-averaged firing rates.
%
% For each neuron, averages firing rate over the context-appropriate peak
% window (twin_obj for object trials, twin_loc for location trials), then
% fits the value × context model on those scalar responses.
%
% Per-neuron p-values use a partial F-test:
%   F = (SS_drop / 1) / (SS_full_resid / df_resid)
%   p = 1 - fcdf(F, 1, df_resid)   where df_resid = n_total - 4
%
% Inputs:
%   psth_hi/lo_obj/loc  {N x 1} cell, each [n_trials x T]
%   twin_obj / twin_loc  [1 x 2] integer index range for each peak window
%
% Outputs:
%   pR2_val / pR2_ctx / pR2_int   scalar — population-mean partial R²
%   pR2_neurons                   [N x 3] — per-neuron partial R²
%   pval_neurons                  [N x 3] — per-neuron F-test p-values

N = numel(psth_hi_obj);
pR2_neurons  = nan(N, 3);
pval_neurons = nan(N, 3);

for n = 1:N
    % Average over context-appropriate peak window → scalar per trial
    y_ho = mean(psth_hi_obj{n}(:, twin_obj(1):twin_obj(2)), 2);
    y_lo = mean(psth_lo_obj{n}(:, twin_obj(1):twin_obj(2)), 2);
    y_hl = mean(psth_hi_loc{n}(:, twin_loc(1):twin_loc(2)), 2);
    y_ll = mean(psth_lo_loc{n}(:, twin_loc(1):twin_loc(2)), 2);

    y = [y_ho; y_lo; y_hl; y_ll];
    n_ho = numel(y_ho);  n_lo = numel(y_lo);
    n_hl = numel(y_hl);  n_ll = numel(y_ll);
    n_total = n_ho + n_lo + n_hl + n_ll;
    df_resid = n_total - 4;   % full model has 4 params

    % Design matrix: [intercept, value(±1), context(±1), V×C]
    X = [ones(n_total, 1), ...
         [+ones(n_ho,1); -ones(n_lo,1); +ones(n_hl,1); -ones(n_ll,1)], ...
         [+ones(n_ho,1); +ones(n_lo,1); -ones(n_hl,1); -ones(n_ll,1)], ...
         [+ones(n_ho,1); -ones(n_lo,1); -ones(n_hl,1); +ones(n_ll,1)]];

    SS_total = sum((y - mean(y)).^2);
    if SS_total < 1e-12, continue; end

    SS_full = sum((y - X*(X\y)).^2);
    MSE     = SS_full / df_resid;   % mean squared error of full model

    for k = 1:3
        Xr        = X(:, [1:4] ~= (k+1));
        SS_drop   = sum((y - Xr*(Xr\y)).^2) - SS_full;   % extra SS explained
        pR2_neurons(n, k)  = SS_drop / SS_total;
        F                  = (SS_drop / 1) / MSE;         % partial F, df_num=1
        pval_neurons(n, k) = 1 - fcdf(max(F, 0), 1, df_resid);
    end
end

pR2_val = nanmean(pR2_neurons(:,1));
pR2_ctx = nanmean(pR2_neurons(:,2));
pR2_int = nanmean(pR2_neurons(:,3));
end

% ═════════════════════════════════════════════════════════════════════════════
%  plot_window_regression
% ═════════════════════════════════════════════════════════════════════════════
function plot_window_regression(R, save_opts)
if nargin < 2, save_opts = struct('enabled', false, 'dir', ''); end
% plot_window_regression  Bar chart + 3-D scatter of window-averaged partial R².
%
% Figure 1: bar chart — population mean ± SEM per predictor, jittered neuron
%   dots (filled = significant at p<0.05 by per-neuron F-test, hollow = n.s.).
%   Wilcoxon population p-value annotated above each bar.
% Figure 2: 3-D scatter — each dot = one neuron, axes = partial R² for
%   value (x), context (y), value×context (z). Colour = dominant predictor.

D  = R.pR2_win_neurons;    % [N × 3]
valid = ~any(isnan(D), 2);
Dv = D(valid, :);

vals   = [R.pR2_win_value, R.pR2_win_context, R.pR2_win_interact];
sems   = [R.pR2_win_sem_value, R.pR2_win_sem_context, R.pR2_win_sem_interact];
p_pop  = R.pval_win_pop;   % [1 × 3] Wilcoxon p-values

col_val = [0.12 0.47 0.71];
col_ctx = [0.84 0.15 0.16];
col_int = [0.58 0.10 0.70];
cols    = [col_val; col_ctx; col_int];
labels  = {'Value', 'Context', 'Value \times Context'};

% ── Figure 1: Bar chart ───────────────────────────────────────────────────────
fig1 = figure('Name', sprintf('Window Regression — %s', R.area), ...
    'Color', 'w', 'Position', [100 100 460 440]);
ax1 = axes(fig1);
hold(ax1, 'on');

for k = 1:3
    bar(ax1, k, vals(k), 0.55, 'FaceColor', cols(k,:), 'EdgeColor', 'none', 'FaceAlpha', 0.75);
    errorbar(ax1, k, vals(k), sems(k), 'k', 'LineWidth', 1.5, 'CapSize', 6);

    % Jittered neuron dots: filled = significant (p<0.05), hollow = n.s.
    y_n = D(valid, k);
    if isfield(R, 'pval_win_neurons')
        p_n = R.pval_win_neurons(valid, k);
    else
        p_n = ones(size(y_n));
    end
    sig = p_n < 0.05;
    jitter = (rand(numel(y_n), 1) - 0.5) * 0.25;
    % Significant neurons: filled
    if any(sig)
        scatter(ax1, k + jitter(sig),  y_n(sig),  22, cols(k,:)*0.7, 'filled', ...
            'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none');
    end
    % Non-significant neurons: hollow
    if any(~sig)
        scatter(ax1, k + jitter(~sig), y_n(~sig), 22, cols(k,:)*0.6, ...
            'MarkerEdgeColor', cols(k,:)*0.6, 'MarkerFaceColor', 'none', 'LineWidth', 0.8);
    end

    % Wilcoxon p-value annotation above the bar
    y_top = vals(k) + sems(k);
    p_str = pval_label(p_pop(k));
    text(ax1, k, y_top * 1.08 + 0.003, p_str, ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

yline(0, 'LineStyle','-','Color','k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(ax1, 'XTick', 1:3, 'XTickLabel', labels, 'FontSize', 10);
ylabel(ax1, 'Partial R²');
title(ax1, sprintf('Window-averaged regression — %s  (N=%d)', R.area, R.n_neurons));
xlim(ax1, [0.4 3.6]);
box(ax1, 'off');

% ── Figure 2: 3-D scatter ─────────────────────────────────────────────────────
% Colour each neuron by its dominant predictor.
[~, dom] = max(Dv, [], 2);
pt_col = cols(dom, :);

fig2 = figure('Name', sprintf('Regression 3D Scatter — %s', R.area), ...
    'Color', 'w', 'Position', [540 100 520 480]);
ax2 = axes(fig2);
scatter3(ax2, Dv(:,1), Dv(:,2), Dv(:,3), 40, pt_col, 'filled', ...
    'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', 'none');
hold(ax2, 'on');

% Population mean marker
scatter3(ax2, vals(1), vals(2), vals(3), 120, 'k', 'd', 'filled');

% Zero-plane reference lines through the mean
xl = xlim(ax2);  yl = ylim(ax2);  zl = zlim(ax2);
plot3(ax2, [0 0], yl, [0 0], 'k:', 'LineWidth', 0.8);
plot3(ax2, xl, [0 0], [0 0], 'k:', 'LineWidth', 0.8);
plot3(ax2, [0 0], [0 0], zl, 'k:', 'LineWidth', 0.8);

xlabel(ax2, 'Partial R²  (Value)',   'FontSize', 10);
ylabel(ax2, 'Partial R²  (Context)', 'FontSize', 10);
zlabel(ax2, 'Partial R²  (V×C)',     'FontSize', 10);
title(ax2, sprintf('Per-neuron partial R² — %s  (N=%d)', R.area, sum(valid)), 'FontSize', 11);

% Manual legend patches
patch(ax2, nan, nan, nan, col_val, 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'DisplayName', 'Value dominant');
patch(ax2, nan, nan, nan, col_ctx, 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'DisplayName', 'Context dominant');
patch(ax2, nan, nan, nan, col_int, 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'DisplayName', 'V×C dominant');
scatter3(ax2, nan, nan, nan, 120, 'k', 'd', 'filled', 'DisplayName', 'Population mean');
legend(ax2, 'Location', 'best', 'FontSize', 8);

grid(ax2, 'on');
box(ax2, 'off');
view(ax2, 35, 25);
save_figure(fig1, save_opts, sprintf('%s_window_regression_bars', R.area));
save_figure(fig2, save_opts, sprintf('%s_window_regression_3d',   R.area));
end

function s = pval_label(p)
% pval_label  Convert p-value to significance stars string.
if p < 0.001
    s = '***';
elseif p < 0.01
    s = '**';
elseif p < 0.05
    s = '*';
else
    s = 'n.s.';
end
end

% ═════════════════════════════════════════════════════════════════════════════
%  classify_cell_types
% ═════════════════════════════════════════════════════════════════════════════
function [sig_flags, group_id, group_names, group_cols] = classify_cell_types(pval_neurons, alpha)
% classify_cell_types  Independent p<0.05 test per predictor; 8 overlap groups.
%
% Each neuron can be significant for any combination of the three predictors.
% group_id encodes the combination as a binary number (val=bit0, ctx=bit1, int=bit2):
%   1 = none  2 = V only  3 = C only  4 = V+C
%   5 = Int only  6 = V+Int  7 = C+Int  8 = V+C+Int
%
% Inputs:
%   pval_neurons  [N × 3]  F-test p-values (cols: value / context / interact)
%   alpha         significance threshold (default 0.05)
%
% Outputs:
%   sig_flags   [N × 3]  logical — which predictors are significant
%   group_id    [N × 1]  integer 1–8
%   group_names {8 × 1}  readable labels
%   group_cols  [8 × 3]  RGB colours per group

if nargin < 2, alpha = 0.05; end

sig_val = pval_neurons(:,1) < alpha;
sig_ctx = pval_neurons(:,2) < alpha;
sig_int = pval_neurons(:,3) < alpha;
sig_flags = [sig_val, sig_ctx, sig_int];

group_id = 1 + sig_val*1 + sig_ctx*2 + sig_int*4;   % 1–8

group_names = { ...
    'None', ...                  % 1: 000
    'Value only', ...            % 2: 100
    'Context only', ...          % 3: 010
    'Value + Context', ...       % 4: 110
    'Interaction only', ...      % 5: 001
    'Value + Interaction', ...   % 6: 101
    'Context + Interaction', ... % 7: 011
    'All three'};                % 8: 111

group_cols = [ ...
    0.75 0.75 0.75; ...   % 1 None        — grey
    0.12 0.47 0.71; ...   % 2 V only      — blue
    0.84 0.15 0.16; ...   % 3 C only      — red
    0.00 0.60 0.60; ...   % 4 V+C         — teal
    0.58 0.10 0.70; ...   % 5 Int only    — purple
    0.80 0.20 0.60; ...   % 6 V+Int       — magenta
    0.93 0.50 0.10; ...   % 7 C+Int       — orange
    0.15 0.15 0.15];      % 8 All         — near-black
end

% ═════════════════════════════════════════════════════════════════════════════
%  plot_cell_type_analysis
% ═════════════════════════════════════════════════════════════════════════════
function plot_cell_type_analysis(R, save_opts)
if nargin < 2, save_opts = struct('enabled', false, 'dir', ''); end
% plot_cell_type_analysis
%
% Figure 1 — 3-circle Venn diagram of significant predictor overlaps.
% Figure 2 — Mean partial R² ± SEM for each of the 8 overlap groups.
% Figure 3 — 3-D scatter coloured by overlap group.

D        = R.pR2_win_neurons;   % [N × 3]
sig      = R.sig_flags;         % [N × 3] logical
gid      = R.cell_type;         % [N × 1] 1–8
gnames   = R.type_names;
gcols    = R.type_cols;

valid    = ~any(isnan(D), 2);
sv = sig(valid,1);  sc = sig(valid,2);  si = sig(valid,3);
gv = gid(valid);
Dv = D(valid,:);
N_valid  = sum(valid);

col_val = [0.12 0.47 0.71];
col_ctx = [0.84 0.15 0.16];
col_int = [0.58 0.10 0.70];

% ── Figure 1: Venn diagram (using venn.m) ────────────────────────────────────
% Zone counts: [V_only, C_only, Int_only, V&C, V&Int, C&Int, V&C&Int]
n_V    = sum(sv & ~sc & ~si);
n_C    = sum(~sv & sc & ~si);
n_I    = sum(~sv & ~sc & si);
n_VC   = sum(sv &  sc & ~si);
n_VI   = sum(sv & ~sc &  si);
n_CI   = sum(~sv & sc &  si);
n_VCI  = sum(sv &  sc &  si);
n_none = sum(~sv & ~sc & ~si);

Z = [n_V, n_C, n_I, n_VC, n_VI, n_CI, n_VCI];

fig1 = figure('Name', sprintf('Venn — %s', R.area), ...
    'Color', 'w', 'Position', [100 100 560 520]);
ax1 = axes(fig1);

% Draw Venn with venn.m (File Exchange)
[hv, venn_s] = venn(Z, 'FaceColor', {col_val, col_ctx, col_int}, ...
    'FaceAlpha', {0.25, 0.25, 0.25}, 'EdgeColor', {col_val, col_ctx, col_int}, ...
    'LineWidth', 2.5);

hold(ax1, 'on');
axis(ax1, 'equal');
axis(ax1, 'off');

% Annotate each zone with count and percentage
% venn_s.ZoneCentroid is [7×2] with [x y] for each zone
% Zone order: [A B C AB AC BC ABC] → [V C Int VC VI CI VCI]
zone_counts = Z;
for k = 1:7
    if zone_counts(k) == 0, continue; end
    xc = venn_s.ZoneCentroid(k, 1);
    yc = venn_s.ZoneCentroid(k, 2);
    pct = 100 * zone_counts(k) / N_valid;
    text(xc, yc, sprintf('%d\n(%.0f%%)', zone_counts(k), pct), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'Parent', ax1);
end

% Circle labels — use the circle patch XData/YData to find centres and radii
% hv(1..3) are the circle patches for V, C, Int
for k = 1:3
    xd = get(hv(k), 'XData');
    yd = get(hv(k), 'YData');
    cx_k = mean(xd);
    cy_k = mean(yd);
    rk   = (max(xd) - min(xd)) / 2;
    labels_k = {'Value','Context','Interaction'};
    label_cols_k = {col_val, col_ctx, col_int};
    % Place label above circle (k=1,2) or below (k=3)
    if k < 3
        text(cx_k, cy_k + rk*1.15, labels_k{k}, ...
            'FontSize', 11, 'FontWeight', 'bold', 'Color', label_cols_k{k}, ...
            'HorizontalAlignment', 'center', 'Parent', ax1);
    else
        text(cx_k, cy_k - rk*1.15, labels_k{k}, ...
            'FontSize', 11, 'FontWeight', 'bold', 'Color', label_cols_k{k}, ...
            'HorizontalAlignment', 'center', 'Parent', ax1);
    end
end

% "None" count in corner
ax_pos = axis(ax1);
text(ax_pos(2), ax_pos(4), sprintf('None: %d (%.0f%%)', n_none, 100*n_none/N_valid), ...
    'FontSize', 9, 'Color', [0.5 0.5 0.5], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Parent', ax1);

title(ax1, sprintf('Significant predictor overlaps — %s  (N=%d)', R.area, N_valid), ...
    'FontSize', 12);

% ── Figure 2: Mean R² per overlap group ──────────────────────────────────────
fig2 = figure('Name', sprintf('R² by overlap group — %s', R.area), ...
    'Color', 'w', 'Position', [640 100 760 420]);
ax2 = axes(fig2);
hold(ax2, 'on');

pred_cols   = [col_val; col_ctx; col_int];
pred_labels = {'Value', 'Context', 'V\timesC'};
bar_w   = 0.22;
offsets = linspace(-bar_w, bar_w, 3);
active_groups = unique(gv)';

for p = 1:3
    for tc = active_groups
        mask = gv == tc;
        mu  = nanmean(Dv(mask, p));
        sem = nanstd(Dv(mask, p)) / sqrt(sum(mask));
        bar(ax2, tc+offsets(p), mu, bar_w*0.9, ...
            'FaceColor', pred_cols(p,:), 'EdgeColor', 'none', 'FaceAlpha', 0.75);
        errorbar(ax2, tc+offsets(p), mu, sem, 'k', 'LineWidth', 1.2, 'CapSize', 4);
    end
end

yline(0,'LineStyle', '-','Color','k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(ax2, 'XTick', active_groups, 'XTickLabel', gnames(active_groups), 'FontSize', 8);
xtickangle(ax2, 25);
ylabel(ax2, 'Partial R²');
title(ax2, sprintf('Mean partial R² per overlap group — %s', R.area), 'FontSize', 11);
for p = 1:3
    patch(ax2, nan, nan, pred_cols(p,:), 'FaceAlpha', 0.75, 'EdgeColor', 'none', ...
        'DisplayName', pred_labels{p});
end
legend(ax2, 'Location', 'best', 'FontSize', 8);
box(ax2, 'off');

% ── Figure 3: 3-D scatter coloured by overlap group ──────────────────────────
fig3 = figure('Name', sprintf('Cell Types 3D — %s', R.area), ...
    'Color', 'w', 'Position', [100 560 520 460]);
ax3 = axes(fig3);
hold(ax3, 'on');

for tc = active_groups
    mask = gv == tc;
    scatter3(ax3, Dv(mask,1), Dv(mask,2), Dv(mask,3), 45, gcols(tc,:), 'filled', ...
        'MarkerFaceAlpha', 0.70, 'MarkerEdgeColor', 'none', 'DisplayName', gnames{tc});
    mu = nanmean(Dv(mask,:), 1);
    scatter3(ax3, mu(1), mu(2), mu(3), 150, gcols(tc,:), 'd', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
end

xl = xlim(ax3);  yl = ylim(ax3);  zl = zlim(ax3);
plot3(ax3, [0 0], yl, [0 0], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot3(ax3, xl, [0 0], [0 0], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot3(ax3, [0 0], [0 0], zl, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

xlabel(ax3, 'Partial R²  (Value)',   'FontSize', 10);
ylabel(ax3, 'Partial R²  (Context)', 'FontSize', 10);
zlabel(ax3, 'Partial R²  (V×C)',     'FontSize', 10);
title(ax3, sprintf('Overlap groups in R² space — %s  (N=%d)', R.area, N_valid), 'FontSize', 11);
legend(ax3, 'Location', 'best', 'FontSize', 8);
grid(ax3, 'on');
box(ax3, 'off');
view(ax3, 35, 25);
save_figure(fig1, save_opts, sprintf('%s_cell_type_venn',     R.area));
save_figure(fig2, save_opts, sprintf('%s_cell_type_r2_groups', R.area));
save_figure(fig3, save_opts, sprintf('%s_cell_type_3d',        R.area));
end

% ═════════════════════════════════════════════════════════════════════════════
%  compute_shared_axes
% ═════════════════════════════════════════════════════════════════════════════
function [w_shared, w_ctxdep, w_obj_win, w_loc_win] = ...
    compute_shared_axes(avg_w_obj, avg_w_loc, twin_obj, twin_loc)
% compute_shared_axes  Derive shared and context-dependent value axes.
%
% Both axes are fixed at their respective peak windows, consistent with the
% projection approach in plot_results.
%
%   Shared axis     = w_obj_win + w_loc_win  (both contexts encoding value same direction)
%   Context-dep axis = w_obj_win − w_loc_win  (difference; positive = more obj-like)
%
% Inputs:
%   avg_w_obj  [N × T]  bootstrap-averaged obj value axis (high_obj − low_obj)
%   avg_w_loc  [N × T]  same for location context
%   twin_obj   [1 × 2]  integer index range [start end] for obj peak window
%   twin_loc   [1 × 2]  integer index range [start end] for loc peak window
%
% Outputs:
%   w_shared   [N × 1]  normalized shared value axis
%   w_ctxdep   [N × 1]  normalized context-dependent value axis
%   w_obj_win  [N × 1]  unnormalized obj window vector
%   w_loc_win  [N × 1]  unnormalized loc window vector

w_obj_win = mean(avg_w_obj(:, twin_obj(1):twin_obj(2)), 2);   % [N × 1]
w_loc_win = mean(avg_w_loc(:, twin_loc(1):twin_loc(2)), 2);   % [N × 1]

raw_shared = w_obj_win + w_loc_win;
raw_ctxdep = w_obj_win - w_loc_win;

assert(norm(raw_shared) > 1e-12, 'compute_shared_axes: shared axis is near-zero — check inputs.');
assert(norm(raw_ctxdep) > 1e-12, 'compute_shared_axes: context-dep axis is near-zero — check inputs.');

w_shared = raw_shared / norm(raw_shared);
w_ctxdep = raw_ctxdep / norm(raw_ctxdep);

fprintf('  cos(shared, w_obj): %.3f   cos(shared, w_loc): %.3f\n', ...
    dot(w_shared, w_obj_win/norm(w_obj_win)), dot(w_shared, w_loc_win/norm(w_loc_win)));
fprintf('  cos(ctxdep, w_obj): %.3f   cos(ctxdep, w_loc): %.3f\n', ...
    dot(w_ctxdep, w_obj_win/norm(w_obj_win)), dot(w_ctxdep, w_loc_win/norm(w_loc_win)));
end

% ═════════════════════════════════════════════════════════════════════════════
%  plot_regression_variance
% ═════════════════════════════════════════════════════════════════════════════
function plot_regression_variance(R, save_opts)
if nargin < 2, save_opts = struct('enabled', false, 'dir', ''); end
% plot_regression_variance  Time-course of partial R² for each regressor.
% Shows population mean ± SEM. Peak windows marked with shading.

ta = R.time_axis;
lw = 2;

col_val = [0.12 0.47 0.71];   % blue  — value
col_ctx = [0.84 0.15 0.16];   % red   — context
col_int = [0.58 0.10 0.70];   % purple — interaction

figA = figure('Name', sprintf('Regression Variance Explained — %s', R.area), ...
    'Color', 'w', 'Position', [100 620 720 340]);
ax = axes(figA);
hold(ax, 'on');

% SEM bands (drawn first so lines sit on top)
fill_band(ax, ta, R.pR2_value,   R.pR2_sem_value,   col_val);
fill_band(ax, ta, R.pR2_context, R.pR2_sem_context, col_ctx);
fill_band(ax, ta, R.pR2_interact, R.pR2_sem_interact, col_int);

% Mean lines
plot(ax, ta, R.pR2_value,    '-', 'Color', col_val, 'LineWidth', lw, 'DisplayName', 'Value');
plot(ax, ta, R.pR2_context,  '-', 'Color', col_ctx, 'LineWidth', lw, 'DisplayName', 'Context');
plot(ax, ta, R.pR2_interact, '-', 'Color', col_int, 'LineWidth', lw, 'DisplayName', 'Value × Context');

% Peak window shading
yl = ylim(ax);
yl(1) = 0;
yl(2) = max(yl(2), 0.02);
ylim(ax, yl);

ta_obj = ta(round(R.twin_obj(1)):round(R.twin_obj(2)));
ta_loc = ta(round(R.twin_loc(1)):round(R.twin_loc(2)));
if numel(ta_obj) >= 2
    patch(ax, [ta_obj(1) ta_obj(end) ta_obj(end) ta_obj(1)], ...
        [yl(1) yl(1) yl(2) yl(2)], col_val, ...
        'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
if numel(ta_loc) >= 2
    patch(ax, [ta_loc(1) ta_loc(end) ta_loc(end) ta_loc(1)], ...
        [yl(1) yl(1) yl(2) yl(2)], col_ctx, ...
        'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

xline( 0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel(ax, 'Time (ms)');
ylabel(ax, 'Partial R²');
title(ax, sprintf('Variance explained per regressor — %s  (N=%d)', R.area, R.n_neurons));
legend(ax, 'Location', 'best', 'FontSize', 8);
box(ax, 'off');
uistack(findobj(ax, 'Type', 'line'), 'top');
save_figure(figA, save_opts, sprintf('%s_regression_variance', R.area));
end

% ═════════════════════════════════════════════════════════════════════════════
%  plot_proj_auc
% ═════════════════════════════════════════════════════════════════════════════
function plot_proj_auc(R, save_opts)
if nargin < 2, save_opts = struct('enabled', false, 'dir', ''); end
% plot_proj_auc  AUC for high-vs-low value discriminability along each axis.
%
% Panel A — obj-value axis: AUC(high_obj vs low_obj) and AUC(high_loc vs low_loc)
% Panel B — loc-value axis: same two contexts
%
% AUC > 0.5 means the axis separates the value conditions in that context.
% Chance = 0.5 (dashed line).

ta      = R.time_axis;
lw      = 2;
col_obj_ctx  = [0.84 0.15 0.16];   % red    — object context
col_loc_ctx  = [0.50 0.15 0.70];   % purple — location context

has_null = isfield(R, 'auc_null_ci_obj_on_obj');

fig = figure('Name', sprintf('Value Discriminability AUC — %s', R.area), ...
    'Color', 'w', 'Position', [100 200 820 380]);

% ── Panel A: projection onto obj-value axis ───────────────────────────────────
ax1 = subplot(1, 2, 1, 'Parent', fig);
hold(ax1, 'on');
if has_null
    fill_ci_lo_hi(ax1, ta, R.auc_null_ci_obj_on_obj(1,:), R.auc_null_ci_obj_on_obj(2,:), col_obj_ctx);  % within-obj null: red
    fill_ci_lo_hi(ax1, ta, R.auc_null_ci_loc_on_obj(1,:), R.auc_null_ci_loc_on_obj(2,:), col_loc_ctx);  % cross (loc on obj) null: purple
end
plot(ax1, ta, R.auc_obj_on_obj, '-', 'Color', col_obj_ctx, 'LineWidth', lw, 'DisplayName', 'Obj value (obj ctx)');
plot(ax1, ta, R.auc_loc_on_obj, '-', 'Color', col_loc_ctx, 'LineWidth', lw, 'DisplayName', 'Loc value (loc ctx)');
yline(0.5,'LineStyle', '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, 'HandleVisibility', 'off');
xline(0,  'LineStyle', '-',  'Color', 'k',           'LineWidth', 1, 'HandleVisibility', 'off');
xlabel(ax1, 'Time (ms)');
ylabel(ax1, 'AUC');
title(ax1, 'Discriminability along obj-value axis');
legend(ax1, 'Location', 'best', 'FontSize', 8);
ylim(ax1, [0.4 0.8]);
xlim(ax1, [ta(1) ta(end)]);
box(ax1, 'off');

% ── Panel B: projection onto loc-value axis ───────────────────────────────────
ax2 = subplot(1, 2, 2, 'Parent', fig);
hold(ax2, 'on');
if has_null
    fill_ci_lo_hi(ax2, ta, R.auc_null_ci_obj_on_loc(1,:), R.auc_null_ci_obj_on_loc(2,:), col_obj_ctx);  % cross (obj on loc) null: red
    fill_ci_lo_hi(ax2, ta, R.auc_null_ci_loc_on_loc(1,:), R.auc_null_ci_loc_on_loc(2,:), col_loc_ctx);  % within-loc null: purple
end
plot(ax2, ta, R.auc_obj_on_loc, '-', 'Color', col_obj_ctx, 'LineWidth', lw, 'DisplayName', 'Obj value (obj ctx)');
plot(ax2, ta, R.auc_loc_on_loc, '-', 'Color', col_loc_ctx, 'LineWidth', lw, 'DisplayName', 'Loc value (loc ctx)');
yline(0.5,'LineStyle', '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, 'HandleVisibility', 'off');
xline(0,  'LineStyle', '-',  'Color', 'k',           'LineWidth', 1, 'HandleVisibility', 'off');
xlabel(ax2, 'Time (ms)');
ylabel(ax2, 'AUC');
title(ax2, 'Discriminability along loc-value axis');
legend(ax2, 'Location', 'best', 'FontSize', 8);
ylim(ax2, [0.4 0.8]);
xlim(ax2, [ta(1) ta(end)]);
box(ax2, 'off');

sgtitle(fig, sprintf('Value discriminability AUC — %s  (N=%d)', R.area, R.n_neurons), ...
    'FontSize', 12);
save_figure(fig, save_opts, sprintf('%s_value_auc', R.area));
end

% ═════════════════════════════════════════════════════════════════════════════
%  auc_timecourse
% ═════════════════════════════════════════════════════════════════════════════
function auc = auc_timecourse(psth_hi, psth_lo, w, nTime)
% auc_timecourse  Per-time-bin AUC for high vs low value along projection axis w.
%
% For each time bin, pools w(n)*r_n(trial,t) across all N neurons and calls
% rocN(high_scores, low_scores).  Single-trial noise per neuron acts as the
% natural reference, giving within-context high / cross-context near-chance
% results when the axes are orthogonal.
N = numel(psth_hi);
whi = cell(N, 1);
wlo = cell(N, 1);
for n = 1:N
    whi{n} = w(n) * psth_hi{n};   % [trials_n × T]
    wlo{n} = w(n) * psth_lo{n};
end
auc = nan(1, nTime);
for t = 1:nTime
    sh = cellfun(@(x) x(:,t), whi, 'UniformOutput', false);
    sh = vertcat(sh{:});
    sl = cellfun(@(x) x(:,t), wlo, 'UniformOutput', false);
    sl = vertcat(sl{:});
    auc(t) = rocN(sh, sl);
end
end

% ═════════════════════════════════════════════════════════════════════════════
%  auc_null_timecourses
% ═════════════════════════════════════════════════════════════════════════════
function [null_oo, null_lo, null_ol, null_ll] = auc_null_timecourses( ...
    psth_hi_obj, psth_lo_obj, psth_hi_loc, psth_lo_loc, ...
    w_obj, w_loc, nTime, n_null)
% auc_null_timecourses  Shuffle-null for value discriminability AUC.
N = numel(psth_hi_obj);
null_oo = nan(n_null, nTime);
null_lo = nan(n_null, nTime);
null_ol = nan(n_null, nTime);
null_ll = nan(n_null, nTime);
for i = 1:n_null
    shi_obj = cell(N,1);  slo_obj = cell(N,1);
    shi_loc = cell(N,1);  slo_loc = cell(N,1);
    for n = 1:N
        [shi_obj{n}, slo_obj{n}] = shuffle_hi_lo(psth_hi_obj{n}, psth_lo_obj{n});
        [shi_loc{n}, slo_loc{n}] = shuffle_hi_lo(psth_hi_loc{n}, psth_lo_loc{n});
    end
    null_oo(i,:) = auc_timecourse(shi_obj, slo_obj, w_obj, nTime);
    null_lo(i,:) = auc_timecourse(shi_loc, slo_loc, w_obj, nTime);
    null_ol(i,:) = auc_timecourse(shi_obj, slo_obj, w_loc, nTime);
    null_ll(i,:) = auc_timecourse(shi_loc, slo_loc, w_loc, nTime);
end
end

function [hi_shuf, lo_shuf] = shuffle_hi_lo(hi, lo)
% shuffle_hi_lo  Pool hi+lo trials, randomly re-split preserving sizes.
n_hi    = size(hi, 1);
n_lo    = size(lo, 1);
all_t   = [hi; lo];
perm    = randperm(n_hi + n_lo);
hi_shuf = all_t(perm(1:n_hi),     :);
lo_shuf = all_t(perm(n_hi+1:end), :);
end

function fill_band(ax, x, y, sem, col)
% fill_band  Draw ±SEM shaded band around a timecourse.
x_fill = [x, fliplr(x)];
y_fill = [y(:)' + sem(:)', fliplr(y(:)' - sem(:)')];
fill(ax, x_fill, y_fill, col, 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
end

function fill_ci_lo_hi(ax, x, lo, hi, col)
% fill_ci_lo_hi  Draw asymmetric CI band from explicit lower/upper bounds.
x_fill = [x(:)', fliplr(x(:)')];
y_fill = [lo(:)', fliplr(hi(:)')];
fill(ax, x_fill, y_fill, col, 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
end

% ═════════════════════════════════════════════════════════════════════════════
%  plot_shared_ctxdep_axes
% ═════════════════════════════════════════════════════════════════════════════
function plot_shared_ctxdep_axes(R, save_opts)
if nargin < 2, save_opts = struct('enabled', false, 'dir', ''); end
% plot_shared_ctxdep_axes  Project 4 conditions onto shared and context-dep axes.
% Two-panel layout matching existing population projection figure style.

ta  = R.time_axis;
lw  = 2;
col_obj = [0.12 0.47 0.71];   % blue  — object context
col_loc = [0.84 0.15 0.16];   % red   — location context

figB = figure('Name', sprintf('Shared & Context-Dep Axes — %s', R.area), ...
    'Color', 'w', 'Position', [550 620 820 380]);

% Panel 1: shared value axis
ax1 = subplot(1, 2, 1, 'Parent', figB);
hold(ax1, 'on');
plot(ax1, ta, R.proj_shared.high_obj, '-',  'Color', col_obj, 'LineWidth', lw, 'DisplayName', 'High obj');
plot(ax1, ta, R.proj_shared.low_obj,  '--', 'Color', col_obj, 'LineWidth', lw, 'DisplayName', 'Low obj');
plot(ax1, ta, R.proj_shared.high_loc, '-',  'Color', col_loc, 'LineWidth', lw, 'DisplayName', 'High loc');
plot(ax1, ta, R.proj_shared.low_loc,  '--', 'Color', col_loc, 'LineWidth', lw, 'DisplayName', 'Low loc');
xline( 0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel(ax1, 'Time (ms)');
ylabel(ax1, 'Projection (a.u.)');
title(ax1, 'Shared value axis  (w_{obj} + w_{loc})');
legend(ax1, 'Location', 'best', 'FontSize', 8);
box(ax1, 'off');

% Panel 2: context-dependent value axis
ax2 = subplot(1, 2, 2, 'Parent', figB);
hold(ax2, 'on');
plot(ax2, ta, R.proj_ctxdep.high_obj, '-',  'Color', col_obj, 'LineWidth', lw, 'DisplayName', 'High obj');
plot(ax2, ta, R.proj_ctxdep.low_obj,  '--', 'Color', col_obj, 'LineWidth', lw, 'DisplayName', 'Low obj');
plot(ax2, ta, R.proj_ctxdep.high_loc, '-',  'Color', col_loc, 'LineWidth', lw, 'DisplayName', 'High loc');
plot(ax2, ta, R.proj_ctxdep.low_loc,  '--', 'Color', col_loc, 'LineWidth', lw, 'DisplayName', 'Low loc');
xline( 0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel(ax2, 'Time (ms)');
ylabel(ax2, 'Projection (a.u.)');
title(ax2, 'Context-dep value axis  (w_{obj} - w_{loc})');
legend(ax2, 'Location', 'best', 'FontSize', 8);
box(ax2, 'off');

sgtitle(figB, sprintf('Shared vs context-dependent value geometry — %s', R.area), ...
    'FontSize', 12);
save_figure(figB, save_opts, sprintf('%s_shared_ctxdep_axes', R.area));
end

% ═════════════════════════════════════════════════════════════════════════════
%  plot_2d_trajectory
% ═════════════════════════════════════════════════════════════════════════════
function plot_2d_trajectory(R, save_opts)
if nargin < 2, save_opts = struct('enabled', false, 'dir', ''); end
% plot_2d_trajectory  2-D state-space trajectory in the orthogonalized
% shared × context-dep plane.
%
% X-axis: projection onto shared value axis (orthogonalized)
% Y-axis: projection onto context-dep value axis (orthogonalized, Gram-Schmidt)
% Each of the 4 conditions traces a path through this 2-D plane as time
% progresses from -500 ms to 300 ms.  Time is encoded by colour (cool→warm)
% and the peak-window epoch is highlighted with thicker line segments.

ta  = R.time_axis;
nT  = numel(ta);
lw  = 1.5;

% Projections onto orthogonalized axes
S = R.proj_orth_shared;    % struct: .high_obj, .low_obj, .high_loc, .low_loc  [T × 1]
C = R.proj_orth_ctxdep;

col_obj = [0.12 0.47 0.71];   % blue  — object context
col_loc = [0.84 0.15 0.16];   % red   — location context

% Peak-window masks (used to draw thicker segments)
t_win_obj = false(1, nT);
t_win_obj(round(R.twin_obj(1)):round(R.twin_obj(2))) = true;
t_win_loc = false(1, nT);
t_win_loc(round(R.twin_loc(1)):round(R.twin_loc(2))) = true;

% Time colormap: cool (early) → warm (late), applied per-segment
cmap = cool(nT);

fig = figure('Name', sprintf('2D Trajectory — %s', R.area), ...
    'Color', 'w', 'Position', [100 100 560 520]);
ax = axes(fig);
hold(ax, 'on');
axis(ax, 'equal');

% Helper: draw a coloured-time trajectory with optional thick-window overlay
    function draw_traj(xs, ys, base_col, t_win_mask, ls)
        xs = xs(:);  ys = ys(:);
        for seg = 1:nT-1
            t_alpha = 0.3 + 0.7 * (seg / nT);   % fade in over time
            seg_col = cmap(seg,:) .* 0.4 + base_col .* 0.6;
            lw_seg  = lw;
            if t_win_mask(seg), lw_seg = lw * 2.5; end
            plot(ax, xs(seg:seg+1), ys(seg:seg+1), ls, ...
                'Color', [seg_col, t_alpha], 'LineWidth', lw_seg, ...
                'HandleVisibility', 'off');
        end
    end

% Draw trajectories for all 4 conditions
draw_traj(S.high_obj, C.high_obj, col_obj, t_win_obj, '-');
draw_traj(S.low_obj,  C.low_obj,  col_obj, t_win_obj, '--');
draw_traj(S.high_loc, C.high_loc, col_loc, t_win_loc, '-');
draw_traj(S.low_loc,  C.low_loc,  col_loc, t_win_loc, '--');

% Start (circle) and end (triangle) markers
t_start = 1;
t_end   = nT;
marker_sz = 80;
scatter(ax, S.high_obj(t_start), C.high_obj(t_start), marker_sz, col_obj, 'o', 'filled');
scatter(ax, S.low_obj(t_start),  C.low_obj(t_start),  marker_sz, col_obj, 'o', 'filled');
scatter(ax, S.high_loc(t_start), C.high_loc(t_start), marker_sz, col_loc, 'o', 'filled');
scatter(ax, S.low_loc(t_start),  C.low_loc(t_start),  marker_sz, col_loc, 'o', 'filled');

scatter(ax, S.high_obj(t_end), C.high_obj(t_end), marker_sz, col_obj, '^', 'filled');
scatter(ax, S.low_obj(t_end),  C.low_obj(t_end),  marker_sz, col_obj, '^', 'filled');
scatter(ax, S.high_loc(t_end), C.high_loc(t_end), marker_sz, col_loc, '^', 'filled');
scatter(ax, S.low_loc(t_end),  C.low_loc(t_end),  marker_sz, col_loc, '^', 'filled');

% Stimulus-onset dot (t=0)
[~, t0] = min(abs(ta));
scatter(ax, S.high_obj(t0), C.high_obj(t0), 55, col_obj, 'pentagram', 'filled');
scatter(ax, S.low_obj(t0),  C.low_obj(t0),  55, col_obj, 'pentagram', 'filled');
scatter(ax, S.high_loc(t0), C.high_loc(t0), 55, col_loc, 'pentagram', 'filled');
scatter(ax, S.low_loc(t0),  C.low_loc(t0),  55, col_loc, 'pentagram', 'filled');

% Axes cross at origin
xl = xlim(ax);  yl = ylim(ax);
plot(ax, xl, [0 0], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(ax, [0 0], yl, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Legend proxies
plot(ax, nan, nan, '-',  'Color', col_obj, 'LineWidth', 2, 'DisplayName', 'High obj');
plot(ax, nan, nan, '--', 'Color', col_obj, 'LineWidth', 2, 'DisplayName', 'Low obj');
plot(ax, nan, nan, '-',  'Color', col_loc, 'LineWidth', 2, 'DisplayName', 'High loc');
plot(ax, nan, nan, '--', 'Color', col_loc, 'LineWidth', 2, 'DisplayName', 'Low loc');
scatter(ax, nan, nan, 40, [0.4 0.4 0.4], 'o', 'filled', 'DisplayName', 'Start (−500 ms)');
scatter(ax, nan, nan, 40, [0.4 0.4 0.4], '^', 'filled', 'DisplayName', 'End (+300 ms)');
scatter(ax, nan, nan, 40, [0.4 0.4 0.4], 'pentagram', 'filled', 'DisplayName', 'Stim onset (0 ms)');

legend(ax, 'Location', 'best', 'FontSize', 8);
xlabel(ax, 'Shared value axis (orth.)',      'FontSize', 11);
ylabel(ax, 'Context-dep value axis (orth.)', 'FontSize', 11);
title(ax, sprintf('Population trajectory — %s', R.area), 'FontSize', 12);
box(ax, 'off');

% Colorbar to show time
colormap(ax, cool);
cb = colorbar(ax);
cb.Label.String = 'Time';
clim(ax, [ta(1) ta(end)]);
cb.Ticks = [ta(1), 0, ta(end)];
cb.TickLabels = {sprintf('%d ms', ta(1)), '0', sprintf('%d ms', ta(end))};
save_figure(fig, save_opts, sprintf('%s_2d_trajectory', R.area));
end

% ═════════════════════════════════════════════════════════════════════════════
%  compute_decoder_accuracy
% ═════════════════════════════════════════════════════════════════════════════
function [acc_wo, acc_wl, acc_o2l, acc_l2o,acc_wo_sem, acc_wl_sem, acc_o2l_sem, acc_l2o_sem] = compute_decoder_accuracy( ...
    psth_hi_obj, psth_lo_obj, psth_hi_loc, psth_lo_loc, ...
    w_obj_dec, w_loc_dec, K_fold)
% compute_decoder_accuracy  Fixed-axis linear decoder tied to the orthogonality axes.
%
% w_obj_dec and w_loc_dec are the peak-window-normalised axis vectors whose
% angle IS theta. Using these (not per-bin avg_w(:,t)) guarantees:
%   dot(w_obj_dec, w_loc_dec) = cos(theta)
% so if theta ≈ 90°, loc class means project equally onto w_obj_dec and
% cross-decoding is near chance.
%
% K-fold CV calibrates the threshold from training folds only — no test
% trial sets its own boundary.
%
% Inputs:
%   psth_hi/lo_obj/loc  {N x 1} cell, each [n_trials x T]
%   w_obj_dec / w_loc_dec  [N x 1] unit vectors (peak-window axes)
%   K_fold              number of CV folds
%
% Outputs:
%   acc_wo   [1 x T]  within-obj:  w_obj_dec axis, K-fold threshold on obj trials
%   acc_wl   [1 x T]  within-loc:  w_loc_dec axis, K-fold threshold on loc trials
%   acc_o2l  [1 x T]  cross obj→loc: w_obj_dec axis, K-fold threshold on loc trials
%   acc_l2o  [1 x T]  cross loc→obj: w_loc_dec axis, K-fold threshold on obj trials

N  = numel(psth_hi_obj);
nT = size(psth_hi_obj{1}, 2);

[Xt_hi_obj, Xt_lo_obj, n_min_obj] = build_trial_array(psth_hi_obj, psth_lo_obj, N, nT);
[Xt_hi_loc, Xt_lo_loc, n_min_loc] = build_trial_array(psth_hi_loc, psth_lo_loc, N, nT);

folds_obj = assign_folds(n_min_obj, K_fold);
folds_loc = assign_folds(n_min_loc, K_fold);

% Pre-project all trials onto both fixed axes: [n_min x T]
% Xt arrays are [n_min x N x T]; pagemtimes applies [N x 1] multiply per page.
proj_hi_o_wo = squeeze(pagemtimes(Xt_hi_obj, w_obj_dec));   % [n_min_obj x T]
proj_lo_o_wo = squeeze(pagemtimes(Xt_lo_obj, w_obj_dec));
proj_hi_l_wl = squeeze(pagemtimes(Xt_hi_loc, w_loc_dec));   % [n_min_loc x T]
proj_lo_l_wl = squeeze(pagemtimes(Xt_lo_loc, w_loc_dec));
proj_hi_l_wo = squeeze(pagemtimes(Xt_hi_loc, w_obj_dec));   % cross: loc trials on obj axis
proj_lo_l_wo = squeeze(pagemtimes(Xt_lo_loc, w_obj_dec));
proj_hi_o_wl = squeeze(pagemtimes(Xt_hi_obj, w_loc_dec));   % cross: obj trials on loc axis
proj_lo_o_wl = squeeze(pagemtimes(Xt_lo_obj, w_loc_dec));

acc_wo  = zeros(1, nT);
acc_wl  = zeros(1, nT);
acc_o2l = zeros(1, nT);
acc_l2o = zeros(1, nT);

for t = 1:nT
    % ── Within obj: w_obj_dec axis, K-fold threshold on obj trials ───────────
    correct_hi_o = false(n_min_obj, 1);
    correct_lo_o = false(n_min_obj, 1);
    for k = 1:K_fold
        test  = (folds_obj == k);
        train = ~test;
        thr = (mean(proj_hi_o_wo(train,t)) + mean(proj_lo_o_wo(train,t))) / 2;
        correct_hi_o(test) = proj_hi_o_wo(test,t) > thr;
        correct_lo_o(test) = proj_lo_o_wo(test,t) < thr;
    end
    acc_wo(t) = mean([correct_hi_o; correct_lo_o]);
    acc_wo_sem(t) = std([correct_hi_o; correct_lo_o])/sqrt(length(correct_hi_o)*2);

    % ── Within loc: w_loc_dec axis, K-fold threshold on loc trials ───────────
    correct_hi_l = false(n_min_loc, 1);
    correct_lo_l = false(n_min_loc, 1);
    for k = 1:K_fold
        test  = (folds_loc == k);
        train = ~test;
        thr = (mean(proj_hi_l_wl(train,t)) + mean(proj_lo_l_wl(train,t))) / 2;
        correct_hi_l(test) = proj_hi_l_wl(test,t) > thr;
        correct_lo_l(test) = proj_lo_l_wl(test,t) < thr;
    end
    acc_wl(t) = mean([correct_hi_l; correct_lo_l]);
    acc_wl_sem(t) = std([correct_hi_l; correct_lo_l])/sqrt(length(correct_hi_l)*2);

    % ── Cross obj→loc: w_obj_dec axis, K-fold threshold on loc trials ─────────
    % cos(theta)≈0 → loc class means project equally → chance.
    correct_hi_l_x = false(n_min_loc, 1);
    correct_lo_l_x = false(n_min_loc, 1);
    for k = 1:K_fold
        test  = (folds_loc == k);
        train = ~(folds_obj == k);
        thr = (mean(proj_hi_o_wo(train,t)) + mean(proj_lo_o_wo(train,t))) / 2;
        
        correct_hi_l_x(test) = proj_hi_l_wo(test,t) > thr;
        correct_lo_l_x(test) = proj_lo_l_wo(test,t) < thr;
    end
    acc_o2l(t) = mean([correct_hi_l_x; correct_lo_l_x]);
    acc_o2l_sem(t) = std([correct_hi_l_x; correct_lo_l_x])/sqrt(length(correct_hi_l_x)*2);

    % ── Cross loc→obj: w_loc_dec axis, K-fold threshold on obj trials ─────────
    correct_hi_o_x = false(n_min_obj, 1);
    correct_lo_o_x = false(n_min_obj, 1);
    for k = 1:K_fold
        test  = (folds_obj == k);
        train = ~(folds_loc == k);
        thr = (mean(proj_hi_l_wl(train,t)) + mean(proj_lo_l_wl(train,t))) / 2;
        
        correct_hi_o_x(test) = proj_hi_o_wl(test,t) > thr;
        correct_lo_o_x(test) = proj_lo_o_wl(test,t) < thr;
    end
    acc_l2o(t) = mean([correct_hi_o_x; correct_lo_o_x]);
    acc_l2o_sem(t) = std([correct_hi_o_x; correct_lo_o_x])/sqrt(length(correct_hi_o_x)*2);
end
end

% ═════════════════════════════════════════════════════════════════════════════
%  build_trial_array
% ═════════════════════════════════════════════════════════════════════════════
function [Xt_hi, Xt_lo, n_min] = build_trial_array(psth_hi, psth_lo, N, nT)
% build_trial_array  Build [n_min x N x T] arrays from {N x 1} cell arrays.
% Samples n_min = min(trials/class across neurons) trials per neuron.

n_min_hi = min(cellfun(@(x) size(x,1), psth_hi));
n_min_lo = min(cellfun(@(x) size(x,1), psth_lo));
%n_min    = min(n_min_hi, n_min_lo);
n_min = 10;

Xt_hi = zeros(n_min, N, nT);
Xt_lo = zeros(n_min, N, nT);
for n = 1:N
    idx_hi = randperm(size(psth_hi{n}, 1), n_min);
    idx_lo = randperm(size(psth_lo{n}, 1), n_min);
    Xt_hi(:, n, :) = psth_hi{n}(idx_hi, :);
    Xt_lo(:, n, :) = psth_lo{n}(idx_lo, :);
end
end

% ═════════════════════════════════════════════════════════════════════════════
%  assign_folds
% ═════════════════════════════════════════════════════════════════════════════
function folds = assign_folds(n_trials, K)
% assign_folds  Return [n_trials x 1] vector of fold indices 1..K,
% balanced so each fold gets floor(n_trials/K) or ceil(n_trials/K) trials.
folds = zeros(n_trials, 1);
idx   = randperm(n_trials);
for k = 1:K
    members = idx(k:K:end);
    folds(members) = k;
end
end

% ═════════════════════════════════════════════════════════════════════════════
%  plot_decoding_accuracy
% ═════════════════════════════════════════════════════════════════════════════
function plot_decoding_accuracy(R, save_opts)
if nargin < 2, save_opts = struct('enabled', false, 'dir', ''); end
% plot_decoding_accuracy  Within- and cross-context value decoding accuracy.
%
% Within-context: K-fold CV — axis re-estimated from training folds at each
%   time bin; tested on held-out fold. Eliminates train-test leakage.
% Cross-context: axis trained on all trials of one context, tested on the
%   other context (disjoint trial sets — no leakage by design).

ta = R.time_axis;
lw = 2;

col_wo  = [0.12 0.47 0.71];   % blue   — within obj
col_o2l = [0.17 0.63 0.17];   % green  — obj axis → loc trials
col_wl  = [0.84 0.15 0.16];   % red    — within loc
col_l2o = [0.58 0.40 0.74];   % purple — loc axis → obj trials

ta_obj = ta(round(R.twin_obj(1)):round(R.twin_obj(2)));
ta_loc = ta(round(R.twin_loc(1)):round(R.twin_loc(2)));

% ── Figure 1: Obj-context — within obj + obj axis → loc trials ───────────────
fig1 = figure('Name', sprintf('Value Decoder — Obj axis (%s)', R.area), ...
    'Color', 'w', 'Position', [100 620 760 390]);
ax1 = axes(fig1);
hold(ax1, 'on');

if isfield(R, 'acc_within_obj_sem')
    fill_band(ax1, ta, R.acc_within_obj    * 100, R.acc_within_obj_sem    * 100, col_wo);
    fill_band(ax1, ta, R.acc_cross_obj2loc * 100, R.acc_cross_obj2loc_sem * 100, col_o2l);
end
plot(ax1, ta, R.acc_within_obj    * 100, '-',  'Color', col_wo,  'LineWidth', lw, ...
    'DisplayName', 'Within obj');
plot(ax1, ta, R.acc_cross_obj2loc * 100, '--', 'Color', col_o2l, 'LineWidth', lw, ...
    'DisplayName', 'Obj axis \rightarrow Loc trials');

yline(50, 'LineStyle',':','Color','k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline( 0, 'LineStyle','-','Color','k', 'LineWidth', 1,   'HandleVisibility', 'off');
yl1 = [35 105];
if numel(ta_obj) >= 2
    patch(ax1, [ta_obj(1) ta_obj(end) ta_obj(end) ta_obj(1)], ...
        [yl1(1) yl1(1) yl1(2) yl1(2)], col_wo, ...
        'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
xlabel(ax1, 'Time (ms)');
ylabel(ax1, 'Decoding accuracy (%)');
title(ax1, sprintf('Obj axis — within obj & cross to loc  (%s, N=%d)', R.area, R.n_neurons));
legend(ax1, 'Location', 'best', 'FontSize', 8);
ylim(ax1, yl1);
xlim(ax1, [ta(1) ta(end)]);
box(ax1, 'off');
uistack(findobj(ax1, 'Type', 'line'), 'top');

% ── Figure 2: Loc-context — within loc + loc axis → obj trials ───────────────
fig2 = figure('Name', sprintf('Value Decoder — Loc axis (%s)', R.area), ...
    'Color', 'w', 'Position', [880 620 760 390]);
ax2 = axes(fig2);
hold(ax2, 'on');

if isfield(R, 'acc_within_loc_sem')
    fill_band(ax2, ta, R.acc_within_loc    * 100, R.acc_within_loc_sem    * 100, col_wl);
    fill_band(ax2, ta, R.acc_cross_loc2obj * 100, R.acc_cross_loc2obj_sem * 100, col_l2o);
end
plot(ax2, ta, R.acc_within_loc    * 100, '-',  'Color', col_wl,  'LineWidth', lw, ...
    'DisplayName', 'Within loc');
plot(ax2, ta, R.acc_cross_loc2obj * 100, '--', 'Color', col_l2o, 'LineWidth', lw, ...
    'DisplayName', 'Loc axis \rightarrow Obj trials');

yline(50, 'LineStyle',':','Color','k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline( 0, 'LineStyle','-','Color','k', 'LineWidth', 1,   'HandleVisibility', 'off');
yl2 = [35 105];
if numel(ta_loc) >= 2
    patch(ax2, [ta_loc(1) ta_loc(end) ta_loc(end) ta_loc(1)], ...
        [yl2(1) yl2(1) yl2(2) yl2(2)], col_wl, ...
        'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
xlabel(ax2, 'Time (ms)');
ylabel(ax2, 'Decoding accuracy (%)');
title(ax2, sprintf('Loc axis — within loc & cross to obj  (%s, N=%d)', R.area, R.n_neurons));
legend(ax2, 'Location', 'best', 'FontSize', 8);
ylim(ax2, yl2);
xlim(ax2, [ta(1) ta(end)]);
box(ax2, 'off');
uistack(findobj(ax2, 'Type', 'line'), 'top');
save_figure(fig1, save_opts, sprintf('%s_decoder_obj_axis', R.area));
save_figure(fig2, save_opts, sprintf('%s_decoder_loc_axis', R.area));
end

% ═════════════════════════════════════════════════════════════════════════════
%  save_figure
% ═════════════════════════════════════════════════════════════════════════════
function save_figure(fig, save_opts, fname)
% save_figure  Export figure as a vectorized EPS file (Adobe Illustrator-compatible).
%
% Uses the painters renderer so all lines, text, and patches are vector objects.
% Transparency (FaceAlpha) is flattened to the background colour in EPS —
% this is an EPS format limitation, not a rendering bug.
if ~isfield(save_opts, 'enabled') || ~save_opts.enabled
    return;
end
if ~isfolder(save_opts.dir)
    mkdir(save_opts.dir);
end
set(fig, 'Renderer', 'painters');
out_path = fullfile(save_opts.dir, [fname '.eps']);
print(fig, out_path, '-depsc', '-painters');
fprintf('  Saved: %s\n', out_path);
end
