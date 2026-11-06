% AUC_ContextDiscrimination.m
%
% Purpose:
%   Compute context-discriminability AUC between Stable and OneDR tasks
%   for RF-location trials, separately for high-value and low-value conditions.
%   Reports:
%     (1) Full AUCPSTH time series (mean ± SEM across neurons)
%     (2) Peak AUC in the 100–200 ms post-target window
%
%   At each time bin, AUC is computed trial-by-trial (rocN) comparing
%   stable-task responses against oneDR-task responses.
%   AUC > 0.5  →  stable firing rate tends to exceed oneDR.
%   Neurons are rectified per-neuron (flip direction set by sign of mean
%   AUC in the 100–200 ms window) so the population average reflects
%   discriminability magnitude rather than direction.
%
% Inputs (loaded from file):
%   PopPSTH_RF_FEF.mat, PopPSTH_RF_SC.mat
%   Expected fields in each file:
%     psth_targeton_raw_stable_RF_high  {N×1}  [n_trials × T_stable]
%     psth_targeton_raw_stable_RF_low   {N×1}
%     psth_targeton_raw_oneDR_RF_high   {N×1}  [n_trials × T_oneDR]
%     psth_targeton_raw_oneDR_RF_low    {N×1}
%     PSTH_Time_TgtOn_stable            [1 × T_stable]  ms from target onset
%     PSTH_Time_TgtOn_OneDR             [1 × T_oneDR]
%
% Outputs:
%   AUC_Context_Results.mat  — saved in ResultsDir
%   AUC_Context_Discrimination.png / .fig  — saved in ResultsDir/Figures/

%% -------------------------------------------------------------------------
%  Paths
% -------------------------------------------------------------------------
ResultsDir = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/ObjectLocationPSTH';
FigDir     = fullfile(ResultsDir, 'Figures');
if ~exist(FigDir, 'dir')
    mkdir(FigDir);
end

% Add Apps to path so PreprocessPSTH and rocN are accessible
AppsRoot = fileparts(fileparts(mfilename('fullpath')));  % .../Apps/
addpath(genpath(AppsRoot));

%% -------------------------------------------------------------------------
%  Parameters
% -------------------------------------------------------------------------
SmoothKernel    = 'gaussian';   % kernel type for PreprocessPSTH
SmoothSigma_ms  = 30;           % Gaussian sigma (ms) — matches CrossDecoder pipeline
PeakWindow_ms   = [100 200];    % target-period window for peak AUC and flip direction
PlotWindow_ms   = [-200 400];   % time range shown in AUCPSTH plots

areas = {'FEF', 'SC'};
nAreas = numel(areas);

%% -------------------------------------------------------------------------
%  Main loop: compute context AUC for each brain area
% -------------------------------------------------------------------------
Results = struct();

for aIdx = 1:nAreas
    area = areas{aIdx};
    fprintf('\n=== %s ===\n', area);

    % --- Load PopPSTH data ---------------------------------------------------
    fname = fullfile(ResultsDir, area, sprintf('PopPSTH_RF_%s.mat', area));
    fprintf('Loading %s ...\n', fname);
    D = load(fname);

    psth_stable_high = D.psth_targeton_raw_stable_RF_high;
    psth_stable_low  = D.psth_targeton_raw_stable_RF_low;
    psth_oneDR_high  = D.psth_targeton_raw_oneDR_RF_high;
    psth_oneDR_low   = D.psth_targeton_raw_oneDR_RF_low;
    time_stable      = D.PSTH_Time_TgtOn_stable;
    time_oneDR       = D.PSTH_Time_TgtOn_OneDR;
    modutype         = D.ModuType_Sel;   % [N×1]: 1=Fast, 2=Slow, 3=Inhib, 4=NoEffect

    nNeurons = numel(psth_stable_high);
    fprintf('  %d neurons loaded\n', nNeurons);

    % --- Smooth PSTHs --------------------------------------------------------
    fprintf('  Smoothing PSTHs (%s %d ms)...\n', SmoothKernel, SmoothSigma_ms);
    [psth_stable_high_s, time_stable_s] = PreprocessPSTH(psth_stable_high, time_stable, ...
        'kernel', SmoothKernel, 'win_ms', SmoothSigma_ms);
    [psth_stable_low_s,  ~             ] = PreprocessPSTH(psth_stable_low,  time_stable, ...
        'kernel', SmoothKernel, 'win_ms', SmoothSigma_ms);
    [psth_oneDR_high_s,  time_oneDR_s  ] = PreprocessPSTH(psth_oneDR_high,  time_oneDR,  ...
        'kernel', SmoothKernel, 'win_ms', SmoothSigma_ms);
    [psth_oneDR_low_s,   ~             ] = PreprocessPSTH(psth_oneDR_low,   time_oneDR,  ...
        'kernel', SmoothKernel, 'win_ms', SmoothSigma_ms);

    % Crop stable to [-500, 300] ms
t_sel_stable = time_stable >= -500 & time_stable <= 300;
time_axis    = time_stable(t_sel_stable);
for n = 1:numel(psth_high_value_object_RF)
    psth_stable_high_s{n} = psth_stable_high_s{n}(:, t_sel_stable);
    psth_stable_low_s{n}  = psth_stable_low_s{n}(:,  t_sel_stable);
end

% Crop oneDR to the same [-500, 300] ms window
t_sel_oneDR = time_oneDR >= -500 & time_oneDR <= 300;
for n = 1:numel(psth_high_value_location_RF)
    psth_oneDR_high_s{n} = psth_oneDR_high_s{n}(:, t_sel_oneDR);
    psth_oneDR_low_s{n}  = psth_oneDR_low_s{n}(:,  t_sel_oneDR);
end



    % --- Find common time axis -----------------------------------------------
    % Both time axes are aligned to target onset; use their intersection so
    % AUC is computed at the same time points for both tasks.
    time_common = time_axis;
    idx_stable = t_sel_stable;
    idx_oneDR = t_sel_oneDR;
    %[time_common, idx_stable, idx_oneDR] = intersect(time_stable_s, time_oneDR_s);
    nT = numel(time_common);
    fprintf('  Common time axis: %d bins (%.0f to %.0f ms)\n', ...
        nT, time_common(1), time_common(end));

    % --- Identify valid neurons ----------------------------------------------
    % A neuron is valid if both stable and oneDR cell arrays are non-empty
    % (all neurons in PopPSTH already passed the ≥10-trial filter, but guard
    % against any unexpected empty cells).
    valid = false(nNeurons, 1);
    for i = 1:nNeurons
        valid(i) = ~isempty(psth_stable_high_s{i}) && ...
                   ~isempty(psth_stable_low_s{i})  && ...
                   ~isempty(psth_oneDR_high_s{i})  && ...
                   ~isempty(psth_oneDR_low_s{i});
    end
    nValid = sum(valid);
    fprintf('  %d / %d neurons valid\n', nValid, nNeurons);

    % --- Compute context AUC per neuron per time bin -------------------------
    % rocN(x, y): x = stable (positive class), y = oneDR (negative class)
    % AUC > 0.5 → stable response > oneDR response
    fprintf('  Computing context AUC PSTH (%d neurons × %d bins)...\n', nValid, nT);
    auc_high = NaN(nNeurons, nT);
    auc_low  = NaN(nNeurons, nT);

    validIdx = find(valid);
    for k = 1:numel(validIdx)
        i = validIdx(k);
        if mod(k, 20) == 0
            fprintf('    neuron %d / %d\n', k, numel(validIdx));
        end
        for t = 1:nT
            x_h = psth_stable_high_s{i}(:, t);
            y_h = psth_oneDR_high_s{i}(:, t);
            auc_high(i, t) = rocN(x_h, y_h);

            x_l = psth_stable_low_s{i}(:, t);
            y_l = psth_oneDR_low_s{i}(:, t);
            auc_low(i, t) = rocN(x_l, y_l);
        end
    end

    % --- Rectify per-neuron --------------------------------------------------
    % Flip direction determined from mean AUC in PeakWindow_ms so that the
    % entire time series of each neuron points in the same direction.
    % Peak discriminability = max(|AUC - 0.5|) in that window.
    peak_mask = time_common >= PeakWindow_ms(1) & time_common <= PeakWindow_ms(2);

    mean_tgt_high = nanmean(auc_high(:, peak_mask), 2);  % [N×1] mean AUC in window
    mean_tgt_low  = nanmean(auc_low(:, peak_mask),  2);

    flip_high = mean_tgt_high < 0.5;  % these neurons: stable < oneDR; flip entire series
    flip_low  = mean_tgt_low  < 0.5;

    auc_high_rect = auc_high;
    auc_low_rect  = auc_low;
    auc_high_rect(flip_high, :) = 1 - auc_high_rect(flip_high, :);
    auc_low_rect(flip_low,  :) = 1 - auc_low_rect(flip_low,  :);

    % --- Peak AUC in target window -------------------------------------------
    % max(|AUC - 0.5|) across time bins in PeakWindow_ms; per neuron
    peak_auc_high = nanmax(abs(auc_high(:, peak_mask) - 0.5), [], 2);
    peak_auc_low  = nanmax(abs(auc_low(:, peak_mask)  - 0.5), [], 2);

    % Valid neurons only (exclude NaN from rocN for empty / all-NaN inputs)
    valid_auc = valid & ~any(isnan(auc_high), 2) & ~any(isnan(auc_low), 2);
    nFinal = sum(valid_auc);
    fprintf('  %d neurons with complete AUC (no NaN)\n', nFinal);

    % --- Population mean ± SEM -----------------------------------------------
    mean_high = nanmean(auc_high_rect(valid_auc, :), 1);    % [1 × nT]
    sem_high  = nanstd( auc_high_rect(valid_auc, :), 0, 1) / sqrt(nFinal);
    mean_low  = nanmean(auc_low_rect( valid_auc, :), 1);
    sem_low   = nanstd( auc_low_rect( valid_auc, :), 0, 1) / sqrt(nFinal);


%{
    figure
    plot(mean_high,'-r')
    hold on 
    plot(mean_low,'-b')
%}
%    keyboard
    % --- Store results -------------------------------------------------------
    R = struct();
    R.time_common    = time_common;
    R.auc_high       = auc_high;          % raw AUC [N × nT]
    R.auc_low        = auc_low;
    R.auc_high_rect  = auc_high_rect;     % rectified [N × nT]
    R.auc_low_rect   = auc_low_rect;
    R.mean_high      = mean_high;         % population mean [1 × nT]
    R.sem_high       = sem_high;
    R.mean_low       = mean_low;
    R.sem_low        = sem_low;
    R.peak_auc_high  = peak_auc_high;     % per-neuron peak [N × 1]
    R.peak_auc_low   = peak_auc_low;
    R.flip_high      = flip_high;
    R.flip_low       = flip_low;
    R.valid_auc      = valid_auc;
    R.n_neurons      = nFinal;

    % --- Targeted Dimensionality Reduction (TDR) ----------------------------
    % Regression model (no interaction):
    %   r = β₀ + β_context * context + β_value * value+ β_interaction *
    %   context*value
    %
    % Condition codes (context: stable=-1, oneDR=+1; value: high=+1, low=-1):
    %   cond 1: stable-high  [-1, +1]
    %   cond 2: stable-low   [-1, -1]
    %   cond 3: oneDR-high   [+1, +1]
    %   cond 4: oneDR-low    [+1, -1]
    %
    % Axes are fit on the key window (PeakWindow_ms), then applied to the
    % full time axis. d_value is Gram-Schmidt orthogonalized w.r.t. d_context.
    fprintf('  Running TDR...\n');

    % Step 1: Condition-averaged PSTHs on common time axis [N × nT × 4]
    cond_avg = NaN(nNeurons, nT, 4);
    for i = find(valid)'
        cond_avg(i,:,1) = nanmean(psth_stable_high_s{i}(:, :), 1);
        cond_avg(i,:,2) = nanmean(psth_stable_low_s{i}(:,  :), 1);
        cond_avg(i,:,3) = nanmean(psth_oneDR_high_s{i}(:,  :),  1);
        cond_avg(i,:,4) = nanmean(psth_oneDR_low_s{i}(:,   :),  1);
    end

    % Step 2: Normalize each neuron by std of condition means in key window
    % so high-FR neurons do not dominate the regression axes
    norm_mask  = time_common >= PeakWindow_ms(1) & time_common <= PeakWindow_ms(2);
    win_means  = squeeze(nanmean(cond_avg(:, norm_mask, :), 2));  % [N × 4]
    neuron_std = std(win_means, 0, 2);                             % [N × 1]
    neuron_std(neuron_std < 1e-6 | isnan(neuron_std)) = 1;        % guard zero/NaN
    cond_avg_norm = cond_avg ./ neuron_std;                        % [N × nT × 4]

    % Step 3: Fit full factorial regression on key-window averages
    % Design matrix: intercept, context, value, context×value (interaction)
    % Conditions: [stable-high, stable-low, oneDR-high, oneDR-low]
    %   context:     [-1, -1, +1, +1]
    %   value:       [+1, -1, +1, -1]
    %   interaction: [-1, +1, +1, -1]  (element-wise product of context × value)
    % All four contrast columns are mutually orthogonal in this balanced design,
    % so the solution is exact (4 parameters, 4 conditions) and each beta is
    % independent — adding the interaction does not change context or value betas.
    X_tdr = [ones(4,1), [-1;-1;1;1], [1;-1;1;-1], [-1;1;1;-1]];  % [4 × 4]
    F_win = squeeze(nanmean(cond_avg_norm(valid_auc, norm_mask, :), 2));  % [nFinal × 4]
    beta  = X_tdr \ F_win';      % [4 × nFinal], exact solution
    beta_context     = beta(2,:)';  % [nFinal × 1]
    beta_value       = beta(3,:)';
    beta_interaction = beta(4,:)';

    % Step 4: Normalize to unit vectors
    d_context   = beta_context     / norm(beta_context);
    d_value_raw = beta_value       / norm(beta_value);
    d_inter_raw = beta_interaction / norm(beta_interaction);

    % Step 5: Gram-Schmidt — sequential orthogonalization
    % d_value ⊥ d_context
    d_value = d_value_raw - (d_context' * d_value_raw) * d_context;
    d_value = d_value / norm(d_value);
    % d_interaction ⊥ d_context and ⊥ d_value
    d_inter = d_inter_raw - (d_context' * d_inter_raw) * d_context ...
                          - (d_value'   * d_inter_raw) * d_value;
    d_inter = d_inter / norm(d_inter);

    fprintf('  TDR axis alignments (|cos| before orth):\n');
    fprintf('    cos(context, value)       = %.3f\n', abs(d_context' * d_value_raw));
    fprintf('    cos(context, interaction) = %.3f\n', abs(d_context' * d_inter_raw));
    fprintf('    cos(value,   interaction) = %.3f\n', abs(d_value_raw' * d_inter_raw));

    % Step 6: Project all conditions onto the three axes [4 × nT]
    cond_valid       = cond_avg_norm(valid_auc, :, :);  % [nFinal × nT × 4]
    proj_context     = NaN(4, nT);
    proj_value       = NaN(4, nT);
    proj_interaction = NaN(4, nT);
    for c = 1:4
        Fc = squeeze(cond_valid(:,:,c));          % [nFinal × nT]
        proj_context(c,:)     = d_context' * Fc;
        proj_value(c,:)       = d_value'   * Fc;
        proj_interaction(c,:) = d_inter'   * Fc;
    end

    % Store TDR results
    R.tdr_d_context        = d_context;
    R.tdr_d_value          = d_value;
    R.tdr_d_interaction    = d_inter;
    R.tdr_proj_context     = proj_context;     % [4 × nT]
    R.tdr_proj_value       = proj_value;
    R.tdr_proj_interaction = proj_interaction;
    R.tdr_cond_labels      = {'Stable-High','Stable-Low','OneDR-High','OneDR-Low'};

    % --- Within-context value axes (cross-context generalization) ------------
    % Fit d_value separately within each context using only that context's
    % two conditions, then measure their alignment via |cos(theta)|.
    % High |cos(theta)| → value coding direction is stable across contexts
    % (context-independent / generalizable).
    X_val_ctx = [ones(2,1), [1;-1]];    % intercept + value code for 2 conditions

    % Stable value axis: stable-high (cond 1) vs stable-low (cond 2)
    F_win_s = squeeze(nanmean(cond_avg_norm(valid_auc, norm_mask, 1:2), 2));  % [nFinal × 2]
    beta_s  = X_val_ctx \ F_win_s';     % [2 × nFinal]
    d_val_s = beta_s(2,:)';
    d_val_s = d_val_s / norm(d_val_s);

    % OneDR value axis: oneDR-high (cond 3) vs oneDR-low (cond 4)
    F_win_d = squeeze(nanmean(cond_avg_norm(valid_auc, norm_mask, 3:4), 2));  % [nFinal × 2]
    beta_d  = X_val_ctx \ F_win_d';
    d_val_d = beta_d(2,:)';
    d_val_d = d_val_d / norm(d_val_d);

    cos_val = abs(dot(d_val_s, d_val_d));   % abs: sign convention is arbitrary
    ang_val = acosd(cos_val);
    fprintf('  Within-context value axes: |cos(theta)| = %.3f  (angle = %.1f deg)\n', ...
        cos_val, ang_val);

    % Project all 4 conditions onto each within-context value axis
    proj_val_sax = NaN(4, nT);   % all conditions onto d_val_stable
    proj_val_dax = NaN(4, nT);   % all conditions onto d_val_oneDR
    for c = 1:4
        Fc = squeeze(cond_valid(:,:,c));          % [nFinal × nT]
        proj_val_sax(c,:) = d_val_s' * Fc;
        proj_val_dax(c,:) = d_val_d' * Fc;
    end

    R.tdr_d_val_stable    = d_val_s;
    R.tdr_d_val_oneDR     = d_val_d;
    R.tdr_cos_val_axes    = cos_val;
    R.tdr_ang_val_axes    = ang_val;
    R.tdr_proj_val_sax    = proj_val_sax;   % [4 × nT]
    R.tdr_proj_val_dax    = proj_val_dax;

    % --- Opto-group TDR: Fast vs Others --------------------------------------
    % Run the same 3-axis TDR separately for ChR2-expressing (Fast) neurons
    % and all other neurons, within each area.
    fast_mask   = valid_auc & (modutype == 1);   % Fast Excitation
    others_mask = valid_auc & (modutype ~= 1);   % Slow + Inhibition + NoEffect
    fprintf('  Opto groups (within valid_auc): Fast=%d, Others=%d\n', ...
        sum(fast_mask), sum(others_mask));

    R.tdr_fast   = computeTDR(cond_avg_norm, fast_mask,   norm_mask, nT);
    R.tdr_others = computeTDR(cond_avg_norm, others_mask, norm_mask, nT);
    fprintf('  Fast:   cos(val_s, val_d) = %.3f (%.1f deg)\n', ...
        R.tdr_fast.cos_val_axes,   R.tdr_fast.ang_val_axes);
    fprintf('  Others: cos(val_s, val_d) = %.3f (%.1f deg)\n', ...
        R.tdr_others.cos_val_axes, R.tdr_others.ang_val_axes);

    Results.(area) = R;
end

%% -------------------------------------------------------------------------
%  Save results
% -------------------------------------------------------------------------
SaveFile = fullfile(ResultsDir, 'AUC_Context_Results.mat');
save(SaveFile, 'Results', 'PeakWindow_ms', 'SmoothKernel', 'SmoothSigma_ms');
fprintf('\nResults saved to %s\n', SaveFile);

%% -------------------------------------------------------------------------
%  Plot
% -------------------------------------------------------------------------
% Layout: 2 rows (FEF, SC) × 3 cols
%   Col 1: High-value AUCPSTH
%   Col 2: Low-value  AUCPSTH
%   Col 3: Peak AUC bar (high vs low in target window)

colorHigh = [0.20 0.40 0.80];   % blue  — high value
colorLow  = [0.80 0.20 0.20];   % red   — low value
shadAlpha = 0.25;

fig = figure('Name', 'Context AUC Discrimination', ...
    'Position', [80 80 1300 700], 'Color', 'w', 'Visible', 'off');

for aIdx = 1:nAreas
    area = areas{aIdx};
    R    = Results.(area);
    t    = R.time_common;
    nN   = R.n_neurons;

    % Restrict to plot window
    tMask = t >= PlotWindow_ms(1) & t <= PlotWindow_ms(2);
    tPlot = t(tMask);

    % ---- Col 1: High-value AUCPSTH ------------------------------------------
    ax1 = subplot(nAreas, 3, (aIdx-1)*3 + 1);
    hold(ax1, 'on');

    mH = R.mean_high(tMask);
    sH = R.sem_high(tMask);
    fill(ax1, [tPlot fliplr(tPlot)], ...
        [mH + sH, fliplr(mH - sH)], colorHigh, ...
        'EdgeColor', 'none', 'FaceAlpha', shadAlpha);
    plot(ax1, tPlot, mH, 'Color', colorHigh, 'LineWidth', 2);

    plot(ax1, PlotWindow_ms,          [0.5 0.5],      'k--', 'LineWidth', 0.8);
    plot(ax1, [0 0],                  [0.35 0.75],    'k:',  'LineWidth', 0.8);
    plot(ax1, [1 1]*PeakWindow_ms(1), [0.35 0.75], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
    plot(ax1, [1 1]*PeakWindow_ms(2), [0.35 0.75], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

    xlabel(ax1, 'Time from Target Onset (ms)');
    ylabel(ax1, 'Context AUC (Stable vs OneDR)');
    title(ax1, sprintf('%s  |  High Value  (n=%d)', area, nN));
    xlim(ax1, PlotWindow_ms);
    ylim(ax1, [0.35 0.75]);
    box(ax1, 'off');

    % ---- Col 2: Low-value AUCPSTH -------------------------------------------
    ax2 = subplot(nAreas, 3, (aIdx-1)*3 + 2);
    hold(ax2, 'on');

    mL = R.mean_low(tMask);
    sL = R.sem_low(tMask);
    fill(ax2, [tPlot fliplr(tPlot)], ...
        [mL + sL, fliplr(mL - sL)], colorLow, ...
        'EdgeColor', 'none', 'FaceAlpha', shadAlpha);
    plot(ax2, tPlot, mL, 'Color', colorLow, 'LineWidth', 2);

    plot(ax2, PlotWindow_ms,          [0.5 0.5],      'k--', 'LineWidth', 0.8);
    plot(ax2, [0 0],                  [0.35 0.75],    'k:',  'LineWidth', 0.8);
    plot(ax2, [1 1]*PeakWindow_ms(1), [0.35 0.75], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
    plot(ax2, [1 1]*PeakWindow_ms(2), [0.35 0.75], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

    xlabel(ax2, 'Time from Target Onset (ms)');
    ylabel(ax2, 'Context AUC (Stable vs OneDR)');
    title(ax2, sprintf('%s  |  Low Value  (n=%d)', area, nN));
    xlim(ax2, PlotWindow_ms);
    ylim(ax2, [0.35 0.75]);
    box(ax2, 'off');

    % ---- Col 3: Peak AUC bar (high vs low, 100–200 ms) ---------------------
    ax3 = subplot(nAreas, 3, (aIdx-1)*3 + 3);
    hold(ax3, 'on');

    peakH_valid = R.peak_auc_high(R.valid_auc);
    peakL_valid = R.peak_auc_low( R.valid_auc);

    barVals = [mean(peakH_valid), mean(peakL_valid)];
    barSEMs = [std(peakH_valid)/sqrt(nN), std(peakL_valid)/sqrt(nN)];

    % Scatter individual neurons (jittered)
    jitter = 0.12;
    xH = 1 + jitter * (rand(nN, 1) - 0.5);
    xL = 2 + jitter * (rand(nN, 1) - 0.5);
    scatter(ax3, xH, peakH_valid, 18, colorHigh, 'filled', 'MarkerFaceAlpha', 0.3);
    scatter(ax3, xL, peakL_valid, 18, colorLow,  'filled', 'MarkerFaceAlpha', 0.3);

    % Bars and error bars on top
    bh = bar(ax3, 1:2, barVals, 0.45);
    bh.FaceColor    = 'flat';
    bh.CData        = [colorHigh; colorLow];
    bh.FaceAlpha    = 0.7;
    bh.EdgeColor    = 'none';
    errorbar(ax3, 1:2, barVals, barSEMs, 'k.', 'LineWidth', 1.5, 'CapSize', 6);

    set(ax3, 'XTick', 1:2, 'XTickLabel', {'High', 'Low'});
    xlabel(ax3, 'Value Condition');
    ylabel(ax3, 'Peak Context AUC  (100–200 ms)');
    title(ax3, sprintf('%s  |  Peak AUC', area));
    xlim(ax3, [0.4 2.6]);
    ylim(ax3, [0 0.5]);
    box(ax3, 'off');
end

sgtitle('Context Discriminability: Stable vs OneDR  (RF location, rectified AUC)', ...
    'FontSize', 13, 'FontWeight', 'bold');

% Save figure
saveas(fig, fullfile(FigDir, 'AUC_Context_Discrimination.fig'));
saveas(fig, fullfile(FigDir, 'AUC_Context_Discrimination.png'));
fprintf('Figures saved to %s\n', FigDir);

%% -------------------------------------------------------------------------
%  TDR Figure
%  3 columns per area:
%    Col 1 — Context axis projections (4 conditions vs time)
%    Col 2 — Value   axis projections (4 conditions vs time)
%    Col 3 — 2D trajectory (context axis × value axis)
% -------------------------------------------------------------------------
% Color / style scheme:
%   Stable = blue,  OneDR = orange
%   High   = solid, Low   = dashed
condColors = [0.20 0.40 0.80;   % stable-high
              0.20 0.40 0.80;   % stable-low
              0.85 0.40 0.10;   % oneDR-high
              0.85 0.40 0.10];  % oneDR-low
condStyles  = {'-', '--', '-', '--'};
condLabels  = {'Stable-High','Stable-Low','OneDR-High','OneDR-Low'};

figTDR = figure('Name', 'TDR: Context & Value Axes', ...
    'Position', [100 50 1200 700], 'Color', 'w', 'Visible', 'off');

for aIdx = 1:nAreas
    area  = areas{aIdx};
    R     = Results.(area);
    t     = R.time_common;
    tMask = t >= PlotWindow_ms(1) & t <= PlotWindow_ms(2);
    tPlot = t(tMask);
    nN    = R.n_neurons;

    % ---- Col 1: Context axis projection ------------------------------------
    ax_c = subplot(nAreas, 3, (aIdx-1)*3 + 1);
    hold(ax_c, 'on');
    for c = 1:4
        plot(ax_c, tPlot, R.tdr_proj_context(c, tMask), ...
            condStyles{c}, 'Color', condColors(c,:), 'LineWidth', 1.8);
    end
    plot(ax_c, [0 0],                    get(ax_c, 'YLim'), 'k:',  'LineWidth', 0.8);
    plot(ax_c, get(ax_c, 'XLim'), [0 0],                   'k--', 'LineWidth', 0.8);
    xlabel(ax_c, 'Time from Target Onset (ms)');
    ylabel(ax_c, 'Projection onto context axis (a.u.)');
    title(ax_c, sprintf('%s  |  Context axis  (n=%d)', area, nN));
    legend(ax_c, condLabels, 'Location', 'northwest', 'FontSize', 8);
    legend(ax_c, 'boxoff');
    box(ax_c, 'off');

    % ---- Col 2: Value axis projection --------------------------------------
    ax_v = subplot(nAreas, 3, (aIdx-1)*3 + 2);
    hold(ax_v, 'on');
    for c = 1:4
        plot(ax_v, tPlot, R.tdr_proj_value(c, tMask), ...
            condStyles{c}, 'Color', condColors(c,:), 'LineWidth', 1.8);
    end
    plot(ax_v, [0 0],                    get(ax_v, 'YLim'), 'k:',  'LineWidth', 0.8);
    plot(ax_v, get(ax_v, 'XLim'), [0 0],                   'k--', 'LineWidth', 0.8);
    xlabel(ax_v, 'Time from Target Onset (ms)');
    ylabel(ax_v, 'Projection onto value axis (a.u.)');
    title(ax_v, sprintf('%s  |  Value axis  (n=%d)', area, nN));
    legend(ax_v, condLabels, 'Location', 'northwest', 'FontSize', 8);
    legend(ax_v, 'boxoff');
    box(ax_v, 'off');

    % ---- Col 3: 2D trajectory (context × value) ----------------------------
    ax_2d = subplot(nAreas, 3, (aIdx-1)*3 + 3);
    hold(ax_2d, 'on');
    t0_idx = find(tPlot == 0, 1);
    for c = 1:4
        xc = R.tdr_proj_context(c, tMask);
        yv = R.tdr_proj_value(c,   tMask);
        plot(ax_2d, xc, yv, condStyles{c}, ...
            'Color', condColors(c,:), 'LineWidth', 1.8);
        if ~isempty(t0_idx)   % filled dot marks target onset
            scatter(ax_2d, xc(t0_idx), yv(t0_idx), 50, ...
                condColors(c,:), 'filled');
        end
    end
    plot(ax_2d, [0 0],                      get(ax_2d, 'YLim'), 'k--', 'LineWidth', 0.8);
    plot(ax_2d, get(ax_2d, 'XLim'), [0 0],                     'k--', 'LineWidth', 0.8);
    xlabel(ax_2d, 'Context axis (a.u.)');
    ylabel(ax_2d, 'Value axis (a.u.)');
    title(ax_2d, sprintf('%s  |  Trajectory', area));
    legend(ax_2d, condLabels, 'Location', 'best', 'FontSize', 8);
    legend(ax_2d, 'boxoff');
    axis(ax_2d, 'square');
    box(ax_2d, 'off');
end

sgtitle(figTDR, ...
    sprintf('TDR — Context & Value axes (fit on %d–%d ms, Gram-Schmidt orth)', ...
    PeakWindow_ms(1), PeakWindow_ms(2)), ...
    'FontSize', 12, 'FontWeight', 'bold');

set(figTDR, 'Visible', 'on');
saveas(figTDR, fullfile(FigDir, 'TDR_ContextValue.fig'));
saveas(figTDR, fullfile(FigDir, 'TDR_ContextValue.png'));
fprintf('TDR figure saved to %s\n', FigDir);

%% -------------------------------------------------------------------------
%  TDR Figure — All three axes (context, value, interaction), de-meaned
%
%  Layout: 2 rows (FEF, SC) × 3 cols
%    Col 1 — Context axis:     stable vs oneDR conditions separate
%    Col 2 — Value main effect: high vs low separate regardless of context
%    Col 3 — Interaction axis: context × value
%            Expected pattern: stable-low ≈ oneDR-high (both code = +1)
%                              stable-high ≈ oneDR-low  (both code = −1)
%            Separation here means value coding reverses between contexts.
% -------------------------------------------------------------------------
figTDR_3ax = figure('Name', 'TDR: Three axes (context, value, interaction)', ...
    'Position', [100 50 1400 700], 'Color', 'w', 'Visible', 'off');

axisLabels   = {'Context axis', 'Value axis (main effect)', 'Interaction axis (ctx \times val)'};
projFields   = {'tdr_proj_context', 'tdr_proj_value', 'tdr_proj_interaction'};
yLabelStrs   = {'Context projection (a.u.)', 'Value projection (a.u.)', 'Interaction projection (a.u.)'};

for aIdx = 1:nAreas
    area  = areas{aIdx};
    R     = Results.(area);
    t     = R.time_common;
    tMask = t >= PlotWindow_ms(1) & t <= PlotWindow_ms(2);
    tPlot = t(tMask);
    nN    = R.n_neurons;

    for col = 1:3
        ax = subplot(nAreas, 3, (aIdx-1)*3 + col);
        hold(ax, 'on');

        proj = R.(projFields{col})(:, tMask);   % [4 × nT_plot]
        proj_dm = proj - mean(proj, 1);          % de-mean across 4 conditions

        for c = 1:4
            plot(ax, tPlot, proj_dm(c,:), condStyles{c}, ...
                'Color', condColors(c,:), 'LineWidth', 1.8);
        end
        plot(ax, [0 0],             get(ax,'YLim'), 'k:',  'LineWidth', 0.8);
        plot(ax, get(ax,'XLim'), [0 0],             'k--', 'LineWidth', 0.8);

        xlabel(ax, 'Time from Target Onset (ms)');
        ylabel(ax, yLabelStrs{col});
        title(ax, sprintf('%s  |  %s  (n=%d)', area, axisLabels{col}, nN));
        if col == 1
            legend(ax, condLabels, 'Location', 'northwest', 'FontSize', 8);
            legend(ax, 'boxoff');
        end
        box(ax, 'off');
    end
end

sgtitle(figTDR_3ax, ...
    sprintf('TDR — Context, Value, and Interaction axes (de-meaned; fit on %d–%d ms)', ...
    PeakWindow_ms(1), PeakWindow_ms(2)), ...
    'FontSize', 12, 'FontWeight', 'bold');

set(figTDR_3ax, 'Visible', 'on');
saveas(figTDR_3ax, fullfile(FigDir, 'TDR_ThreeAxes.fig'));
saveas(figTDR_3ax, fullfile(FigDir, 'TDR_ThreeAxes.png'));
fprintf('Three-axis TDR figure saved to %s\n', FigDir);

%% -------------------------------------------------------------------------
%  TDR Figure — De-meaned projections
%  Same layout as above, but each time bin is mean-subtracted across the
%  4 conditions so that shared activity (the common offset) is removed.
%  This centers the trajectories on zero and makes the *direction* of
%  separation (not its absolute sign) the salient feature.
% -------------------------------------------------------------------------
figTDR_dm = figure('Name', 'TDR: Context & Value Axes (de-meaned)', ...
    'Position', [150 30 1200 700], 'Color', 'w', 'Visible', 'off');

for aIdx = 1:nAreas
    area  = areas{aIdx};
    R     = Results.(area);
    t     = R.time_common;
    tMask = t >= PlotWindow_ms(1) & t <= PlotWindow_ms(2);
    tPlot = t(tMask);
    nN    = R.n_neurons;

    % De-mean: subtract cross-condition mean at each time point
    pc    = R.tdr_proj_context(:, tMask);   % [4 × nT_plot]
    pv    = R.tdr_proj_value(:,   tMask);
    pc_dm = pc - mean(pc, 1);               % subtract condition mean
    pv_dm = pv - mean(pv, 1);

    % ---- Col 1: De-meaned context axis projection --------------------------
    ax_c = subplot(nAreas, 3, (aIdx-1)*3 + 1);
    hold(ax_c, 'on');
    for c = 1:4
        plot(ax_c, tPlot, pc_dm(c,:), ...
            condStyles{c}, 'Color', condColors(c,:), 'LineWidth', 1.8);
    end
    yc = get(ax_c, 'YLim');
    plot(ax_c, [0 0],              yc,     'k:',  'LineWidth', 0.8);
    plot(ax_c, get(ax_c,'XLim'),  [0 0],   'k--', 'LineWidth', 0.8);
    xlabel(ax_c, 'Time from Target Onset (ms)');
    ylabel(ax_c, 'Context projection, de-meaned (a.u.)');
    title(ax_c, sprintf('%s  |  Context axis  (n=%d)', area, nN));
    legend(ax_c, condLabels, 'Location', 'northwest', 'FontSize', 8);
    legend(ax_c, 'boxoff');
    box(ax_c, 'off');

    % ---- Col 2: De-meaned value axis projection ----------------------------
    ax_v = subplot(nAreas, 3, (aIdx-1)*3 + 2);
    hold(ax_v, 'on');
    for c = 1:4
        plot(ax_v, tPlot, pv_dm(c,:), ...
            condStyles{c}, 'Color', condColors(c,:), 'LineWidth', 1.8);
    end
    yv = get(ax_v, 'YLim');
    plot(ax_v, [0 0],              yv,     'k:',  'LineWidth', 0.8);
    plot(ax_v, get(ax_v,'XLim'),  [0 0],   'k--', 'LineWidth', 0.8);
    xlabel(ax_v, 'Time from Target Onset (ms)');
    ylabel(ax_v, 'Value projection, de-meaned (a.u.)');
    title(ax_v, sprintf('%s  |  Value axis  (n=%d)', area, nN));
    legend(ax_v, condLabels, 'Location', 'northwest', 'FontSize', 8);
    legend(ax_v, 'boxoff');
    box(ax_v, 'off');

    % ---- Col 3: 2D trajectory (de-meaned context × de-meaned value) --------
    ax_2d = subplot(nAreas, 3, (aIdx-1)*3 + 3);
    hold(ax_2d, 'on');
    t0_idx = find(tPlot == 0, 1);
    for c = 1:4
        plot(ax_2d, pc_dm(c,:), pv_dm(c,:), condStyles{c}, ...
            'Color', condColors(c,:), 'LineWidth', 1.8);
        if ~isempty(t0_idx)
            scatter(ax_2d, pc_dm(c, t0_idx), pv_dm(c, t0_idx), 50, ...
                condColors(c,:), 'filled');
        end
    end
    plot(ax_2d, [0 0],                      get(ax_2d,'YLim'), 'k--', 'LineWidth', 0.8);
    plot(ax_2d, get(ax_2d,'XLim'), [0 0],                      'k--', 'LineWidth', 0.8);
    xlabel(ax_2d, 'Context axis, de-meaned (a.u.)');
    ylabel(ax_2d, 'Value axis, de-meaned (a.u.)');
    title(ax_2d, sprintf('%s  |  Trajectory (de-meaned)', area));
    legend(ax_2d, condLabels, 'Location', 'best', 'FontSize', 8);
    legend(ax_2d, 'boxoff');
    axis(ax_2d, 'square');
    box(ax_2d, 'off');
end

sgtitle(figTDR_dm, ...
    sprintf('TDR — De-meaned projections (shared activity removed; fit on %d–%d ms)', ...
    PeakWindow_ms(1), PeakWindow_ms(2)), ...
    'FontSize', 12, 'FontWeight', 'bold');

set(figTDR_dm, 'Visible', 'on');
saveas(figTDR_dm, fullfile(FigDir, 'TDR_ContextValue_Demeaned.fig'));
saveas(figTDR_dm, fullfile(FigDir, 'TDR_ContextValue_Demeaned.png'));
fprintf('De-meaned TDR figure saved to %s\n', FigDir);

%% -------------------------------------------------------------------------
%  TDR Figure — Within-context value axes & cross-context generalization
%
%  Layout: 2 rows (FEF, SC) × 3 cols
%    Col 1 — All 4 conditions projected onto d_val_stable (de-meaned)
%             Within-context: stable-high/low should separate (by construction)
%             Cross-context:  does oneDR-high/low also separate here?
%    Col 2 — Same, but projected onto d_val_oneDR
%    Col 3 — |cos(theta)| summary bar comparing FEF vs SC
%             High |cos| → value axis is consistent across task contexts
% -------------------------------------------------------------------------
figTDR_gen = figure('Name', 'TDR: Cross-context value generalization', ...
    'Position', [200 30 1300 700], 'Color', 'w', 'Visible', 'off');

for aIdx = 1:nAreas
    area  = areas{aIdx};
    R     = Results.(area);
    t     = R.time_common;
    tMask = t >= PlotWindow_ms(1) & t <= PlotWindow_ms(2);
    tPlot = t(tMask);
    nN    = R.n_neurons;

    % De-mean across 4 conditions at each time point
    ps    = R.tdr_proj_val_sax(:, tMask);
    pd    = R.tdr_proj_val_dax(:, tMask);
    ps_dm = ps - mean(ps, 1);
    pd_dm = pd - mean(pd, 1);

    % ---- Col 1: Project onto stable value axis ------------------------------
    ax_s = subplot(nAreas, 3, (aIdx-1)*3 + 1);
    hold(ax_s, 'on');
    for c = 1:4
        plot(ax_s, tPlot, ps_dm(c,:), condStyles{c}, ...
            'Color', condColors(c,:), 'LineWidth', 1.8);
    end
    plot(ax_s, [0 0],              get(ax_s,'YLim'), 'k:',  'LineWidth', 0.8);
    plot(ax_s, get(ax_s,'XLim'),  [0 0],             'k--', 'LineWidth', 0.8);
    xlabel(ax_s, 'Time from Target Onset (ms)');
    ylabel(ax_s, 'Projection, de-meaned (a.u.)');
    title(ax_s, sprintf('%s  |  Stable value axis', area));
    legend(ax_s, condLabels, 'Location', 'northwest', 'FontSize', 8);
    legend(ax_s, 'boxoff');
    box(ax_s, 'off');

    % ---- Col 2: Project onto oneDR value axis -------------------------------
    ax_d = subplot(nAreas, 3, (aIdx-1)*3 + 2);
    hold(ax_d, 'on');
    for c = 1:4
        plot(ax_d, tPlot, pd_dm(c,:), condStyles{c}, ...
            'Color', condColors(c,:), 'LineWidth', 1.8);
    end
    plot(ax_d, [0 0],              get(ax_d,'YLim'), 'k:',  'LineWidth', 0.8);
    plot(ax_d, get(ax_d,'XLim'),  [0 0],             'k--', 'LineWidth', 0.8);
    xlabel(ax_d, 'Time from Target Onset (ms)');
    ylabel(ax_d, 'Projection, de-meaned (a.u.)');
    title(ax_d, sprintf('%s  |  OneDR value axis', area));
    legend(ax_d, condLabels, 'Location', 'northwest', 'FontSize', 8);
    legend(ax_d, 'boxoff');
    box(ax_d, 'off');
end

% ---- Col 3 (spans both rows): |cos(theta)| summary -------------------------
ax_cos = subplot(nAreas, 3, [3 6]);
hold(ax_cos, 'on');
cos_vals = [Results.FEF.tdr_cos_val_axes, Results.SC.tdr_cos_val_axes];
ang_vals = [Results.FEF.tdr_ang_val_axes, Results.SC.tdr_ang_val_axes];
bh = bar(ax_cos, 1:2, cos_vals, 0.5);
bh.FaceColor = 'flat';
bh.CData     = [0.20 0.40 0.80; 0.85 0.40 0.10];
bh.FaceAlpha = 0.8;
bh.EdgeColor = 'none';
for k = 1:2
    text(ax_cos, k, cos_vals(k) + 0.03, ...
        sprintf('%.3f\n(%.1f\\circ)', cos_vals(k), ang_vals(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end
plot(ax_cos, [0.4 2.6], [1 1], 'k--', 'LineWidth', 0.8);  % perfect alignment
set(ax_cos, 'XTick', 1:2, 'XTickLabel', {'FEF', 'SC'}, 'YLim', [0 1.2]);
xlabel(ax_cos, 'Area');
ylabel(ax_cos, '|cos(\theta)| between stable & oneDR value axes');
title(ax_cos, {'Cross-context value axis alignment'; ...
    '1 = identical axes,  0 = orthogonal'});
box(ax_cos, 'off');

sgtitle(figTDR_gen, ...
    'TDR — Within-context value axes & cross-context generalization (de-meaned)', ...
    'FontSize', 12, 'FontWeight', 'bold');

set(figTDR_gen, 'Visible', 'on');
saveas(figTDR_gen, fullfile(FigDir, 'TDR_CrossContextGeneralization.fig'));
saveas(figTDR_gen, fullfile(FigDir, 'TDR_CrossContextGeneralization.png'));
fprintf('Cross-context generalization figure saved to %s\n', FigDir);

%% -------------------------------------------------------------------------
%  TDR Figure — Fast vs Others comparison (FEF and SC)
%
%  Layout: 4 rows × 3 cols
%    Rows: FEF-Fast, FEF-Others, SC-Fast, SC-Others
%    Cols: Context axis, Value axis (main effect), Interaction axis
%  All projections are de-meaned across the 4 conditions.
%  Title of each row includes group label and n neurons.
%
%  Companion summary figure: |cos(theta)| for 4 groups (2 areas × 2 groups).
% -------------------------------------------------------------------------
grpFields  = {'tdr_fast', 'tdr_others'};
grpNames   = {'Fast (ChR2)', 'Others'};
grpColors  = {[0.85 0.30 0.10], [0.40 0.40 0.40]};   % orange-red, grey

rowLabels = {};
rowTDRs   = {};
rowClrs   = {};
for aIdx = 1:nAreas
    for gi = 1:2
        rowLabels{end+1} = sprintf('%s — %s', areas{aIdx}, grpNames{gi}); %#ok<SAGROW>
        rowTDRs{end+1}   = Results.(areas{aIdx}).(grpFields{gi});          %#ok<SAGROW>
        rowClrs{end+1}   = grpColors{gi};                                  %#ok<SAGROW>
    end
end
nRows = numel(rowLabels);

figTDR_grp = figure('Name', 'TDR: Fast vs Others', ...
    'Position', [80 20 1400 900], 'Color', 'w', 'Visible', 'off');

axTitles3 = {'Context axis', 'Value axis (main effect)', 'Interaction axis'};
projFlds3  = {'proj_context', 'proj_value', 'proj_interaction'};
yLbls3     = {'Context proj. (a.u.)', 'Value proj. (a.u.)', 'Interaction proj. (a.u.)'};

for row = 1:nRows
    T    = rowTDRs{row};
    clr  = rowClrs{row};
    t    = Results.(areas{ceil(row/2)}).time_common;
    tMask = t >= PlotWindow_ms(1) & t <= PlotWindow_ms(2);
    tPlot = t(tMask);

    for col = 1:3
        ax = subplot(nRows, 3, (row-1)*3 + col);
        hold(ax, 'on');

        if T.n_neurons >= 4
            proj    = T.(projFlds3{col})(:, tMask);
            proj_dm = proj - mean(proj, 1);
            for c = 1:4
                % Conditions: 1=stable-high, 2=stable-low, 3=oneDR-high, 4=oneDR-low
                % Darken or lighten for high vs low within the group color
                lw = 1.5;
                ls = condStyles{c};
                plot(ax, tPlot, proj_dm(c,:), ls, 'Color', condColors(c,:), 'LineWidth', lw);
            end
        else
            text(ax, 0.5, 0.5, 'n < 4', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
        end

        plot(ax, [0 0],            get(ax,'YLim'), 'k:',  'LineWidth', 0.8);
        plot(ax, get(ax,'XLim'), [0 0],            'k--', 'LineWidth', 0.8);
        box(ax, 'off');

        if row == 1
            title(ax, axTitles3{col}, 'FontSize', 10);
        end
        if col == 1
            cosStr = '';
            if T.n_neurons >= 4
                cosStr = sprintf('  |cos|=%.2f (%.0f°)', T.cos_val_axes, T.ang_val_axes);
            end
            ylabel(ax, sprintf('%s (n=%d)%s', rowLabels{row}, T.n_neurons, cosStr), ...
                'FontSize', 8);
        else
            ylabel(ax, yLbls3{col}, 'FontSize', 8);
        end
        xlabel(ax, 'Time (ms)', 'FontSize', 8);

        if row == 1 && col == 1
            legend(ax, condLabels, 'Location', 'northwest', 'FontSize', 7);
            legend(ax, 'boxoff');
        end
    end
end

sgtitle(figTDR_grp, ...
    'TDR — Fast vs Others  (de-meaned; context, value, interaction axes)', ...
    'FontSize', 12, 'FontWeight', 'bold');
set(figTDR_grp, 'Visible', 'on');
saveas(figTDR_grp, fullfile(FigDir, 'TDR_GroupComparison.fig'));
saveas(figTDR_grp, fullfile(FigDir, 'TDR_GroupComparison.png'));
fprintf('Group comparison TDR figure saved to %s\n', FigDir);

% ---- Summary: |cos(theta)| across 4 groups ---------------------------------
figCos = figure('Name', 'TDR: cos(val_s, val_d) by group', ...
    'Position', [200 200 500 420], 'Color', 'w', 'Visible', 'off');
ax_cs = axes(figCos);
hold(ax_cs, 'on');

cos4 = [Results.FEF.tdr_fast.cos_val_axes,   Results.FEF.tdr_others.cos_val_axes, ...
        Results.SC.tdr_fast.cos_val_axes,     Results.SC.tdr_others.cos_val_axes];
ang4 = [Results.FEF.tdr_fast.ang_val_axes,   Results.FEF.tdr_others.ang_val_axes, ...
        Results.SC.tdr_fast.ang_val_axes,     Results.SC.tdr_others.ang_val_axes];
n4   = [Results.FEF.tdr_fast.n_neurons,      Results.FEF.tdr_others.n_neurons, ...
        Results.SC.tdr_fast.n_neurons,        Results.SC.tdr_others.n_neurons];
xpos = [1, 1.6, 3, 3.6];
barClrs = [0.85 0.30 0.10; 0.55 0.55 0.55; 0.85 0.30 0.10; 0.55 0.55 0.55];

for k = 1:4
    bar(ax_cs, xpos(k), cos4(k), 0.45, 'FaceColor', barClrs(k,:), ...
        'FaceAlpha', 0.85, 'EdgeColor', 'none');
    text(ax_cs, xpos(k), cos4(k) + 0.03, ...
        sprintf('%.2f\n(%.0f°)\nn=%d', cos4(k), ang4(k), n4(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end
plot(ax_cs, [0.5 4.1], [1 1], 'k--', 'LineWidth', 0.8);
set(ax_cs, 'XTick', [1.3 3.3], 'XTickLabel', {'FEF', 'SC'}, ...
    'YLim', [0 1.25], 'XLim', [0.5 4.1]);
legend(ax_cs, {'Fast (ChR2)', 'Others'}, 'Location', 'northeast', 'FontSize', 9);
legend(ax_cs, 'boxoff');
ylabel(ax_cs, '|cos(\theta)| between stable & oneDR value axes');
title(ax_cs, {'Cross-context value axis alignment'; 'by opto group and area'});
box(ax_cs, 'off');

set(figCos, 'Visible', 'on');
saveas(figCos, fullfile(FigDir, 'TDR_GroupCosComparison.fig'));
saveas(figCos, fullfile(FigDir, 'TDR_GroupCosComparison.png'));
fprintf('Group cos comparison figure saved to %s\n', FigDir);

%% -------------------------------------------------------------------------
%  TDR Figure — 3D trajectories  (context × value × interaction)
%
%  Layout: 2 rows (FEF, SC) × 3 cols (All neurons, Fast, Others)
%  Each panel: 4 condition trajectories de-meaned across conditions.
%  Filled dot marks target onset (t = 0); open circle marks trial end.
%  Colors / line styles follow the existing convention:
%    blue  = stable,  orange = oneDR,  solid = high,  dashed = low
% -------------------------------------------------------------------------
grp3D_keys  = {'all', 'tdr_fast', 'tdr_others'};
grp3D_names = {'All neurons', 'Fast (ChR2)', 'Others'};

fig3D = figure('Name', 'TDR: 3D trajectories', ...
    'Position', [100 60 1350 820], 'Color', 'w', 'Visible', 'off');

for aIdx = 1:nAreas
    area  = areas{aIdx};
    R     = Results.(area);
    t     = R.time_common;
    tMask = t >= PlotWindow_ms(1) & t <= PlotWindow_ms(2);
    tPlot = t(tMask);
    t0_idx = find(tPlot == 0, 1);

    for gi = 1:3
        ax3 = subplot(nAreas, 3, (aIdx-1)*3 + gi);
        hold(ax3, 'on');

        % Retrieve projections for this group
        if gi == 1   % All neurons
            pc   = R.tdr_proj_context(:,     tMask);
            pv   = R.tdr_proj_value(:,       tMask);
            pint = R.tdr_proj_interaction(:, tMask);
            nG   = R.n_neurons;
        else
            T = R.(grp3D_keys{gi});
            if T.n_neurons < 4
                text(ax3, 0.5, 0.5, 0.5, 'n < 4', 'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
                title(ax3, sprintf('%s  |  %s', area, grp3D_names{gi}), 'FontSize', 9);
                box(ax3, 'on'); continue;
            end
            pc   = T.proj_context(:,     tMask);
            pv   = T.proj_value(:,       tMask);
            pint = T.proj_interaction(:, tMask);
            nG   = T.n_neurons;
        end

        % De-mean across 4 conditions at each time point
        pc   = pc   - mean(pc,   1);
        pv   = pv   - mean(pv,   1);
        pint = pint - mean(pint, 1);

        % Plot 4 condition trajectories
        for c = 1:4
            plot3(ax3, pc(c,:), pv(c,:), pint(c,:), ...
                condStyles{c}, 'Color', condColors(c,:), 'LineWidth', 1.6);
            if ~isempty(t0_idx)
                % Filled dot at target onset
                scatter3(ax3, pc(c,t0_idx), pv(c,t0_idx), pint(c,t0_idx), ...
                    55, condColors(c,:), 'filled');
            end
            % Open circle at trajectory end
            scatter3(ax3, pc(c,end), pv(c,end), pint(c,end), ...
                35, condColors(c,:), 'o');
        end

        xlabel(ax3, 'Context (a.u.)',     'FontSize', 8);
        ylabel(ax3, 'Value (a.u.)',        'FontSize', 8);
        zlabel(ax3, 'Interaction (a.u.)', 'FontSize', 8);
        title(ax3, sprintf('%s  |  %s  (n=%d)', area, grp3D_names{gi}, nG), 'FontSize', 9);
        grid(ax3, 'on');
        view(ax3, 35, 25);

        if aIdx == 1 && gi == 1
            legend(ax3, condLabels, 'Location', 'best', 'FontSize', 7);
            legend(ax3, 'boxoff');
        end
    end
end

sgtitle(fig3D, ...
    'TDR — 3D trajectories: context \times value \times interaction  (de-meaned, filled dot = target onset)', ...
    'FontSize', 11, 'FontWeight', 'bold');
set(fig3D, 'Visible', 'on');
saveas(fig3D, fullfile(FigDir, 'TDR_3DTrajectory.fig'));
saveas(fig3D, fullfile(FigDir, 'TDR_3DTrajectory.png'));
fprintf('3D trajectory figure saved to %s\n', FigDir);

%% =========================================================================
%  Local function: computeTDR
%  Run 3-axis TDR (context, value, interaction) for a subset of neurons.
% ==========================================================================
function T = computeTDR(cond_avg_norm, sel, norm_mask, nT)
% computeTDR  Compute TDR axes and projections for a selected neuron subset.
%
% Inputs:
%   cond_avg_norm : [N × nT_full × 4]  normalized condition-averaged PSTHs
%   sel           : [N × 1] logical    neurons to include
%   norm_mask     : [1 × nT_full] logical  key-window time bins for regression
%   nT            : scalar             total number of time bins
%
% Output:  struct T with fields:
%   n_neurons, d_context, d_value, d_inter,
%   proj_context, proj_value, proj_interaction,
%   cos_val_axes, ang_val_axes

T = struct();
T.n_neurons = sum(sel);

% Guard: need at least 4 neurons (= number of conditions)
if T.n_neurons < 4
    empties = {'d_context','d_value','d_inter', ...
               'proj_context','proj_value','proj_interaction', ...
               'cos_val_axes','ang_val_axes'};
    for k = 1:numel(empties), T.(empties{k}) = []; end
    return;
end

% Regression design matrix (orthogonal contrasts)
X = [ones(4,1), [-1;-1;1;1], [1;-1;1;-1], [-1;1;1;-1]];
F = squeeze(nanmean(cond_avg_norm(sel, norm_mask, :), 2));  % [nSel × 4]
beta = X \ F';   % [4 × nSel]

d_ctx_raw   = beta(2,:)';  d_ctx_raw   = d_ctx_raw   / norm(d_ctx_raw);
d_val_raw   = beta(3,:)';
d_inter_raw = beta(4,:)';

% Gram-Schmidt: d_value ⊥ d_context, d_inter ⊥ both
d_ctx = d_ctx_raw;
d_val = d_val_raw - (d_ctx' * d_val_raw) * d_ctx;
d_val = d_val / norm(d_val);
d_int = d_inter_raw - (d_ctx' * d_inter_raw) * d_ctx ...
                    - (d_val' * d_inter_raw) * d_val;
d_int = d_int / norm(d_int);

% Project all 4 conditions [4 × nT]
csel = cond_avg_norm(sel, :, :);
p_ctx = NaN(4, nT);  p_val = NaN(4, nT);  p_int = NaN(4, nT);
for c = 1:4
    Fc = squeeze(csel(:,:,c));
    p_ctx(c,:) = d_ctx' * Fc;
    p_val(c,:) = d_val' * Fc;
    p_int(c,:) = d_int' * Fc;
end

% Within-context value axes for cos(theta)
Xv  = [ones(2,1), [1;-1]];
Fs  = squeeze(nanmean(cond_avg_norm(sel, norm_mask, 1:2), 2));
bvs = Xv \ Fs';  dvs = bvs(2,:)';  dvs = dvs / norm(dvs);
Fd  = squeeze(nanmean(cond_avg_norm(sel, norm_mask, 3:4), 2));
bvd = Xv \ Fd';  dvd = bvd(2,:)';  dvd = dvd / norm(dvd);

T.d_context        = d_ctx;
T.d_value          = d_val;
T.d_inter          = d_int;
T.proj_context     = p_ctx;
T.proj_value       = p_val;
T.proj_interaction = p_int;
T.cos_val_axes     = abs(dot(dvs, dvd));
T.ang_val_axes     = acosd(T.cos_val_axes);
end
