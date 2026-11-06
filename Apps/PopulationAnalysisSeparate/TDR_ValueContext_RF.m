
function TDR = TDR_ValueContext_RF(psth_high_value_object_RF,psth_low_value_object_RF,psth_high_value_location_RF,psth_low_value_location_RF,time_axis,use_pseudo)

% TDR_ValueContext_RF  Targeted Dimensionality Reduction — RF trials only.
%
% Identical pipeline to TDR_ValueContext.m but restricted to into-RF
% conditions only (4 conditions instead of 8). The RF position axis is
% dropped; three predictors are retained: value, context, value×context.
%
% This version is the primary analysis for the decoder-consistent comparison
% (decoder also uses RF only). Use TDR_ValueContext.m when you want the
% full 8-condition geometry including the RF position axis.
%
% Inputs (expected as [nNeurons × nTimeBins] in workspace, OR pass directly):
%   psth_high_value_object_RF    [nNeurons × nTimeBins]
%   psth_low_value_object_RF
%   psth_high_value_location_RF
%   psth_low_value_location_RF
%
% Optional name-value:
%   'time_axis'   [1×nTimeBins]  seconds re stimulus onset (default linspace)
%   'var_thresh'  0.90           PCA variance threshold for denoising
%   'use_pseudo'  true           generate pseudo-PSTHs for testing
%
% Output:
%   TDR struct (also assigned to base workspace as 'TDR_RF')
%
% Reference: Tafazoli et al. (2026) Nature 650:164-172.
%{
%% ── Parse options ────────────────────────────────────────────────────────────
p = inputParser;
addOptional(p, 'psth_HV_obj', [], @isnumeric);
addOptional(p, 'psth_LV_obj', [], @isnumeric);
addOptional(p, 'psth_HV_loc', [], @isnumeric);
addOptional(p, 'psth_LV_loc', [], @isnumeric);
addParameter(p, 'time_axis',  [],   @isnumeric);
addParameter(p, 'var_thresh', 0.90, @isnumeric);
addParameter(p, 'use_pseudo', true, @islogical);
parse(p, varargin{:});
opts = p.Results;
%}
var_thresh = 0.9;
%% =====================================================================
%  SECTION 0 — DATA INPUT
%  Set use_pseudo=false and supply real PSTHs in the workspace.
% =====================================================================
%if ~isempty(opts.psth_HV_obj)
if use_pseudo == 0
    % Use directly passed PSTH matrices
    %{
    psth_high_value_object_RF   = opts.psth_HV_obj;
    psth_low_value_object_RF    = opts.psth_LV_obj;
    psth_high_value_location_RF = opts.psth_HV_loc;
    psth_low_value_location_RF  = opts.psth_LV_loc;
    %}
    [nNeurons, nTimeBins] = size(psth_high_value_object_RF);
    %{
    if ~isempty(opts.time_axis)
        time_axis = opts.time_axis(:)';
    else
        time_axis = linspace(-0.2, 0.8, nTimeBins);
    end
    %}
    fprintf('  nNeurons=%d  nTimeBins=%d\n\n', nNeurons, nTimeBins);
else%if opts.use_pseudo
    fprintf('=== Generating pseudo-PSTHs (RF only) for testing ===\n');
    nNeurons  = 120;
    nTimeBins = 60;
    time_axis = linspace(-0.2, 0.8, nTimeBins);
    rng(42);

    value_w       = randn(nNeurons, 1) * 1.0;
    context_w     = randn(nNeurons, 1) * 0.5;
    interaction_w = randn(nNeurons, 1) * 0.3;
    noise_level   = 0.4;

    t_peak   = 0.20;
    sigma_t  = 0.10;
    temp_prof = exp(-(time_axis - t_peak).^2 / (2*sigma_t^2));

    mk_psth = @(v, ctx) ...
        (value_w*v + context_w*ctx + interaction_w*(v*ctx)) ...
        .* temp_prof + noise_level * randn(nNeurons, nTimeBins);

    psth_high_value_object_RF   = mk_psth(+1, +1);
    psth_low_value_object_RF    = mk_psth(-1, +1);
    psth_high_value_location_RF = mk_psth(+1, -1);
    psth_low_value_location_RF  = mk_psth(-1, -1);

    fprintf('  nNeurons=%d  nTimeBins=%d\n\n', nNeurons, nTimeBins);
 %{
else
    % Load from base workspace
    vars = {'psth_high_value_object_RF', 'psth_low_value_object_RF', ...
            'psth_high_value_location_RF', 'psth_low_value_location_RF'};
    for k = 1:numel(vars)
        if ~evalin('base', sprintf('exist(''%s'',''var'')', vars{k}))
            error('Variable %s not found in workspace.', vars{k});
        end
        eval(sprintf('%s = evalin(''base'',''%s'');', vars{k}, vars{k}));
    end
    [nNeurons, nTimeBins] = size(psth_high_value_object_RF);

    if ~isempty(opts.time_axis)
        time_axis = opts.time_axis(:)';
    elseif evalin('base', 'exist(''time_axis'',''var'')')
        time_axis = evalin('base', 'time_axis');
    else
        warning('time_axis not found — using default linspace(-0.2, 0.8).');
        time_axis = linspace(-0.2, 0.8, nTimeBins);
    end
 %}
end

%% =====================================================================
%  SECTION 1 — PACK CONDITIONS
%
%  Condition index map:
%    1 = obj  HV RF
%    2 = obj  LV RF
%    3 = loc  HV RF
%    4 = loc  LV RF
%
%  Predictors (3):
%    value:        +1=high, -1=low
%    context:      +1=object, -1=location
%    val x ctx:    product
% =====================================================================
fprintf('Step 1: Packing RF-only conditions...\n');

all_psths = cat(3, ...
    psth_high_value_object_RF,   psth_low_value_object_RF, ...
    psth_high_value_location_RF, psth_low_value_location_RF);
% all_psths: [nNeurons × nTimeBins × 4]

nConditions = 4;

design_matrix = [...
%  val  ctx
   +1,  +1;   % obj HV RF
   -1,  +1;   % obj LV RF
   +1,  -1;   % loc HV RF
   -1,  -1];  % loc LV RF

predictor_names = {'Value', 'Context'};
nPredictors     = 2;

%% =====================================================================
%  SECTION 2 — Z-SCORE EACH NEURON
% =====================================================================
fprintf('Step 2: Z-scoring each neuron...\n');
%{
all_flat  = reshape(all_psths, nNeurons, nTimeBins * nConditions);
mu_n      = mean(all_flat, 2);
sigma_n   = std(all_flat, 0, 2) + 1e-8;
all_psths_z = (all_psths - mu_n) ./ sigma_n;
%}
all_psths_z=all_psths;
%% =====================================================================
%  SECTION 3 — PCA DENOISING
% =====================================================================
fprintf('Step 3: PCA denoising...\n');

X_concat = reshape(permute(all_psths_z, [3 2 1]), nConditions * nTimeBins, nNeurons);
[V_pca, ~, latent] = pca(X_concat);
var_expl = cumsum(latent) / sum(latent);
%nPCA     = find(var_expl >= opts.var_thresh, 1, 'first');
nPCA     = find(var_expl >= var_thresh, 1, 'first');
nPCA     = max(nPCA, 2);

fprintf('  Keeping %d PCs (%.1f%% variance)\n', nPCA, var_expl(nPCA)*100);

D = V_pca(:, 1:nPCA) * V_pca(:, 1:nPCA)';

%% =====================================================================
%  SECTION 4 — GLM: REGRESSION COEFFICIENTS AT EACH TIMEPOINT
% =====================================================================
fprintf('Step 4: GLM regression at each timepoint...\n');

X_design = [ones(nConditions, 1), design_matrix];  % [4 × 4]
beta_raw = zeros(nPredictors, nTimeBins, nNeurons);

for t = 1:nTimeBins
    Y_t = squeeze(all_psths_z(:, t, :))';   % [4 × nNeurons]
    B_t = X_design \ Y_t;                    % [4 × nNeurons]
    beta_raw(:, t, :) = B_t(2:end, :);       % drop intercept
end

%% =====================================================================
%  SECTION 5 — DENOISE REGRESSION VECTORS; FIND PEAK TIME
% =====================================================================
fprintf('Step 5: Denoising regression vectors...\n');

beta_denoised  = zeros(nPredictors, nTimeBins, nNeurons);
beta_max       = zeros(nNeurons, nPredictors);
t_max_idx      = zeros(1, nPredictors);
beta_norms_all = zeros(nPredictors, nTimeBins);

for v = 1:nPredictors
    beta_v     = squeeze(beta_raw(v, :, :))';   % [nNeurons × nTimeBins]
    beta_v_pca = D * beta_v;

    beta_denoised(v, :, :) = beta_v_pca';

    norms = sqrt(sum(beta_v_pca.^2, 1));
    beta_norms_all(v, :) = norms;
    [~, t_max_idx(v)] = max(norms);
    beta_max(:, v) = beta_v_pca(:, t_max_idx(v));

    fprintf('  %-25s  t_max = %.3f s\n', predictor_names{v}, time_axis(t_max_idx(v)));
end

%{
%% =====================================================================
%  SECTION 5b — REGRESSION-BASED VARIANCE ACCOUNTED FOR (VAF)
%
%  At each timepoint, the denoised fitted contribution of predictor v is:
%    Y_hat_v(t) = design_matrix(:,v) * beta_denoised(v,t,:)   [4 × nNeurons]
%  The design columns are zero-mean across conditions so Y_hat_v is already
%  mean-centered.  VAF is the fraction of total cross-condition variance
%  captured by that predictor's fitted signal:
%    VAF_v(t) = ||Y_hat_v||² / ||Y_t_centered||² × 100
%  NOTE: context tends to dominate SS_total if object vs location stimuli
%  drive different visual responses, washing out value/interaction VAF.
%  Use beta_norms_all for a cleaner per-predictor comparison.
% =====================================================================
fprintf('Step 5b: Computing regression-based VAF per predictor...\n');

vaf_reg_timewise = zeros(nPredictors, nTimeBins);

for t = 1:nTimeBins
    Y_t   = squeeze(all_psths_z(:, t, :))';          % [nConditions × nNeurons]
    Y_t_c = Y_t - mean(Y_t, 1);                      % mean-center across conditions
    SS_total = sum(Y_t_c(:).^2);
    if SS_total < 1e-12, continue; end
    for v = 1:nPredictors
        beta_vt = squeeze(beta_denoised(v, t, :))';  % [1 × nNeurons]
        Y_hat_v = design_matrix(:, v) * beta_vt;     % [nConditions × nNeurons]
        vaf_reg_timewise(v, t) = sum(Y_hat_v(:).^2) / SS_total * 100;
    end
end

vaf_reg_overall = mean(vaf_reg_timewise, 2);     % [nPredictors × 1]
vaf_reg_peak    = max(vaf_reg_timewise,  [], 2); % [nPredictors × 1]

fprintf('  Regression VAF per predictor (mean | peak):\n');
for v = 1:nPredictors
    fprintf('    %-25s  mean = %5.2f%%  |  peak = %5.2f%%\n', ...
        predictor_names{v}, vaf_reg_overall(v), vaf_reg_peak(v));
end
%}
%{
%% =====================================================================
%  SECTION 5b — PERMUTATION TEST: SIGNIFICANCE OF EACH PREDICTOR
%  Null: shuffle neuron identities (columns of Y_t) at each timepoint.
%  This destroys the population-level structure linking neural responses
%  to the design, while preserving each neuron's marginal distribution
%  and keeping the design matrix intact. With only 4 conditions, shuffling
%  condition labels gives only 24 unique permutations and keeps val*ctx
%  coupled — neuron shuffling avoids both problems.
%    p_peak(v)    — fraction of perms where peak null norm ≥ observed peak
%    sig_timewise — logical [nPredictors × nTimeBins], p < 0.05 at each bin
% =====================================================================
nPerms = 1000;
fprintf('Step 5b: Permutation test (%d permutations, neuron shuffle)...\n', nPerms);

beta_norms_null = zeros(nPredictors, nTimeBins, nPerms);

for perm = 1:nPerms
    beta_perm = zeros(nPredictors, nTimeBins, nNeurons);
    for t = 1:nTimeBins
        Y_t      = squeeze(all_psths_z(:, t, :))';      % [4 × nNeurons]
        Y_perm   = Y_t(:, randperm(nNeurons));           % shuffle neuron columns
        B_t      = X_design \ Y_perm;                   % [4 × nNeurons]
        beta_perm(:, t, :) = B_t(2:end, :);
    end

    for v = 1:nPredictors
        bv     = squeeze(beta_perm(v, :, :))';           % [nNeurons × nTimeBins]
        bv_pca = D * bv;
        beta_norms_null(v, :, perm) = sqrt(sum(bv_pca.^2, 1));
    end
end

% Timewise threshold: 95th percentile of null distribution at each bin.
% Reshape to 2D [P*T × nPerms] so prctile operates reliably across perms,
% then reshape back — avoids MATLAB version inconsistencies with 3D prctile.
null_flat      = reshape(beta_norms_null, nPredictors*nTimeBins, nPerms);
thresh_flat    = prctile(null_flat, 95, 2);                          % [P*T × 1]
null_thresh_95 = reshape(thresh_flat, nPredictors, nTimeBins);       % [nPredictors × nTimeBins]
sig_timewise   = beta_norms_all > null_thresh_95;                    % logical [nPredictors × nTimeBins]

% Single p-value per predictor: fraction of perms where peak null ≥ observed peak
null_peak = squeeze(max(beta_norms_null, [], 2));    % [nPredictors × nPerms]
obs_peak  = max(beta_norms_all, [], 2);              % [nPredictors × 1]
p_peak    = mean(null_peak >= obs_peak, 2);          % [nPredictors × 1]

fprintf('  Peak-time p-values:\n');
for v = 1:nPredictors
    sig_str = '';
    if p_peak(v) < 0.05, sig_str = '*'; end
    fprintf('    %-25s  p = %.4f  %s\n', predictor_names{v}, p_peak(v), sig_str);
end
%}
%% =====================================================================
%  SECTION 5c — CROSS-CONTEXT ANGLE BETWEEN VALUE AXES
%
%  Computed here (before QR) because this uses raw context-specific
%  regression vectors — measuring the natural geometry between the two
%  contexts' value directions before orthogonalization distorts them.
% =====================================================================
fprintf('Step 5c: Cross-context angle between value axes...\n');

% Context-specific value regression vector at value peak time, denoised.
% t_max_idx(1) = value predictor peak time, already found in Section 5.
X_obj_HV_z = (psth_high_value_object_RF   - mu_n) ./ sigma_n;
X_obj_LV_z = (psth_low_value_object_RF    - mu_n) ./ sigma_n;
X_loc_HV_z = (psth_high_value_location_RF - mu_n) ./ sigma_n;
X_loc_LV_z = (psth_low_value_location_RF  - mu_n) ./ sigma_n;

t_pk_val = t_max_idx(1);   % value-predictor peak time index
beta_val_obj = D * (X_obj_HV_z(:, t_pk_val) - X_obj_LV_z(:, t_pk_val));
beta_val_loc = D * (X_loc_HV_z(:, t_pk_val) - X_loc_LV_z(:, t_pk_val));

cos_theta = (beta_val_obj' * beta_val_loc) / ...
            (norm(beta_val_obj) * norm(beta_val_loc));
cos_theta = min(max(cos_theta, -1), 1);
angle_deg = acosd(abs(cos_theta));

fprintf('  Angle between value axes (obj vs loc): %.2f degrees\n', angle_deg);
fprintf('  cos(θ)=%.4f  |  Shared variance=%.1f%%\n', cos_theta, cos_theta^2*100);

%% =====================================================================
%  SECTION 6 — QR DECOMPOSITION → ORTHOGONAL TDR AXES
% =====================================================================
fprintf('Step 6: QR orthogonalization...\n');

% Context anchored first (col 1) so the value axis is orthogonalized against
% context — this asks "is there value coding independent of which task?"
[Q_tdr, ~] = qr(beta_max(:, [2 1]), 0);   % [nNeurons × 2]
TDR_axes   = Q_tdr;
% Column 1: context axis (anchored — untouched by QR)
% Column 2: value axis (orthogonalized against context)

% Align sign: TDR col 1 (context) aligns with beta_max col 2 (context predictor)
%             TDR col 2 (value)   aligns with beta_max col 1 (value predictor)
axis_to_predictor = [2, 1];
for k = 1:nPredictors
    v = axis_to_predictor(k);
    if dot(TDR_axes(:,k), beta_max(:,v)) < 0
        TDR_axes(:,k) = -TDR_axes(:,k);
    end
end

fprintf('  3 orthogonal axes in %d-dimensional space\n', nNeurons);

%% =====================================================================
%  SECTION 7 — PROJECT ALL CONDITION PSTHs ONTO TDR AXES
% =====================================================================
fprintf('Step 7: Projecting PSTHs onto TDR axes...\n');

% projections: [nConditions × nTimeBins × nPredictors]
projections = zeros(nConditions, nTimeBins, nPredictors);
for c = 1:nConditions
    X_c = squeeze(all_psths_z(:, :, c));    % [nNeurons × nTimeBins]
    for v = 1:nPredictors
        projections(c, :, v) = TDR_axes(:, v)' * X_c;
    end
end

P = @(cond, ax) squeeze(projections(cond, :, ax));

%% =====================================================================
%  SECTION 8 — VALUE ENCODING DISTANCES
%  |projection_HV - projection_LV| along the value axis, per context
% =====================================================================
fprintf('Step 8: Computing value encoding distances...\n');

% Object context (cond 1=objHV, 2=objLV) — value axis is col 2
val_dist_obj = abs(P(1,2) - P(2,2));

% Location context (cond 3=locHV, 4=locLV)
val_dist_loc = abs(P(3,2) - P(4,2));

%{
%% =====================================================================
%  SECTION 9b — VARIANCE ACCOUNTED FOR (VAF) BY EACH TDR AXIS
%
%  At each timepoint, mean-center activity across the 4 conditions (per
%  neuron), then compute what fraction of the total cross-condition
%  variance is captured by projecting onto each TDR axis.
%
%  vaf_timewise(v,t): VAF (%) for predictor v at timepoint t
%  vaf_overall(v):    mean across all timebins
%  vaf_peak(v):       peak across all timebins
% =====================================================================
fprintf('Step 9b: Computing variance accounted for (VAF) per predictor...\n');

vaf_timewise = zeros(nPredictors, nTimeBins);

for t = 1:nTimeBins
    Y_t   = squeeze(all_psths_z(:, t, :));   % [nNeurons × nConditions]
    Y_t_c = Y_t - mean(Y_t, 2);              % mean-center across conditions
    total_var_t = sum(Y_t_c(:).^2);
    if total_var_t < 1e-12, continue; end
    for v = 1:nPredictors
        proj_v = TDR_axes(:, v)' * Y_t_c;   % [1 × nConditions]
        vaf_timewise(v, t) = sum(proj_v.^2) / total_var_t * 100;
    end
end

vaf_overall = mean(vaf_timewise, 2);     % [nPredictors × 1]  mean across time
vaf_peak    = max(vaf_timewise,  [], 2); % [nPredictors × 1]  peak across time

fprintf('  VAF per predictor (mean across time | peak):\n');
for v = 1:nPredictors
    fprintf('    %-25s  mean = %5.2f%%  |  peak = %5.2f%%\n', ...
        predictor_names{v}, vaf_overall(v), vaf_peak(v));
end
%}

%% =====================================================================
%  SECTION 10 — PACK RESULTS
% =====================================================================
TDR = struct();
TDR.time_axis           = time_axis;
TDR.predictor_names     = predictor_names;
TDR.axes                = TDR_axes;           % [nNeurons × 3]
TDR.beta_max            = beta_max;
TDR.beta_denoised       = beta_denoised;
TDR.beta_norms          = beta_norms_all;
TDR.t_max_idx           = t_max_idx;
TDR.projections         = projections;        % [4 × nTimeBins × 3]
TDR.val_dist_obj        = val_dist_obj;
TDR.val_dist_loc        = val_dist_loc;
TDR.beta_val_obj        = beta_val_obj;
TDR.beta_val_loc        = beta_val_loc;
TDR.cos_theta           = cos_theta;
TDR.angle_value_axes_deg = angle_deg;
TDR.nPCA                = nPCA;
TDR.var_explained       = var_expl;
TDR.D                   = D;
TDR.nNeurons            = nNeurons;
TDR.nTimeBins           = nTimeBins;
% VAF fields commented out — re-enable Section 5b and 9b first
% TDR.vaf_reg_timewise    = vaf_reg_timewise;
% TDR.vaf_reg_overall     = vaf_reg_overall;
% TDR.vaf_reg_peak        = vaf_reg_peak;
% TDR.vaf_timewise        = vaf_timewise;
% TDR.vaf_overall         = vaf_overall;
% TDR.vaf_peak            = vaf_peak;
%TDR.p_peak              = p_peak;          % [nPredictors × 1] permutation p-value at peak
%TDR.sig_timewise        = sig_timewise;    % [nPredictors × nTimeBins] p<0.05 mask
%TDR.null_thresh_95      = null_thresh_95;  % [nPredictors × nTimeBins] 95th pct null norm
%TDR.nPerms              = nPerms;

assignin('base', 'TDR_RF', TDR);

%% =====================================================================
%  SECTION 11 — VISUALISATION
% =====================================================================
fprintf('\nStep 10: Plotting...\n');

c_HV   = [0.85 0.15 0.15];
c_LV   = [0.15 0.35 0.85];
c_obj  = [0.15 0.65 0.30];
c_loc  = [0.80 0.45 0.10];
lw     = 2;
t_win  = time_axis >= 0.1 & time_axis <= 0.5;

figure('Name','TDR (RF only): Value × Context', ...
       'Position',[40 40 1400 800], 'Color','w');

%% Panel 1: Regression coefficient norms
ax1 = subplot(2,4,1);
hold on;
clrs = lines(nPredictors);
for v = 1:nPredictors
    %{
    % Shade significant periods along the x-axis
    sig_t = sig_timewise(v,:);
    if any(sig_t)
        ymax_shade = max(beta_norms_all(:)) * 1.15;
        fill([time_axis(sig_t), fliplr(time_axis(sig_t))], ...
            [zeros(1,sum(sig_t)), ymax_shade*ones(1,sum(sig_t))], ...
            clrs(v,:), 'EdgeColor','none', 'FaceAlpha', 0.12, ...
            'HandleVisibility','off');
    end
    % 95th percentile null threshold
    plot(time_axis, null_thresh_95(v,:), '--', ...
        'Color', clrs(v,:), 'LineWidth', 0.8, 'HandleVisibility','off');
    %}
    % Observed norm
    %{
    plot(time_axis, beta_norms_all(v,:), '-', ...
        'Color', clrs(v,:), 'LineWidth', lw, ...
        'DisplayName', sprintf('%s (p=%.3f)', strrep(predictor_names{v},'_',' '), p_peak(v)));
    %}
     plot(time_axis, beta_norms_all(v,:), '-', ...
        'Color', clrs(v,:), 'LineWidth', lw, ...
        'DisplayName', sprintf('%s ', strrep(predictor_names{v},'_',' ')));

    xline(time_axis(t_max_idx(v)), 'LineStyle', '--', 'Color', [clrs(v,:), 0.5], ...
        'HandleVisibility','off');
end
xline(0,'LineStyle',':','Color','k','LineWidth',1,'HandleVisibility','off');
xlabel('Time (s)'); ylabel('‖β‖ (denoised)');
title('Regression coeff norms (dashed = 95% null)');
legend('Location','northwest','FontSize',7); box off;

%% Panel 2: Value axis — object context
ax2 = subplot(2,4,2);
hold on;
plot(time_axis, P(1,2), '-', 'Color', c_HV, 'LineWidth', lw, 'DisplayName','HV');
plot(time_axis, P(2,2), '-', 'Color', c_LV,  'LineWidth', lw, 'DisplayName','LV');
xline(0,'LineStyle',':','Color','k','LineWidth',1,'HandleVisibility','off');
xlabel('Time (s)'); ylabel('Projection (a.u.)');
title('Value axis ⊥ context — Object ctx');
legend('Location','northwest','FontSize',7); box off;

%% Panel 3: Value axis — location context
ax3 = subplot(2,4,3);
hold on;
plot(time_axis, P(3,2), '-', 'Color', c_HV, 'LineWidth', lw, 'DisplayName','HV');
plot(time_axis, P(4,2), '-', 'Color', c_LV,  'LineWidth', lw, 'DisplayName','LV');
xline(0,'LineStyle',':','Color','k','LineWidth',1,'HandleVisibility','off');
xlabel('Time (s)'); ylabel('Projection (a.u.)');
title('Value axis ⊥ context — Location ctx');
legend('Location','northwest','FontSize',7); box off;

%% Panel 4: Value encoding distance over time
ax4 = subplot(2,4,4);
hold on;
plot(time_axis, val_dist_obj, '-', 'Color', c_obj, 'LineWidth', lw, ...
    'DisplayName','Object ctx');
plot(time_axis, val_dist_loc, '-', 'Color', c_loc, 'LineWidth', lw, ...
    'DisplayName','Location ctx');
xline(0,'LineStyle',':','Color','k','LineWidth',1,'HandleVisibility','off');
patch([time_axis(t_win), fliplr(time_axis(t_win))], ...
    [zeros(1,sum(t_win)), ones(1,sum(t_win))*max([val_dist_obj,val_dist_loc])*1.2], ...
    [0.9 0.9 0.9],'EdgeColor','none','FaceAlpha',0.3,'HandleVisibility','off');
xlabel('Time (s)'); ylabel('|HV − LV|');
title('Value encoding distance');
legend('Location','northwest','FontSize',7); box off;
ylim([0, max([val_dist_obj, val_dist_loc])*1.15]);

%% Panel 5: Context axis projections
ax5 = subplot(2,4,5);
hold on;
plot(time_axis, P(1,1), '-',  'Color', c_obj,     'LineWidth', lw, 'DisplayName','Obj HV');
plot(time_axis, P(2,1), '--', 'Color', c_obj,     'LineWidth', lw, 'DisplayName','Obj LV');
plot(time_axis, P(3,1), '-',  'Color', c_loc,     'LineWidth', lw, 'DisplayName','Loc HV');
plot(time_axis, P(4,1), '--', 'Color', c_loc,     'LineWidth', lw, 'DisplayName','Loc LV');
xline(0,'LineStyle',':','Color','k','LineWidth',1,'HandleVisibility','off');
xlabel('Time (s)'); ylabel('Projection (a.u.)');
title('Context axis (anchored)');
legend('Location','northwest','FontSize',7); box off;

% Panel 6: Interaction axis — removed (interaction term dropped)

%% Panel 7: 2D trajectory — value axis vs context axis
ax7 = subplot(2,4,7);
t_traj = time_axis >= 0 & time_axis <= 0.5;
hold on;
cond_info = {1, c_HV, '-', 'Obj HV'; 2, c_LV, '-', 'Obj LV'; ...
             3, c_HV*0.6, '--', 'Loc HV'; 4, c_LV*0.6, '--', 'Loc LV'};
for k = 1:4
    ci  = cond_info{k,1};
    col = cond_info{k,2};
    sty = cond_info{k,3};
    nm  = cond_info{k,4};
    xd  = squeeze(projections(ci, t_traj, 2));   % value axis (col 2, ⊥ context)
    yd  = squeeze(projections(ci, t_traj, 1));   % context axis (col 1, anchored)
    plot(xd, yd, sty, 'Color', col, 'LineWidth', lw, 'DisplayName', nm);
    plot(xd(1),   yd(1),   'o', 'Color', col, 'MarkerFaceColor','w',  'MarkerSize',7,'HandleVisibility','off');
    plot(xd(end), yd(end), 's', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize',7,'HandleVisibility','off');
end
xlabel('Value axis (\perp context)'); ylabel('Context axis (anchored)');
title('2D trajectory: value vs context');
legend('Location','best','FontSize',7);
axis equal; box off;

%% Panel 8: Cross-context angle summary
ax8 = subplot(2,4,8);
ang_rad = acos(abs(cos_theta));
hold on;
quiver(0,0,1,0,0,'Color',c_obj,'LineWidth',3,'MaxHeadSize',0.5,'DisplayName','Object ctx');
quiver(0,0,cos(ang_rad),sin(ang_rad),0,'Color',c_loc,'LineWidth',3,'MaxHeadSize',0.5,'DisplayName','Location ctx');
arc_t = linspace(0, ang_rad, 50);
plot(0.3*cos(arc_t), 0.3*sin(arc_t), 'k-', 'LineWidth',1.5,'HandleVisibility','off');
text(0.35, 0.12, sprintf('θ=%.1f°', angle_deg), 'FontSize', 11);
xlim([-0.2 1.3]); ylim([-0.3 1.3]);
legend('Location','northwest','FontSize',7);
title(sprintf('Value axis angle: %.1f°\nShared variance: %.0f%%', ...
    angle_deg, cos_theta^2*100));
axis off;

sgtitle(sprintf('TDR (RF only) — Value × Context  (N=%d, nPCA=%d)', ...
    nNeurons, nPCA), 'FontSize',13,'FontWeight','bold');

%% =====================================================================
%  SUMMARY
% =====================================================================
fprintf('\n======== TDR RF-only Summary ========\n');
fprintf('  Neurons:  %d  |  TimeBins: %d  |  PCs: %d (%.1f%%)\n', ...
    nNeurons, nTimeBins, nPCA, var_expl(nPCA)*100);
fprintf('\n  Peak times of regression norms:\n');
for v = 1:nPredictors
    fprintf('    %-25s  t = %.3f s\n', predictor_names{v}, time_axis(t_max_idx(v)));
end
fprintf('\n  Mean value encoding distance (100–500 ms):\n');
fprintf('    Object context:   %.4f\n', mean(val_dist_obj(t_win)));
fprintf('    Location context: %.4f\n', mean(val_dist_loc(t_win)));
fprintf('\n  Cross-context value axis angle:\n');
fprintf('    θ = %.2f deg  |  cos(θ) = %.4f\n', angle_deg, cos_theta);
fprintf('    Shared variance between value axes: %.1f%%\n', cos_theta^2*100);
% VAF summary commented out
% fprintf('\n  Regression VAF per predictor (mean | peak):\n');
% for v = 1:nPredictors
%     fprintf('    %-25s  %9.2f  %9.2f\n', predictor_names{v}, vaf_reg_overall(v), vaf_reg_peak(v));
% end
% fprintf('\n  TDR axis VAF (mean | peak):\n');
% for v = 1:nPredictors
%     fprintf('    %-25s  %9.2f  %9.2f\n', predictor_names{v}, vaf_overall(v), vaf_peak(v));
% end
fprintf('======================================\n');

end  % TDR_ValueContext_RF