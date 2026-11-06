function TDR_ValueContext(psth_high_value_object_RF,psth_low_value_object_RF,psth_high_value_object_nonRF,psth_low_value_object_nonRF,psth_high_value_location_RF,psth_low_value_location_RF,psth_high_value_location_nonRF,psth_low_value_location_nonRF )
% TDR_ValueContext — Targeted Dimensionality Reduction for value × context analysis
%
% Identifies orthogonal neural subspaces encoding value (high vs low),
% task context (object vs location), RF position (RF vs nonRF), and their
% interaction, following Tafazoli et al. (2026) Nature 650:164-172.
%
% Pipeline:
%   1. Z-score each neuron across all conditions and time
%   2. PCA denoising of population response matrix
%   3. GLM regression at each timepoint → coefficient vectors β per predictor
%   4. Denoise β vectors, find peak time t_max per predictor
%   5. QR decompose stacked β_max → orthogonal encoding axes
%   6. Project all condition PSTHs onto TDR axes
%   7. Compute value encoding distances and cross-context axis angle
%   8. Visualize
%
% Inputs (expected as [nNeurons × nTimeBins] matrices in workspace):
%   psth_high_value_object_RF, psth_low_value_object_RF
%   psth_high_value_object_nonRF, psth_low_value_object_nonRF
%   psth_high_value_location_RF, psth_low_value_location_RF
%   psth_high_value_location_nonRF, psth_low_value_location_nonRF
%
% Output:
%   TDR struct (saved to workspace) containing axes, projections, distances
%
% Reference: Tafazoli S et al. (2026) Building compositional tasks with
%            shared neural subspaces. Nature 650:164-172.
%            https://doi.org/10.1038/s41586-025-09805-2

%% =====================================================================
%  SECTION 0 — PSEUDO-PSTH GENERATION FOR TESTING
%  Replace this entire section with real data when available.
%  Expected format: each psth_X is [nNeurons × nTimeBins]
% =====================================================================

USE_PSEUDO_DATA = false;  % set to false to use real workspace PSTHs

if USE_PSEUDO_DATA
    fprintf('=== Generating pseudo-PSTHs for testing ===\n');

    nNeurons  = 120;   % pseudopopulation size
    nTimeBins = 60;    % time bins
    time_axis = linspace(-0.2, 0.8, nTimeBins);  % seconds re stimulus onset

    rng(42);

    % Ground-truth per-neuron weight vectors for each predictor
    value_w       = randn(nNeurons, 1) * 1.0;   % strong value encoding
    context_w     = randn(nNeurons, 1) * 0.6;   % moderate context encoding
    RF_w          = randn(nNeurons, 1) * 0.8;   % RF position encoding
    interaction_w = randn(nNeurons, 1) * 0.3;   % value × context (weaker)

    % Temporal profile: Gaussian bump peaking at ~200 ms post-stimulus
    t_peak      = 0.20;   % seconds
    sigma_t     = 0.10;
    temp_prof   = exp(-(time_axis - t_peak).^2 / (2*sigma_t^2));  % [1×nTimeBins]

    noise_level = 0.4;

    % Condition coding: value (+1=high, -1=low), context (+1=obj, -1=loc), RF (+1=RF, -1=nonRF)
    mk_psth = @(v, ctx, rf) ...
        (value_w*v + context_w*ctx + RF_w*rf + interaction_w*(v*ctx)) ...
        .* temp_prof + noise_level * randn(nNeurons, nTimeBins);

    psth_high_value_object_RF      = mk_psth(+1, +1, +1);
    psth_low_value_object_RF       = mk_psth(-1, +1, +1);
    psth_high_value_object_nonRF   = mk_psth(+1, +1, -1);
    psth_low_value_object_nonRF    = mk_psth(-1, +1, -1);
    psth_high_value_location_RF    = mk_psth(+1, -1, +1);
    psth_low_value_location_RF     = mk_psth(-1, -1, +1);
    psth_high_value_location_nonRF = mk_psth(+1, -1, -1);
    psth_low_value_location_nonRF  = mk_psth(-1, -1, -1);

    fprintf('  nNeurons=%d  nTimeBins=%d\n\n', nNeurons, nTimeBins);
else
    %{
    % ---- Retrieve real data from base workspace ----
    vars = {'psth_high_value_object_RF',    'psth_low_value_object_RF', ...
            'psth_high_value_object_nonRF', 'psth_low_value_object_nonRF', ...
            'psth_high_value_location_RF',  'psth_low_value_location_RF', ...
            'psth_high_value_location_nonRF','psth_low_value_location_nonRF'};
    for k = 1:numel(vars)
        if ~evalin('base', sprintf('exist(''%s'',''var'')', vars{k}))
            error('Variable %s not found in workspace.', vars{k});
        end
        eval(sprintf('%s = evalin(''base'',''%s'');', vars{k}, vars{k}));
    end
    %}
    [nNeurons, nTimeBins] = size(psth_high_value_object_RF);
    % User must supply time_axis separately, or define it here:
    % time_axis = linspace(-0.2, 0.8, nTimeBins);
    if ~evalin('base','exist(''time_axis'',''var'')')
        warning('time_axis not found in workspace; using default linspace(-0.2,0.8).');
        time_axis = linspace(-0.2, 0.8, nTimeBins);
    else
        time_axis = evalin('base','time_axis');
    end
end

%% =====================================================================
%  SECTION 1 — PACK CONDITIONS INTO A SINGLE ARRAY
% =====================================================================
%
% Condition index map (used throughout):
%   1 = obj  HV RF     5 = loc HV RF
%   2 = obj  LV RF     6 = loc LV RF
%   3 = obj  HV nonRF  7 = loc HV nonRF
%   4 = obj  LV nonRF  8 = loc LV nonRF
%
% all_psths: [nNeurons × nTimeBins × 8]

fprintf('Step 1: Packing conditions...\n');

all_psths = cat(3, ...
    psth_high_value_object_RF,       psth_low_value_object_RF, ...
    psth_high_value_object_nonRF,    psth_low_value_object_nonRF, ...
    psth_high_value_location_RF,     psth_low_value_location_RF, ...
    psth_high_value_location_nonRF,  psth_low_value_location_nonRF);

nConditions = size(all_psths, 3);   % 8

% Design matrix rows = conditions, columns = [value, context, RF, val×ctx]
%   value:    +1 = high,    -1 = low
%   context:  +1 = object,  -1 = location
%   RF:       +1 = in-RF,   -1 = non-RF
%   val×ctx:  product of value and context columns
design_matrix = [...
  %  val  ctx   RF  val*ctx
    +1,  +1,  +1,   +1;   % obj  HV RF
    -1,  +1,  +1,   -1;   % obj  LV RF
    +1,  +1,  -1,   +1;   % obj  HV nonRF
    -1,  +1,  -1,   -1;   % obj  LV nonRF
    +1,  -1,  +1,   -1;   % loc  HV RF
    -1,  -1,  +1,   +1;   % loc  LV RF
    +1,  -1,  -1,   -1;   % loc  HV nonRF
    -1,  -1,  -1,   +1];  % loc  LV nonRF

predictor_names = {'Value', 'Context', 'RF_position', 'Value_x_Context'};
nPredictors     = size(design_matrix, 2);   % 4

%% =====================================================================
%  SECTION 2 — Z-SCORE EACH NEURON
%  Following Tafazoli et al.: z-score each unit by subtracting the mean
%  and dividing by s.d. computed across all conditions and timepoints.
% =====================================================================

fprintf('Step 2: Z-scoring each neuron...\n');

all_flat  = reshape(all_psths, nNeurons, nTimeBins * nConditions);
mu_n      = mean(all_flat, 2);           % [nNeurons × 1]
sigma_n   = std(all_flat, 0, 2) + 1e-8; % [nNeurons × 1]

all_psths_z = (all_psths - mu_n) ./ sigma_n;  % [nNeurons × nTimeBins × nConditions]

%% =====================================================================
%  SECTION 3 — PCA DENOISING
%  Build population response matrix X: [(nConditions·nTimeBins) × nNeurons]
%  Keep top N_PCA components explaining ≥90% of variance.
%  Denoising matrix D = V_PCA * V_PCA'  [nNeurons × nNeurons]
% =====================================================================

fprintf('Step 3: PCA denoising...\n');

% Reshape: stack all (condition, time) rows → [nConditions*nTimeBins × nNeurons]
X_concat = reshape(permute(all_psths_z, [3 2 1]), nConditions * nTimeBins, nNeurons);

[V_pca, ~, latent] = pca(X_concat);           % V_pca: [nNeurons × nNeurons]
var_expl  = cumsum(latent) / sum(latent);
nPCA      = find(var_expl >= 0.90, 1, 'first');
nPCA      = max(nPCA, 2);                      % always keep ≥2 components

fprintf('  Keeping %d PCs (%.1f%% variance explained)\n', nPCA, var_expl(nPCA)*100);

D = V_pca(:, 1:nPCA) * V_pca(:, 1:nPCA)';    % denoising projector [nNeurons × nNeurons]

%% =====================================================================
%  SECTION 4 — GLM: REGRESSION COEFFICIENTS AT EACH TIMEPOINT
%  For each time bin t: fit y = [1, design] * B  via OLS
%  β_{v,t}: regression coefficient of predictor v for each neuron at time t
%  beta_raw: [nPredictors × nTimeBins × nNeurons]
% =====================================================================

fprintf('Step 4: GLM regression at each timepoint...\n');

X_design = [ones(nConditions, 1), design_matrix];  % [8 × 5] (intercept + 4 preds)
beta_raw  = zeros(nPredictors, nTimeBins, nNeurons);

for t = 1:nTimeBins
    Y_t  = squeeze(all_psths_z(:, t, :))';      % [nConditions × nNeurons]
    B_t  = X_design \ Y_t;                       % [5 × nNeurons], OLS
    beta_raw(:, t, :) = B_t(2:end, :);           % drop intercept row
end

%% =====================================================================
%  SECTION 5 — DENOISE REGRESSION VECTORS; FIND PEAK TIME
%  β_{v,t}^PCA = D · β_{v,t}
%  t_max_v = argmax_t ‖β_{v,t}^PCA‖
%  β_v^max = β_{v,t_max_v}^PCA   (time-independent, denoised vector)
% =====================================================================

fprintf('Step 5: Denoising regression vectors and finding peak times...\n');

beta_denoised = zeros(nPredictors, nTimeBins, nNeurons);
beta_max      = zeros(nNeurons, nPredictors);   % [nNeurons × nPredictors]
t_max_idx     = zeros(1, nPredictors);
beta_norms_all = zeros(nPredictors, nTimeBins);

for v = 1:nPredictors
    beta_v     = squeeze(beta_raw(v, :, :))';        % [nNeurons × nTimeBins]
    beta_v_pca = D * beta_v;                          % denoised [nNeurons × nTimeBins]

    beta_denoised(v, :, :) = beta_v_pca';

    norms = sqrt(sum(beta_v_pca.^2, 1));              % [1 × nTimeBins]
    beta_norms_all(v, :) = norms;
    [~, t_max_idx(v)] = max(norms);
    beta_max(:, v) = beta_v_pca(:, t_max_idx(v));    % peak regression vector

    fprintf('  %-25s  t_max = %.3f s\n', predictor_names{v}, time_axis(t_max_idx(v)));
end

%% =====================================================================
%  SECTION 6 — QR DECOMPOSITION → ORTHOGONAL TDR AXES
%  B_max = [β_value | β_context | β_RF | β_interaction]  [nNeurons × 4]
%  B_max = Q · R  (economy QR)
%  Columns of Q are orthonormal encoding axes in the order of predictors.
%  The first column is most aligned with the value regression vector, etc.
% =====================================================================

fprintf('Step 6: QR orthogonalization → TDR axes...\n');

[Q_tdr, ~] = qr(beta_max, 0);   % Q_tdr: [nNeurons × nPredictors]
TDR_axes   = Q_tdr;             % columns: [value_ax | context_ax | RF_ax | interaction_ax]

fprintf('  %d orthogonal axes computed in %d-dimensional space\n', nPredictors, nNeurons);

%% =====================================================================
%  SECTION 7 — PROJECT ALL CONDITION PSTHs ONTO TDR AXES
%  p_{v,c}(t) = TDR_axes[:,v]' · X_c   where X_c is [nNeurons × nTimeBins]
%  projections: [nConditions × nTimeBins × nPredictors]
% =====================================================================

fprintf('Step 7: Projecting PSTHs onto TDR axes...\n');

projections = zeros(nConditions, nTimeBins, nPredictors);

for c = 1:nConditions
    X_c = squeeze(all_psths_z(:, :, c));     % [nNeurons × nTimeBins]
    for v = 1:nPredictors
        projections(c, :, v) = TDR_axes(:, v)' * X_c;  % [1 × nTimeBins]
    end
end

% Convenience aliases for the value axis (axis index 1)
P = @(cond, ax) squeeze(projections(cond, :, ax));  % [1 × nTimeBins]

%% =====================================================================
%  SECTION 8 — VALUE ENCODING DISTANCES
%  Measures separability of high vs low value along the value axis,
%  analogous to the "colour encoding distance" in Tafazoli et al.
%
%  value_enc_dist_X(t) = mean over RF positions of |proj_HV(t) - proj_LV(t)|
%  for context X.
% =====================================================================

fprintf('Step 8: Computing value encoding distances...\n');

% Object context (cond 1=objHVrf, 2=objLVrf, 3=objHVnRF, 4=objLVnRF)
val_dist_obj = ( abs(P(1,1) - P(2,1)) + abs(P(3,1) - P(4,1)) ) / 2;

% Location context (cond 5=locHVrf, 6=locLVrf, 7=locHVnRF, 8=locLVnRF)
val_dist_loc = ( abs(P(5,1) - P(6,1)) + abs(P(7,1) - P(8,1)) ) / 2;

%% =====================================================================
%  SECTION 9 — CROSS-CONTEXT ANGLE BETWEEN VALUE AXES
%  Fit a context-specific value regression vector (difference of means,
%  denoised) for object and location contexts separately, then measure
%  the cosine angle between the two vectors.
%  θ ≈ 0°  → shared value subspace
%  θ ≈ 90° → orthogonal (distinct) subspaces
% =====================================================================

fprintf('Step 9: Cross-context angle between value axes...\n');

% Average across RF positions within each context
X_obj_HV_z = (bsxfun(@minus, psth_high_value_object_RF,    mu_n) ./ sigma_n + ...
               bsxfun(@minus, psth_high_value_object_nonRF, mu_n) ./ sigma_n) / 2;
X_obj_LV_z = (bsxfun(@minus, psth_low_value_object_RF,     mu_n) ./ sigma_n + ...
               bsxfun(@minus, psth_low_value_object_nonRF,  mu_n) ./ sigma_n) / 2;
X_loc_HV_z = (bsxfun(@minus, psth_high_value_location_RF,    mu_n) ./ sigma_n + ...
               bsxfun(@minus, psth_high_value_location_nonRF, mu_n) ./ sigma_n) / 2;
X_loc_LV_z = (bsxfun(@minus, psth_low_value_location_RF,     mu_n) ./ sigma_n + ...
               bsxfun(@minus, psth_low_value_location_nonRF,  mu_n) ./ sigma_n) / 2;

% Context-specific value vectors: difference of means, denoised, time-averaged
beta_val_obj = D * mean(X_obj_HV_z - X_obj_LV_z, 2);   % [nNeurons × 1]
beta_val_loc = D * mean(X_loc_HV_z - X_loc_LV_z, 2);

cos_theta = (beta_val_obj' * beta_val_loc) / ...
            (norm(beta_val_obj) * norm(beta_val_loc));
cos_theta = min(max(cos_theta, -1), 1);    % numerical safety
angle_deg = acosd(abs(cos_theta));         % wrap to [0,90]

fprintf('  Angle between value axes (obj vs loc): %.2f degrees\n', angle_deg);
fprintf('  cos(θ) = %.4f  |  Shared variance = %.1f%%\n', cos_theta, cos_theta^2*100);

%% =====================================================================
%  SECTION 10 — PACK RESULTS INTO STRUCT
% =====================================================================

TDR = struct();
TDR.time_axis              = time_axis;
TDR.predictor_names        = predictor_names;
TDR.axes                   = TDR_axes;           % [nNeurons × nPredictors]
TDR.beta_max               = beta_max;           % [nNeurons × nPredictors]
TDR.beta_denoised          = beta_denoised;      % [nPred × nTimeBins × nNeurons]
TDR.beta_norms             = beta_norms_all;     % [nPredictors × nTimeBins]
TDR.t_max_idx              = t_max_idx;          % peak time index per predictor
TDR.projections            = projections;        % [nConditions × nTimeBins × nPredictors]
TDR.val_dist_object        = val_dist_obj;       % [1 × nTimeBins]
TDR.val_dist_location      = val_dist_loc;       % [1 × nTimeBins]
TDR.beta_val_obj           = beta_val_obj;       % context-specific value vector
TDR.beta_val_loc           = beta_val_loc;
TDR.cos_theta              = cos_theta;
TDR.angle_value_axes_deg   = angle_deg;
TDR.nPCA                   = nPCA;
TDR.var_explained          = var_expl;
TDR.D                      = D;
TDR.design_matrix          = design_matrix;
TDR.nNeurons               = nNeurons;
TDR.nTimeBins              = nTimeBins;

assignin('base', 'TDR', TDR);

%% =====================================================================
%  SECTION 11 — VISUALIZATION
% =====================================================================

fprintf('\nStep 10: Plotting...\n');

% Color palette
c_HV   = [0.85 0.15 0.15];   % red  — high value
c_LV   = [0.15 0.35 0.85];   % blue — low value
c_obj  = [0.15 0.65 0.30];   % green — object context
c_loc  = [0.80 0.45 0.10];   % orange — location context
lw = 2.0;
t_win  = time_axis >= 0.1 & time_axis <= 0.5;   % analysis window for scalar stats

figure('Name', 'TDR: Value × Context', ...
       'Position', [40 40 1500 900], ...
       'Color', 'w');

% ── Panel 1: Regression coefficient norms over time ──────────────────
ax1 = subplot(3, 4, 1);
hold on;
clrs_pred = lines(nPredictors);
for v = 1:nPredictors
    plot(time_axis, beta_norms_all(v, :), '-', ...
        'Color', clrs_pred(v,:), 'LineWidth', lw, ...
        'DisplayName', strrep(predictor_names{v},'_',' '));
    xline(time_axis(t_max_idx(v)), 'LineStyle', '--', 'Color', [clrs_pred(v,:), 0.5], ...
          'HandleVisibility', 'off');
end
xline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time from stim onset (s)'); ylabel('‖β‖  (denoised)');
title('Regression coefficient norms');
legend('Location', 'northwest', 'FontSize', 7);
box off;

% ── Panel 2: Value axis — object context ─────────────────────────────
ax2 = subplot(3, 4, 2);
hold on;
plot(time_axis, P(1,1), '-',  'Color', c_HV,       'LineWidth', lw, 'DisplayName', 'HV  RF');
plot(time_axis, P(2,1), '-',  'Color', c_LV,        'LineWidth', lw, 'DisplayName', 'LV  RF');
plot(time_axis, P(3,1), '--', 'Color', c_HV*0.65,  'LineWidth', lw, 'DisplayName', 'HV  nonRF');
plot(time_axis, P(4,1), '--', 'Color', c_LV*0.65,  'LineWidth', lw, 'DisplayName', 'LV  nonRF');
xline(0,'LineStyle',':','Color','k','LineWidth',1,'HandleVisibility','off');
xlabel('Time (s)'); ylabel('Projection (a.u.)');
title('Value axis — Object context');
legend('Location', 'northwest', 'FontSize', 7);
box off;

% ── Panel 3: Value axis — location context ───────────────────────────
ax3 = subplot(3, 4, 3);
hold on;
plot(time_axis, P(5,1), '-',  'Color', c_HV,       'LineWidth', lw, 'DisplayName', 'HV  RF');
plot(time_axis, P(6,1), '-',  'Color', c_LV,        'LineWidth', lw, 'DisplayName', 'LV  RF');
plot(time_axis, P(7,1), '--', 'Color', c_HV*0.65,  'LineWidth', lw, 'DisplayName', 'HV  nonRF');
plot(time_axis, P(8,1), '--', 'Color', c_LV*0.65,  'LineWidth', lw, 'DisplayName', 'LV  nonRF');
xline(0,'LineStyle',':','Color','k','LineWidth',1,'HandleVisibility','off');
xlabel('Time (s)'); ylabel('Projection (a.u.)');
title('Value axis — Location context');
legend('Location', 'northwest', 'FontSize', 7);
box off;

% ── Panel 4: Value encoding distance over time ───────────────────────
ax4 = subplot(3, 4, 4);
hold on;
plot(time_axis, val_dist_obj, '-',  'Color', c_obj, 'LineWidth', lw, 'DisplayName', 'Object ctx');
plot(time_axis, val_dist_loc, '-',  'Color', c_loc, 'LineWidth', lw, 'DisplayName', 'Location ctx');
patch([time_axis(t_win), fliplr(time_axis(t_win))], ...
      [zeros(1,sum(t_win))*-1, ones(1,sum(t_win))*2], ...
      [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
xline(0,'LineStyle',':','Color','k','LineWidth',1,'HandleVisibility','off');
xlabel('Time (s)'); ylabel('|HV − LV| projection');
title('Value encoding distance');
legend('Location', 'northwest', 'FontSize', 7);
ylim([0, max([val_dist_obj, val_dist_loc])*1.15]);
box off;

% ── Panel 5: Context axis projections ────────────────────────────────
ax5 = subplot(3, 4, 5);
hold on;
plot(time_axis, P(1,2), '-',  'Color', c_obj,      'LineWidth', lw, 'DisplayName', 'Obj HV RF');
plot(time_axis, P(2,2), '--', 'Color', c_obj,      'LineWidth', lw, 'DisplayName', 'Obj LV RF');
plot(time_axis, P(5,2), '-',  'Color', c_loc,      'LineWidth', lw, 'DisplayName', 'Loc HV RF');
plot(time_axis, P(6,2), '--', 'Color', c_loc,      'LineWidth', lw, 'DisplayName', 'Loc LV RF');
xline(0,'LineStyle',':','Color','k','LineWidth',1,'HandleVisibility','off');
xlabel('Time (s)'); ylabel('Projection (a.u.)');
title('Context axis');
legend('Location', 'northwest', 'FontSize', 7);
box off;

% ── Panel 6: RF position axis ────────────────────────────────────────
ax6 = subplot(3, 4, 6);
hold on;
c_rf  = [0.55 0.10 0.75];
c_nrf = [0.10 0.75 0.75];
plot(time_axis, P(1,3), '-',  'Color', c_rf,  'LineWidth', lw, 'DisplayName', 'Obj RF');
plot(time_axis, P(3,3), '--', 'Color', c_nrf, 'LineWidth', lw, 'DisplayName', 'Obj nonRF');
plot(time_axis, P(5,3), '-',  'Color', c_rf*0.6,  'LineWidth', lw, 'DisplayName', 'Loc RF');
plot(time_axis, P(7,3), '--', 'Color', c_nrf*0.6, 'LineWidth', lw, 'DisplayName', 'Loc nonRF');
xline(0,'LineStyle',':','Color','k','LineWidth',1,'HandleVisibility','off');
xlabel('Time (s)'); ylabel('Projection (a.u.)');
title('RF position axis');
legend('Location', 'northwest', 'FontSize', 7);
box off;

% ── Panel 7: Interaction axis ────────────────────────────────────────
ax7 = subplot(3, 4, 7);
hold on;
plot(time_axis, P(1,4), '-',  'Color', c_HV,      'LineWidth', lw, 'DisplayName', 'Obj HV RF');
plot(time_axis, P(2,4), '--', 'Color', c_LV,       'LineWidth', lw, 'DisplayName', 'Obj LV RF');
plot(time_axis, P(5,4), '-',  'Color', c_HV*0.6,  'LineWidth', lw, 'DisplayName', 'Loc HV RF');
plot(time_axis, P(6,4), '--', 'Color', c_LV*0.6,  'LineWidth', lw, 'DisplayName', 'Loc LV RF');
xline(0,'LineStyle',':','Color','k','LineWidth',1,'HandleVisibility','off');
xlabel('Time (s)'); ylabel('Projection (a.u.)');
title('Value × Context interaction axis');
legend('Location', 'northwest', 'FontSize', 7);
box off;

% ── Panel 8: 2D state-space trajectory (value axis vs context axis) ──
ax8 = subplot(3, 4, 8);
t_traj = time_axis >= 0 & time_axis <= 0.6;
hold on;
cond_labels = {'Obj HV RF','Obj LV RF','Loc HV RF','Loc LV RF'};
cond_idx    = [1 2 5 6];
traj_colors = {c_HV, c_LV, c_HV*0.6, c_LV*0.6};
traj_styles = {'-', '-', '--', '--'};
for k = 1:4
    ci = cond_idx(k);
    xdata = squeeze(projections(ci, t_traj, 1));
    ydata = squeeze(projections(ci, t_traj, 2));
    plot(xdata, ydata, traj_styles{k}, 'Color', traj_colors{k}, ...
         'LineWidth', lw, 'DisplayName', cond_labels{k});
    plot(xdata(1), ydata(1), 'o', 'Color', traj_colors{k}, ...
         'MarkerFaceColor', 'w', 'MarkerSize', 7, 'HandleVisibility', 'off');  % start
    plot(xdata(end), ydata(end), 's', 'Color', traj_colors{k}, ...
         'MarkerFaceColor', traj_colors{k}, 'MarkerSize', 7, 'HandleVisibility', 'off'); % end
end
xlabel('Value axis'); ylabel('Context axis');
title('2D trajectory: value vs context');
legend('Location', 'best', 'FontSize', 7);
axis equal; box off;

% ── Panel 9: Time-averaged value encoding distance (bar) ─────────────
ax9 = subplot(3, 4, 9);
d_obj = mean(val_dist_obj(t_win));
d_loc = mean(val_dist_loc(t_win));
b = bar([1 2], [d_obj, d_loc], 0.5);
b.FaceColor = 'flat';
b.CData = [c_obj; c_loc];
b.EdgeColor = 'none';
set(gca, 'XTick', [1 2], 'XTickLabel', {'Object ctx','Location ctx'});
ylabel('Mean |HV − LV| (100–500 ms)');
title('Value separability (time-avg)');
box off;

% ── Panel 10: Cross-context angle summary ────────────────────────────
ax10 = subplot(3, 4, 10);
% Draw compass-style: two unit vectors with angle between them
ang_rad = acos(abs(cos_theta));
vec1 = [1; 0];                         % object value axis (reference)
vec2 = [cos(ang_rad); sin(ang_rad)];   % location value axis
hold on;
quiver(0, 0, vec1(1), vec1(2), 0, 'Color', c_obj, 'LineWidth', 3, ...
       'MaxHeadSize', 0.5, 'DisplayName', 'Object ctx');
quiver(0, 0, vec2(1), vec2(2), 0, 'Color', c_loc, 'LineWidth', 3, ...
       'MaxHeadSize', 0.5, 'DisplayName', 'Location ctx');
% Arc showing angle
arc_t = linspace(0, ang_rad, 50);
plot(0.3*cos(arc_t), 0.3*sin(arc_t), 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(0.35, 0.15, sprintf('θ=%.1f°', angle_deg), 'FontSize', 11);
xlim([-0.2 1.3]); ylim([-0.3 1.3]);
legend('Location', 'northwest', 'FontSize', 7);
title(sprintf('Value axis angle: %.1f°\nShared variance: %.0f%%', ...
      angle_deg, cos_theta^2*100));
axis off;

% ── Panel 11: PCA scree plot ──────────────────────────────────────────
ax11 = subplot(3, 4, 11);
n_show = min(30, nNeurons);
plot(1:n_show, var_expl(1:n_show)*100, 'k-o', ...
     'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', lw);
xline(nPCA, 'LineStyle', '--', 'Color', 'r');
yline(90, 'LineStyle', ':', 'Color', 'b');
yl = ylim;
text(nPCA - 0.3, yl(2)*0.88, sprintf('nPCA=%d', nPCA), ...
    'Color', 'r', 'FontSize', 9, 'HorizontalAlignment', 'right');
text(n_show*0.97, 91, '90%', 'Color', 'b', 'FontSize', 9, 'HorizontalAlignment', 'right');
xlabel('# PCs'); ylabel('Cumulative variance (%)');
title('PCA denoising'); box off;

% ── Panel 12: Context-specific value axis similarity matrix ──────────
ax12 = subplot(3, 4, 12);
BV = [beta_val_obj, beta_val_loc];
nrm = sqrt(sum(BV.^2, 1));
corr_mat = (BV' * BV) ./ (nrm' * nrm);
imagesc(corr_mat, [-1 1]);
colorbar;
colormap(ax12, redblue_cmap(64));
set(gca, 'XTick', [1 2], 'XTickLabel', {'Object','Location'}, ...
         'YTick', [1 2], 'YTickLabel', {'Object','Location'});
title(sprintf('Value axis similarity\ncos(θ)=%.3f', cos_theta));
axis square;

sgtitle(sprintf('TDR Analysis — Value × Context  (N=%d neurons, nPCA=%d)', ...
        nNeurons, nPCA), 'FontSize', 13, 'FontWeight', 'bold');

%% =====================================================================
%  SUMMARY PRINTOUT
% =====================================================================

fprintf('\n======== TDR Analysis Summary ========\n');
fprintf('  Neurons:      %d\n', nNeurons);
fprintf('  Time bins:    %d\n', nTimeBins);
fprintf('  PCA components retained: %d (%.1f%% var)\n', nPCA, var_expl(nPCA)*100);
fprintf('\n  Peak times of regression norms:\n');
for v = 1:nPredictors
    fprintf('    %-25s  t = %.3f s\n', predictor_names{v}, time_axis(t_max_idx(v)));
end
fprintf('\n  Value encoding distance (100–500 ms avg):\n');
fprintf('    Object context:   %.4f\n', mean(val_dist_obj(t_win)));
fprintf('    Location context: %.4f\n', mean(val_dist_loc(t_win)));
fprintf('\n  Cross-context value axis angle:\n');
fprintf('    θ = %.2f deg   cos(θ) = %.4f\n', angle_deg, cos_theta);
fprintf('    → %.1f%% shared variance between value axes\n', cos_theta^2*100);
fprintf('======================================\n');

end  % function TDR_ValueContext

%% =====================================================================
%  LOCAL HELPER: Red-white-blue colormap
% =====================================================================
function cmap = redblue_cmap(n)
% n-level diverging red-white-blue colormap
half = ceil(n/2);
r = [linspace(0.90, 1.0, half), linspace(1.0, 0.15, half)];
g = [linspace(0.20, 1.0, half), linspace(1.0, 0.30, half)];
b = [linspace(0.20, 1.0, half), linspace(1.0, 0.85, half)];
cmap = [r(1:n)', g(1:n)', b(1:n)'];
end
