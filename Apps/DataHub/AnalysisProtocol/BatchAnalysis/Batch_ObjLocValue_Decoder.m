% Batch_ObjLocValue_Decoder.m
% Standalone value decoder script. Loads pre-computed population PSTH matrices
% saved by Batch_ObjLocValue_Organize (OutputFlag=true) and runs the full value
% axis geometry and decoding analysis.
%
% Prerequisites: run Batch_ObjLocValue_Organize with OutputFlag=true first.
%
% Loaded variables:
%   psth_high_value_object_RF   [N × T] — stable context, high-value trials, z-scored
%   psth_low_value_object_RF    [N × T] — stable context, low-value trials, z-scored
%   psth_high_value_location_RF [N × T] — oneDR context, high-value trials, z-scored
%   psth_low_value_location_RF  [N × T] — oneDR context, low-value trials, z-scored
%   time_axis_TDR               [1 × T] — time axis in ms, aligned to target onset

%% Load saved inputs
AreaLabel = 'FEF';   % 'FEF' or 'SC'
SaveRoot  = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/ObjectLocationPSTH';
LoadFile  = fullfile(SaveRoot, AreaLabel, sprintf('ValueDecoder_Input_%s.mat', AreaLabel));
load(LoadFile);
fprintf('Loaded ValueDecoder inputs from %s\n', LoadFile);

%% ── Value decoder — axis construction ───────────────────────────────────────
% ValueDecoder_[context] [N × T]: value selectivity per neuron per time bin
%   = mean(high value trials) − mean(low value trials)
% ValueAxis [N × 1]: scalar weight per neuron, averaged over ±5 bins around
%   each neuron's peak selectivity time → direction of value in state space

T_tdr     = size(psth_high_value_object_RF, 2);
N_neurons = size(psth_high_value_object_RF, 1);

%psth_high_value_object_RF= all_psths_z()
%{
all_psths = cat(3, ...
    psth_high_value_object_RF,   psth_low_value_object_RF, ...
    psth_high_value_location_RF, psth_low_value_location_RF);


nNeurons = size(psth_high_value_object_RF,1);
nTimeBins= size(psth_high_value_object_RF,2);
nConditions= 4;

all_flat  = reshape(all_psths, nNeurons, nTimeBins * nConditions);
mu_n      = mean(all_flat, 2);
sigma_n   = std(all_flat, 0, 2) + 1e-8;

psth_high_value_object_RF = (psth_high_value_object_RF - mu_n) ./ sigma_n;
psth_low_value_object_RF = (psth_low_value_object_RF - mu_n) ./ sigma_n;

psth_high_value_location_RF = (psth_high_value_location_RF - mu_n) ./ sigma_n;
psth_low_value_location_RF = (psth_low_value_location_RF - mu_n) ./ sigma_n;
%}

ValueDecoder_stable = psth_high_value_object_RF   - psth_low_value_object_RF;    % [N × T]
ValueDecoder_oneDR  = psth_high_value_location_RF - psth_low_value_location_RF;  % [N × T]

% Population-level peak: find the time of maximum average value difference
% across neurons. All neurons share this single window.
[~, peak_index]       = max(abs(mean(ValueDecoder_stable, 1)));   % scalar
[~, peak_index_onedr] = max(abs(mean(ValueDecoder_oneDR,  1)));   % scalar

win_s = max(peak_index-5,       1) : min(peak_index+5,       T_tdr);
win_d = max(peak_index_onedr-5, 1) : min(peak_index_onedr+5, T_tdr);

ValueAxis_stable = mean(ValueDecoder_stable(:, win_s), 2);   % [N × 1]
ValueAxis_oneDR  = mean(ValueDecoder_oneDR(:,  win_d), 2);   % [N × 1]

% Unit vectors along each value axis in N-dimensional neural state space
w_stable = ValueAxis_stable / norm(ValueAxis_stable);   % [N × 1]
w_oneDR  = ValueAxis_oneDR  / norm(ValueAxis_oneDR);    % [N × 1]

% Gram-Schmidt: component of w_oneDR orthogonal to w_stable
% Needed as the y-axis of the 2D visualization basis — defined here so it is
% available for both the bootstrap CI and the VAF/geometry sections below.
w_oneDR_orth_raw = w_oneDR - dot(w_oneDR, w_stable) * w_stable;
w_oneDR_orth     = w_oneDR_orth_raw / norm(w_oneDR_orth_raw);   % [N × 1]

% Angle between axes — decoder is always high-low, so sign is meaningful
cos_theta = dot(w_stable, w_oneDR);
theta_deg = acosd(cos_theta);
fprintf('\nAngle between stable and oneDR value axes: %.2f deg\n', theta_deg);

% perm_p placeholder — replaced by bootstrap CI (see below).
% A random-axis permutation test is not used: in N-dimensional space the null
% concentrates so tightly at 90° that even a 1° deviation reaches p < 0.05,
% making the p-value sensitive to N rather than to scientific effect size.
% The bootstrap CI directly captures sampling uncertainty in the angle.
perm_p = NaN;

%% ── Bootstrap CI for each value axis direction ───────────────────────────────
% Resample neurons with replacement. Each resample gives a new w_stable_boot
% and w_oneDR_boot. Project onto the 2D basis {w_stable, w_oneDR_orth} to get
% the angle of each boot axis in the vector diagram plane.
% Sign-align each boot axis with the original so CI doesn't span a full 360°.
N_boot = 500;
boot_angle_stable = zeros(1, N_boot);
boot_angle_oneDR  = zeros(1, N_boot);
boot_theta        = zeros(1, N_boot);

for b = 1:N_boot
    idx = randi(N_neurons, N_neurons, 1);

    % Stable axis — population peak time on resampled neurons
    VD_s = psth_high_value_object_RF(idx,:) - psth_low_value_object_RF(idx,:);
    [~, pk_s] = max(abs(mean(VD_s, 1)));
    win_bs = max(pk_s-5,1) : min(pk_s+5, T_tdr);
    w_s = mean(VD_s(:, win_bs), 2);
    w_s = w_s / norm(w_s);
    % Sign-align w_s to w_stable: resampling can flip the sign when more
    % low-value neurons are drawn, which would give a spurious ~180°-θ angle.
    if dot(w_s, w_stable) < 0, w_s = -w_s; end

    % oneDR axis — population peak time on resampled neurons
    VD_d = psth_high_value_location_RF(idx,:) - psth_low_value_location_RF(idx,:);
    [~, pk_d] = max(abs(mean(VD_d, 1)));
    win_bd = max(pk_d-5,1) : min(pk_d+5, T_tdr);
    w_d = mean(VD_d(:, win_bd), 2);
    w_d = w_d / norm(w_d);
    % Sign-align w_d to w_oneDR for the same reason.
    if dot(w_d, w_oneDR) < 0, w_d = -w_d; end

    % 2D angles in the {w_stable, w_oneDR_orth} basis
    boot_angle_stable(b) = atan2d(dot(w_s, w_oneDR_orth), dot(w_s, w_stable));
    boot_angle_oneDR(b)  = atan2d(dot(w_d, w_oneDR_orth), dot(w_d, w_stable));
    boot_theta(b)        = acosd(dot(w_s, w_d));
end

ci_s = prctile(boot_angle_stable, [5 95]);   % 90% CI of stable axis angle in 2D
ci_d = prctile(boot_angle_oneDR,  [5 95]);   % 90% CI of oneDR  axis angle in 2D
ci_t = prctile(boot_theta,        [5 95]);   % 90% CI of the angle between axes

fprintf('\n--- Bootstrap CI (n=%d resamples) ---\n', N_boot);
fprintf('  Angle between value axes:  %.1f° [%.1f°, %.1f°]  (90%% CI)\n', theta_deg, ci_t(1), ci_t(2));

% Bootstrap CI figure for the angle between the two axes
figure('Name', 'Bootstrap CI: angle between value axes', 'Position', [100 100 500 380]);
hold on;
histogram(boot_theta, 30, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');
xline(theta_deg, 'color', 'r','LineStyle','-',  'LineWidth', 2,   'DisplayName', sprintf('Observed: %.1f°', theta_deg));
xline(ci_t(1),   'color','b','LineStyle','--', 'LineWidth', 1.5, 'DisplayName', sprintf('90%% CI:  [%.1f°, %.1f°]', ci_t(1), ci_t(2)));
xline(ci_t(2),   'color','b','LineStyle','--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('Angle between value axes (deg)');
ylabel('Bootstrap count');
title(sprintf('Bootstrap CI for axis angle  (n=%d)', N_boot));
legend('Location', 'best'); box off;

%% ── Project all conditions onto each value axis ──────────────────────────────
% Each projection is a [1 × T] timecourse: how much the population moves
% along the value axis when conditions are varied across the two contexts.
proj_stable_high_obj = w_stable' * psth_high_value_object_RF;     % stable axis, stable ctx high
proj_stable_low_obj  = w_stable' * psth_low_value_object_RF;      % stable axis, stable ctx low
proj_stable_high_loc = w_stable' * psth_high_value_location_RF;   % stable axis, oneDR ctx high
proj_stable_low_loc  = w_stable' * psth_low_value_location_RF;    % stable axis, oneDR ctx low

proj_oneDR_high_obj  = w_oneDR'  * psth_high_value_object_RF;     % oneDR axis, stable ctx high
proj_oneDR_low_obj   = w_oneDR'  * psth_low_value_object_RF;
proj_oneDR_high_loc  = w_oneDR'  * psth_high_value_location_RF;   % oneDR axis, oneDR ctx high
proj_oneDR_low_loc   = w_oneDR'  * psth_low_value_location_RF;

figure('Name', 'Projections onto value axes', 'Position', [100 100 1000 420]);

subplot(1,2,1); hold on;
plot(time_axis_TDR, proj_stable_high_obj, 'b-',  'LineWidth', 1.5, 'DisplayName', 'Stable high');
plot(time_axis_TDR, proj_stable_low_obj,  'b--', 'LineWidth', 1.5, 'DisplayName', 'Stable low');
plot(time_axis_TDR, proj_stable_high_loc, 'r-',  'LineWidth', 1.5, 'DisplayName', 'OneDR high');
plot(time_axis_TDR, proj_stable_low_loc,  'r--', 'LineWidth', 1.5, 'DisplayName', 'OneDR low');
xline(0, 'LineStyle', '--', 'Color', 'k');
xlabel('Time (ms)'); ylabel('Projection (a.u.)');
title('Stable value axis'); legend('Location', 'best'); box off;

subplot(1,2,2); hold on;
plot(time_axis_TDR, proj_oneDR_high_obj,  'b-',  'LineWidth', 1.5, 'DisplayName', 'Stable high');
plot(time_axis_TDR, proj_oneDR_low_obj,   'b--', 'LineWidth', 1.5, 'DisplayName', 'Stable low');
plot(time_axis_TDR, proj_oneDR_high_loc,  'r-',  'LineWidth', 1.5, 'DisplayName', 'OneDR high');
plot(time_axis_TDR, proj_oneDR_low_loc,   'r--', 'LineWidth', 1.5, 'DisplayName', 'OneDR low');
xline(0, 'LineStyle', '--', 'Color', 'k');
xlabel('Time (ms)'); ylabel('Projection (a.u.)');
title('OneDR value axis'); legend('Location', 'best'); box off;

%% ── Value-difference VAF: how much of the high−low signal does each axis capture ──
% The value signal is the difference in population response between high and low value:
%   D_stable = ValueDecoder_stable  [N × T]  (already computed above)
%   D_oneDR  = ValueDecoder_oneDR   [N × T]
%
% For axis w, the fraction of the value-difference signal captured in context c:
%   VAF_c(w) = ||w' * D_c||²_F / ||D_c||²_F
%
% No centering needed — D is already a difference, so it is naturally mean-free
% across the two conditions. The denominator is the total power in the value signal,
% not the full data variance, so temporal PSTH dynamics do not inflate the denominator.

pow_stable = sum(ValueDecoder_stable(:).^2);   % total value-difference power, stable
pow_oneDR  = sum(ValueDecoder_oneDR(:).^2);    % total value-difference power, oneDR
pow_total  = pow_stable + pow_oneDR;

% --- Global summary ---
fprintf('\n--- Value-difference VAF: %% of high-low signal captured by each axis ---\n');
fprintf('  %-22s  %12s  %12s  %14s\n', 'Axis', 'Stable ctx', 'OneDR ctx', 'Both combined');

axes_labels = {'Stable axis',  'OneDR axis',  '2D subspace'};
axes_w1     = { w_stable,       w_oneDR,        w_stable    };
axes_w2     = { [],             [],             w_oneDR_orth };  % second basis for 2D

for k = 1:3
    w1 = axes_w1{k};
    num_s = sum((w1' * ValueDecoder_stable).^2);
    num_d = sum((w1' * ValueDecoder_oneDR ).^2);
    if ~isempty(axes_w2{k})   % 2D subspace: add orthogonal component
        w2 = axes_w2{k};
        num_s = num_s + sum((w2' * ValueDecoder_stable).^2);
        num_d = num_d + sum((w2' * ValueDecoder_oneDR ).^2);
    end
    fprintf('  %-22s  %10.1f%%  %10.1f%%  %12.1f%%\n', axes_labels{k}, ...
        num_s/pow_stable*100, num_d/pow_oneDR*100, (num_s+num_d)/pow_total*100);
end

%% ── Shared and context-differential value axes ───────────────────────────────
% The 2D subspace of {w_stable, w_oneDR} admits a symmetric orthonormal basis:
%   w_shared — bisector (sum of unit vectors): common value direction across contexts
%   w_diff   — anti-bisector (difference): direction along which contexts disagree
%
% Orthogonality is guaranteed:
%   dot(w_stable + w_oneDR, w_stable - w_oneDR) = ||w_stable||² - ||w_oneDR||² = 0
% Together they span the same 2D subspace, so VAF_shared_ax + VAF_diff_ax = VAF_2D.

% Note: if dot(w_stable,w_oneDR) < 0, neurons predominantly have opposite
% preferences across contexts. w_shared_raw will be near zero in that case —
% the shared axis is ill-defined, and the geometry is dominated by w_diff.
w_shared_raw = w_stable + w_oneDR;
w_diff_raw   = w_stable - w_oneDR;
if norm(w_shared_raw) < 1e-6
    warning('w_shared is near-zero: axes are nearly anti-parallel (opposite value preferences dominate). Shared axis is not meaningful.');
end
w_shared     = w_shared_raw / norm(w_shared_raw);   % [N × 1]
w_diff       = w_diff_raw   / norm(w_diff_raw);     % [N × 1]

fprintf('\n--- Shared / Context-Differential Axes ---\n');
fprintf('  Orthogonality check (should be ~0): %.2e\n', dot(w_shared, w_diff));
fprintf('  Angle: w_shared vs w_stable = %.1f deg\n', acosd(dot(w_shared, w_stable)));
fprintf('  Angle: w_shared vs w_oneDR  = %.1f deg\n', acosd(dot(w_shared, w_oneDR)));

%% Project all conditions onto shared and context-differential axes
proj_shared_high_obj = w_shared' * psth_high_value_object_RF;     % [1 × T]
proj_shared_low_obj  = w_shared' * psth_low_value_object_RF;
proj_shared_high_loc = w_shared' * psth_high_value_location_RF;
proj_shared_low_loc  = w_shared' * psth_low_value_location_RF;

proj_diff_high_obj   = w_diff' * psth_high_value_object_RF;
proj_diff_low_obj    = w_diff' * psth_low_value_object_RF;
proj_diff_high_loc   = w_diff' * psth_high_value_location_RF;
proj_diff_low_loc    = w_diff' * psth_low_value_location_RF;

% On w_shared: both contexts should show high > low separation in the same direction.
% On w_diff:  contexts differ — stable high and oneDR low (or vice versa) may
%             both project high, revealing context-specific value coding.
figure('Name', 'Shared vs differential value axes', 'Position', [100 100 1000 420]);

subplot(1,2,1); hold on;
plot(time_axis_TDR, proj_shared_high_obj, 'b-',  'LineWidth', 1.5, 'DisplayName', 'Stable high');
plot(time_axis_TDR, proj_shared_low_obj,  'b--', 'LineWidth', 1.5, 'DisplayName', 'Stable low');
plot(time_axis_TDR, proj_shared_high_loc, 'r-',  'LineWidth', 1.5, 'DisplayName', 'OneDR high');
plot(time_axis_TDR, proj_shared_low_loc,  'r--', 'LineWidth', 1.5, 'DisplayName', 'OneDR low');
xline(0, 'LineStyle', '--', 'Color', 'k');
xlabel('Time (ms)'); ylabel('Projection (a.u.)');
title('Shared value axis'); legend('Location', 'best'); box off;

subplot(1,2,2); hold on;
plot(time_axis_TDR, proj_diff_high_obj,  'b-',  'LineWidth', 1.5, 'DisplayName', 'Stable high');
plot(time_axis_TDR, proj_diff_low_obj,   'b--', 'LineWidth', 1.5, 'DisplayName', 'Stable low');
plot(time_axis_TDR, proj_diff_high_loc,  'r-',  'LineWidth', 1.5, 'DisplayName', 'OneDR high');
plot(time_axis_TDR, proj_diff_low_loc,   'r--', 'LineWidth', 1.5, 'DisplayName', 'OneDR low');
xline(0, 'LineStyle', '--', 'Color', 'k');
xlabel('Time (ms)'); ylabel('Projection (a.u.)');
title('Context-differential value axis'); legend('Location', 'best'); box off;

%% 2D trajectory: population path in shared × context-differential plane
% Each curve is (proj_shared(t), proj_diff(t)) as time advances; arrows show
% direction of increasing time near target onset.
figure('Name', 'Trajectory: shared × context-differential (2D)', ...
    'Position', [130 130 620 540]);
ax_traj = axes; hold(ax_traj, 'on');

traj_specs = {
    proj_shared_high_obj, proj_diff_high_obj, [0.15 0.45 0.85], '-', 2.0, 'Stable, high value';
    proj_shared_low_obj,  proj_diff_low_obj,  [0.45 0.65 0.95], '--', 1.8, 'Stable, low value';
    proj_shared_high_loc, proj_diff_high_loc, [0.85 0.25 0.12], '-', 2.0, 'OneDR, high value';
    proj_shared_low_loc,  proj_diff_low_loc,  [0.95 0.55 0.35], '--', 1.8, 'OneDR, low value'
};

for tr = 1:size(traj_specs, 1)
    xs = traj_specs{tr, 1}(:)';
    ys = traj_specs{tr, 2}(:)';
    col = traj_specs{tr, 3};
    lst = traj_specs{tr, 4};
    lw  = traj_specs{tr, 5};
    leg = traj_specs{tr, 6};
    plot(ax_traj, xs, ys, lst, 'Color', col, 'LineWidth', lw, 'DisplayName', leg);
end

% Target onset (t = 0 ms): index on shared time axis
[~, i_t0] = min(abs(time_axis_TDR(:)));
for tr = 1:size(traj_specs, 1)
    xs = traj_specs{tr, 1}(:)';
    ys = traj_specs{tr, 2}(:)';
    col = traj_specs{tr, 3};
    scatter(ax_traj, xs(i_t0), ys(i_t0), 55, col, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.6, 'HandleVisibility', 'off');
end

% Small arrows: direction of time (t0 -> t0+3 bins) on each trajectory
dt_arrow = min(3, T_tdr - i_t0);
if dt_arrow >= 1 && i_t0 + dt_arrow <= T_tdr
    for tr = 1:size(traj_specs, 1)
        xs = traj_specs{tr, 1}(:)';
        ys = traj_specs{tr, 2}(:)';
        col = traj_specs{tr, 3};
        quiver(ax_traj, xs(i_t0), ys(i_t0), xs(i_t0+dt_arrow)-xs(i_t0), ys(i_t0+dt_arrow)-ys(i_t0), ...
            0, 'Color', col, 'LineWidth', 1.2, 'MaxHeadSize', 0.35, 'HandleVisibility', 'off');
    end
end

xlabel(ax_traj, 'Shared value axis (projection, a.u.)');
ylabel(ax_traj, 'Context-differential value axis (projection, a.u.)');
title(ax_traj, 'Population trajectories: shared vs context-specific value');
legend(ax_traj, 'Location', 'best');
grid(ax_traj, 'on');
axis(ax_traj, 'equal');
box(ax_traj, 'off');

%% ── Vector diagram: geometry of value axes in 2D ────────────────────────────
% All four axes (w_stable, w_oneDR, w_shared, w_diff) live in the same 2D plane.
% Project onto the orthonormal basis {w_stable, w_oneDR_orth} so the coordinates
% are exact — not an approximation: w_stable=(1,0), w_oneDR=(cosθ, sinθ).

% 2D projections (dot products with the two basis vectors)
v_stable = [dot(w_stable, w_stable),  dot(w_stable, w_oneDR_orth)];  % = [1, 0]
v_oneDR  = [dot(w_oneDR,  w_stable),  dot(w_oneDR,  w_oneDR_orth)];  % = [cosθ, sinθ]
v_shared = [dot(w_shared, w_stable),  dot(w_shared, w_oneDR_orth)];  % bisector
v_diff   = [dot(w_diff,   w_stable),  dot(w_diff,   w_oneDR_orth)];  % anti-bisector

% Standardise orientation: keep w_diff in upper half-plane so both derived
% axes appear between w_stable and w_oneDR (just a display convention).
if v_diff(2) < 0,  v_diff   = -v_diff;   end
if v_shared(2) < 0, v_shared = -v_shared; end

col_stable = [0.2 0.5 0.9];
col_oneDR  = [0.9 0.4 0.1];
col_shared = [0.1 0.7 0.2];
col_diff   = [0.7 0.2 0.7];

vecs   = {v_stable, v_oneDR, v_shared, v_diff};
cols   = {col_stable, col_oneDR, col_shared, col_diff};
lws    = {2.5, 2.5, 1.8, 1.8};
lstys  = {'-', '-', '--', '--'};
vlbls  = {'w_{stable}', 'w_{oneDR}', 'w_{shared}', 'w_{diff}'};

scale = 0.88;   % vector length (leaves room for labels)
ahl   = 0.055;  % arrowhead length
ahw   = 0.030;  % arrowhead half-width

figure('Name', 'Value axis geometry', 'Position', [100 100 500 520]);
ax_geo = axes; hold on; axis equal;


% ── Bootstrap CI wedges for each axis ─────────────────────────────────────
% Shaded cone around each arrow showing the 90% bootstrap CI of that axis
% direction (uncertainty from finite neuron sampling). Narrow = stable estimate;
% wide = sensitive to which neurons are included.
arc_bs = linspace(deg2rad(ci_s(1)), deg2rad(ci_s(2)), 60);
fill(ax_geo, [0, scale*cos(arc_bs), 0], [0, scale*sin(arc_bs), 0], col_stable, ...
    'FaceAlpha', 0.20, 'EdgeColor', 'none', ...
    'DisplayName', sprintf('w_{stable} 90%% CI (%.1f°–%.1f°)', ci_s(1), ci_s(2)));

arc_bd = linspace(deg2rad(ci_d(1)), deg2rad(ci_d(2)), 60);
fill(ax_geo, [0, scale*cos(arc_bd), 0], [0, scale*sin(arc_bd), 0], col_oneDR, ...
    'FaceAlpha', 0.20, 'EdgeColor', 'none', ...
    'DisplayName', sprintf('w_{oneDR} 90%% CI (%.1f°–%.1f°)', ci_d(1), ci_d(2)));

for k = 1:4
    vx  = vecs{k}(1) * scale;
    vy  = vecs{k}(2) * scale;
    ang = atan2(vy, vx);

    % Shaft — stops just before tip so arrowhead sits cleanly
    plot([0, vx - ahl*cos(ang)], [0, vy - ahl*sin(ang)], ...
        lstys{k}, 'Color', cols{k}, 'LineWidth', lws{k}, 'HandleVisibility','off');

    % Filled arrowhead triangle
    perp   = [-sin(ang), cos(ang)];
    tip    = [vx, vy];
    base_l = tip - ahl*[cos(ang), sin(ang)] + ahw*perp;
    base_r = tip - ahl*[cos(ang), sin(ang)] - ahw*perp;
    fill([tip(1), base_l(1), base_r(1)], [tip(2), base_l(2), base_r(2)], ...
        cols{k}, 'EdgeColor', cols{k}, 'HandleVisibility','off');

    % Invisible proxy for the legend
    plot(NaN, NaN, lstys{k}, 'Color', cols{k}, 'LineWidth', lws{k}, 'DisplayName', vlbls{k});

    % Label at vector tip (offset outward along the vector direction)
    lbl_off = 0.10;
    text(vx + lbl_off*cos(ang), vy + lbl_off*sin(ang), vlbls{k}, ...
        'Color', cols{k}, 'FontSize', 11, 'FontWeight','bold', 'HorizontalAlignment','center');
end

% Arc marking the full angle θ between w_stable and w_oneDR
arc_r    = 0.28;
theta_rad = atan2(v_oneDR(2), v_oneDR(1));
arc_pts  = linspace(0, theta_rad, 80);
plot(arc_r*cos(arc_pts), arc_r*sin(arc_pts), 'k-', 'LineWidth', 1.2, 'HandleVisibility','off');
theta_mid_rad = theta_rad / 2;
text(0.37*cos(theta_mid_rad), 0.37*sin(theta_mid_rad), sprintf('\\theta=%.1f°', theta_deg), ...
    'FontSize', 10, 'HorizontalAlignment','center', 'FontWeight','bold');

% Small arc marking θ/2 (the bisector angle)
arc_r2   = 0.16;
ang_sh   = atan2(v_shared(2), v_shared(1));
arc_pts2 = linspace(0, ang_sh, 40);
plot(arc_r2*cos(arc_pts2), arc_r2*sin(arc_pts2), ':', 'Color',[0.55 0.55 0.55], ...
    'LineWidth',0.9,'HandleVisibility','off');
text(0.22*cos(ang_sh/2), 0.22*sin(ang_sh/2), sprintf('%.1f°', theta_deg/2), ...
    'FontSize', 9, 'HorizontalAlignment','center', 'Color',[0.45 0.45 0.45]);

% Right-angle marker showing w_shared ⊥ w_diff
sq = 0.07;
u_sh = v_shared / norm(v_shared);
u_df = v_diff   / norm(v_diff);
sq_pts = [sq*u_sh; sq*u_sh + sq*u_df; sq*u_df];
plot(sq_pts(:,1), sq_pts(:,2), 'k-', 'LineWidth', 0.8, 'HandleVisibility','off');

% Origin dot
plot(0, 0, 'k.', 'MarkerSize', 12, 'HandleVisibility','off');

xlim([-0.25 1.55]);  ylim([-0.45 1.35]);
set(ax_geo, 'XTick',[], 'YTick',[]);
xlabel('w_{stable} direction');
ylabel('w_{oneDR\_orth} direction');
title(sprintf('Value axis geometry  \\theta=%.1f° [%.1f°, %.1f°]', theta_deg, ci_t(1), ci_t(2)));
legend('Location','northwest'); box off;

%% Value-difference VAF: shared vs context-differential axes
% Same logic — denominator is the value-difference power per context.
% w_shared and w_diff are orthogonal, so their VAFs add up to the 2D subspace VAF.

num_sh_s = sum((w_shared' * ValueDecoder_stable).^2);
num_sh_d = sum((w_shared' * ValueDecoder_oneDR ).^2);
num_df_s = sum((w_diff'   * ValueDecoder_stable).^2);
num_df_d = sum((w_diff'   * ValueDecoder_oneDR ).^2);

VAF_sh_stable = num_sh_s / pow_stable;
VAF_sh_oneDR  = num_sh_d / pow_oneDR;
VAF_df_stable = num_df_s / pow_stable;
VAF_df_oneDR  = num_df_d / pow_oneDR;
VAF_sh_total  = (num_sh_s + num_sh_d) / pow_total;
VAF_df_total  = (num_df_s + num_df_d) / pow_total;

fprintf('  %-22s  %10.1f%%  %10.1f%%  %12.1f%%\n', 'Shared axis',       VAF_sh_stable*100, VAF_sh_oneDR*100, VAF_sh_total*100);
fprintf('  %-22s  %10.1f%%  %10.1f%%  %12.1f%%\n', 'Context-diff axis', VAF_df_stable*100, VAF_df_oneDR*100, VAF_df_total*100);
fprintf('  %-22s  %10.1f%%  %10.1f%%  %12.1f%%\n', 'Shared + diff (=2D)', ...
    (VAF_sh_stable+VAF_df_stable)*100, (VAF_sh_oneDR+VAF_df_oneDR)*100, (VAF_sh_total+VAF_df_total)*100);
fprintf('  Fraction of 2D that is shared:  %.1f%%\n', VAF_sh_total/(VAF_sh_total+VAF_df_total)*100);
fprintf('  Fraction of 2D that is context-diff: %.1f%%\n', VAF_df_total/(VAF_sh_total+VAF_df_total)*100);

%% Value-difference VAF timecourse — shared vs context-differential, per context
% At each time t: what fraction of the value-difference vector (high-low at that bin)
% is captured by the shared axis vs the context-differential axis?
VAF_sh_stable_t = zeros(1, T_tdr);
VAF_sh_oneDR_t  = zeros(1, T_tdr);
VAF_df_stable_t = zeros(1, T_tdr);
VAF_df_oneDR_t  = zeros(1, T_tdr);

for t = 1:T_tdr
    d_s = ValueDecoder_stable(:, t);   % [N × 1] value difference at time t
    d_d = ValueDecoder_oneDR(:,  t);
    ps_t = d_s' * d_s;
    pd_t = d_d' * d_d;
    if ps_t > 1e-12
        VAF_sh_stable_t(t) = (w_shared' * d_s)^2 / ps_t;
        VAF_df_stable_t(t) = (w_diff'   * d_s)^2 / ps_t;
    end
    if pd_t > 1e-12
        VAF_sh_oneDR_t(t)  = (w_shared' * d_d)^2 / pd_t;
        VAF_df_oneDR_t(t)  = (w_diff'   * d_d)^2 / pd_t;
    end
end

figure('Name', 'Value VAF: shared vs context-diff', 'Position', [100 100 1000 420]);

subplot(1,2,1); hold on;
fill([time_axis_TDR, fliplr(time_axis_TDR)], ...
     [(VAF_sh_stable_t+VAF_df_stable_t)*100, zeros(1,T_tdr)], ...
     [0.85 0.85 0.85], 'FaceAlpha',0.4,'EdgeColor','none','HandleVisibility','off');
plot(time_axis_TDR, (VAF_sh_stable_t+VAF_df_stable_t)*100, 'k-',  'LineWidth',2,   'DisplayName','2D subspace');
plot(time_axis_TDR, VAF_sh_stable_t*100, 'Color',[0.1 0.7 0.2], 'LineWidth',1.5, 'DisplayName','Shared axis');
plot(time_axis_TDR, VAF_df_stable_t*100, 'Color',[0.7 0.2 0.7], 'LineWidth',1.5, 'DisplayName','Context-diff axis');
xline(0, 'LineStyle','--','Color','k');
xlabel('Time (ms)'); ylabel('% value difference captured');
title(sprintf('Stable context  (global: shared %.0f%%, diff %.0f%%)', VAF_sh_stable*100, VAF_df_stable*100));
legend('Location','best'); box off; ylim([0 100]);

subplot(1,2,2); hold on;
fill([time_axis_TDR, fliplr(time_axis_TDR)], ...
     [(VAF_sh_oneDR_t+VAF_df_oneDR_t)*100, zeros(1,T_tdr)], ...
     [0.85 0.85 0.85], 'FaceAlpha',0.4,'EdgeColor','none','HandleVisibility','off');
plot(time_axis_TDR, (VAF_sh_oneDR_t+VAF_df_oneDR_t)*100,  'k-',  'LineWidth',2,   'DisplayName','2D subspace');
plot(time_axis_TDR, VAF_sh_oneDR_t*100, 'Color',[0.1 0.7 0.2], 'LineWidth',1.5, 'DisplayName','Shared axis');
plot(time_axis_TDR, VAF_df_oneDR_t*100, 'Color',[0.7 0.2 0.7], 'LineWidth',1.5, 'DisplayName','Context-diff axis');
xline(0, 'LineStyle','--','Color','k');
xlabel('Time (ms)'); ylabel('% value difference captured');
title(sprintf('OneDR context  (global: shared %.0f%%, diff %.0f%%)', VAF_sh_oneDR*100, VAF_df_oneDR*100));
legend('Location','best'); box off; ylim([0 100]);

%% ── Pseudo-population decoding: within-context and cross-context ─────────────
% At each time bin, add independent Gaussian noise (SD = 1, appropriate for
% z-scored PSTHs) to the mean population response to generate pseudo-trials.
% Project onto the value axis and classify by comparing to the threshold
% (midpoint of the two class means computed from the TRAINING context).
%
%   Within-context:  axis + threshold trained and tested on the same context
%   Cross-context:   axis + threshold from context A, tested on context B
%
% Chance level = 50%. Accuracy is the fraction of pseudo-trials correctly
% classified, averaged over N_pseudo trials per class.

N_pseudo = 1000;   % pseudo-trials per class per time bin

acc_within_stable_t = zeros(1, T_tdr);
acc_within_oneDR_t  = zeros(1, T_tdr);
acc_cross_s_on_d_t  = zeros(1, T_tdr);   % stable axis + stable threshold → oneDR data
acc_cross_d_on_s_t  = zeros(1, T_tdr);   % oneDR  axis + oneDR  threshold → stable data

for t = 1:T_tdr
    r_hi_s = psth_high_value_object_RF(:, t);     % [N × 1] mean response at time t
    r_lo_s = psth_low_value_object_RF(:, t);
    r_hi_d = psth_high_value_location_RF(:, t);
    r_lo_d = psth_low_value_location_RF(:, t);

    % Pseudo-trials: mean + Gaussian noise, [N × N_pseudo]
    pt_hi_s = r_hi_s + randn(N_neurons, N_pseudo);
    pt_lo_s = r_lo_s + randn(N_neurons, N_pseudo);
    pt_hi_d = r_hi_d + randn(N_neurons, N_pseudo);
    pt_lo_d = r_lo_d + randn(N_neurons, N_pseudo);

    % Thresholds from training-context true means (no leakage from pseudo-trials)
    thr_s = (w_stable' * r_hi_s + w_stable' * r_lo_s) / 2;
    thr_d = (w_oneDR'  * r_hi_d + w_oneDR'  * r_lo_d) / 2;

    % Within stable
    acc_within_stable_t(t) = mean([w_stable'*pt_hi_s > thr_s, w_stable'*pt_lo_s < thr_s]);

    % Within oneDR
    acc_within_oneDR_t(t)  = mean([w_oneDR'*pt_hi_d > thr_d, w_oneDR'*pt_lo_d < thr_d]);

    % Cross: stable axis + stable threshold → oneDR context
    acc_cross_s_on_d_t(t)  = mean([w_stable'*pt_hi_d > thr_s, w_stable'*pt_lo_d < thr_s]);

    % Cross: oneDR axis + oneDR threshold → stable context
    acc_cross_d_on_s_t(t)  = mean([w_oneDR'*pt_hi_s > thr_d, w_oneDR'*pt_lo_s < thr_d]);
end

figure('Name', 'Pseudo-population decoding accuracy', 'Position', [100 100 950 420]);

subplot(1,2,1); hold on;
plot(time_axis_TDR, acc_within_stable_t*100, 'b-', 'LineWidth', 2, 'DisplayName', 'Stable context');
plot(time_axis_TDR, acc_within_oneDR_t *100, 'r-', 'LineWidth', 2, 'DisplayName', 'OneDR context');
yline(50, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(0,  'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time (ms)'); ylabel('Accuracy (%)');
title('Within-context decoding');
legend('Location', 'best'); box off; ylim([40 100]);

subplot(1,2,2); hold on;
plot(time_axis_TDR, acc_cross_s_on_d_t*100, 'b-', 'LineWidth', 2, 'DisplayName', 'Stable axis → OneDR');
plot(time_axis_TDR, acc_cross_d_on_s_t*100, 'r-', 'LineWidth', 2, 'DisplayName', 'OneDR axis → Stable');
yline(50, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(0,  'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time (ms)'); ylabel('Accuracy (%)');
title('Cross-context decoding');
legend('Location', 'best'); box off; ylim([40 100]);
