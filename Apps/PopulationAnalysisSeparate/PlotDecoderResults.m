function PlotDecoderResults(Results, varargin)
% PlotDecoderResults  Visualise CrossDecoder_ValueContext output.
%
% Inputs:
%   Results   struct from CrossDecoder_ValueContext, or path to saved .mat
%
% Optional name-value:
%   'area_label'   ''       e.g. 'FEF' or 'SC' — added to figure title
%   'clim'         [0.4 1]  colour axis limits for T×T matrices
%   'smooth_ms'    0        additional smoothing of timecourses (ms), 0=none
%
% Example:
%   PlotDecoderResults(Results, 'area_label', 'FEF');
%   PlotDecoderResults('Results/decoder_fef.mat', 'area_label', 'FEF');

%% Parse
if ischar(Results)
    d = load(Results);
    Results = d.Results;
end

p = inputParser;
addParameter(p, 'area_label', '',       @ischar);
addParameter(p, 'clim',       [0.4 1],  @isnumeric);
addParameter(p, 'smooth_ms',  0,        @isnumeric);
parse(p, varargin{:});
opts = p.Results;

%% Unpack
ta       = Results.time_axis;
nT       = numel(ta);
dt_ms    = round((ta(2)-ta(1))*1000);
chance   = Results.chance;

acc_wo   = Results.acc_within_obj(:)';
acc_wl   = Results.acc_within_loc(:)';
acc_o2l  = Results.acc_cross_obj2loc(:)';
acc_l2o  = Results.acc_cross_loc2obj(:)';
acc_ct_o = Results.acc_crosstemporal_obj;   % may be stored as acc_ct_obj
acc_ct_l = Results.acc_crosstemporal_loc;
acc_ct_x = Results.acc_crosstemporal_cross;

% Handle both field name conventions
if isfield(Results,'acc_ct_obj')
    acc_ct_o = Results.acc_ct_obj;
    acc_ct_l = Results.acc_ct_loc;
    acc_ct_x = Results.acc_ct_cross;
end

%% Optional smoothing of timecourses
if opts.smooth_ms > 0
    k = opts.smooth_ms / dt_ms;
    win = gausswin(round(6*k), round(6*k)/(2*k));
    win = win / sum(win);
    smth = @(v) filtfilt(win, 1, v(:)');
    acc_wo  = smth(acc_wo);
    acc_wl  = smth(acc_wl);
    acc_o2l = smth(acc_o2l);
    acc_l2o = smth(acc_l2o);
end

%% Colours
c_wo  = [0.15 0.65 0.30];   % green  — within object
c_wl  = [0.80 0.45 0.10];   % orange — within location
c_o2l = [0.15 0.35 0.85];   % blue   — train obj test loc
c_l2o = [0.85 0.15 0.35];   % red    — train loc test obj
lw    = 2;
clim  = opts.clim;

%% Figure
fig_title = 'Cross-Decoder Results';
if ~isempty(opts.area_label)
    fig_title = sprintf('Cross-Decoder Results — %s', opts.area_label);
end
if Results.quick_test
    fig_title = [fig_title, '  [QUICK TEST — not for analysis]'];
end

figure('Name', fig_title, 'Position', [40 40 1500 900], 'Color', 'w');

%% ── Panel A: Same-time accuracy timecourses ──────────────────────────────────
ax1 = subplot(2, 3, 1);
hold on;
plot(ta, acc_wo,  '-',  'Color', c_wo,  'LineWidth', lw, 'DisplayName', 'Within obj');
plot(ta, acc_wl,  '-',  'Color', c_wl,  'LineWidth', lw, 'DisplayName', 'Within loc');
plot(ta, acc_o2l, '--', 'Color', c_o2l, 'LineWidth', lw, 'DisplayName', 'Obj→Loc');
plot(ta, acc_l2o, '--', 'Color', c_l2o, 'LineWidth', lw, 'DisplayName', 'Loc→Obj');
yline(chance, 'k:', 'Chance', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
xline(0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time from stim onset (s)');
ylabel('Decoding accuracy');
title('(A & B) Same-time accuracy');
legend('Location', 'northwest', 'FontSize', 8);
ylim([max(0.3, chance-0.15), 1.02]);
xlim([ta(1) ta(end)]);
box off;

%% ── Panel B: Cross-temporal matrix — within object context ──────────────────
ax2 = subplot(2, 3, 2);
imagesc(ta, ta, acc_ct_o, clim);
colorbar; axis xy square;
hold on;
xline(0, 'w--', 'LineWidth', 1.2);
yline(0, 'w--', 'LineWidth', 1.2);
xlabel('Test time (s)'); ylabel('Train time (s)');
title('(C) Cross-temporal: within object');
colormap(ax2, hot(256));

%% ── Panel C: Cross-temporal matrix — within location context ────────────────
ax3 = subplot(2, 3, 3);
imagesc(ta, ta, acc_ct_l, clim);
colorbar; axis xy square;
hold on;
xline(0, 'w--', 'LineWidth', 1.2);
yline(0, 'w--', 'LineWidth', 1.2);
xlabel('Test time (s)'); ylabel('Train time (s)');
title('(C) Cross-temporal: within location');
colormap(ax3, hot(256));

%% ── Panel D: Cross-temporal cross-context matrix ─────────────────────────────
ax4 = subplot(2, 3, 4);
imagesc(ta, ta, acc_ct_x, clim);
colorbar; axis xy square;
hold on;
xline(0, 'w--', 'LineWidth', 1.2);
yline(0, 'w--', 'LineWidth', 1.2);
xlabel('Test time — location ctx (s)');
ylabel('Train time — object ctx (s)');
title('(D) Cross-temporal cross-context: obj→loc');
colormap(ax4, hot(256));

%% ── Panel E: Diagonals overlaid ──────────────────────────────────────────────
ax5 = subplot(2, 3, 5);
hold on;
d_obj  = diag(acc_ct_o);
d_loc  = diag(acc_ct_l);
d_xctx = diag(acc_ct_x);
plot(ta, d_obj,  '-',  'Color', c_wo,  'LineWidth', lw, 'DisplayName', 'CT within obj');
plot(ta, d_loc,  '-',  'Color', c_wl,  'LineWidth', lw, 'DisplayName', 'CT within loc');
plot(ta, d_xctx, '--', 'Color', c_o2l, 'LineWidth', lw, 'DisplayName', 'CT cross obj→loc');
yline(chance, 'k:', 'Chance', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
xline(0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time (s)');
ylabel('Accuracy (diagonal)');
title('(E) Cross-temporal diagonals');
legend('Location', 'northwest', 'FontSize', 8);
ylim([max(0.3, chance-0.15), 1.02]);
xlim([ta(1) ta(end)]);
box off;

%% ── Panel F: Summary bar — peak accuracy in analysis window ──────────────────
ax6 = subplot(2, 3, 6);
t_win   = ta >= 0.0 & ta <= 0.4;
metrics = [max(acc_wo(t_win)), max(acc_wl(t_win)), ...
           max(acc_o2l(t_win)), max(acc_l2o(t_win)), ...
           max(d_obj(t_win)),   max(d_xctx(t_win))];
labels  = {'Within\nobj','Within\nloc','Obj→Loc','Loc→Obj','CT-diag\nobj','CT-diag\ncross'};
bar_colors = [c_wo; c_wl; c_o2l; c_l2o; c_wo; c_o2l];

b = bar(1:6, metrics, 0.6);
b.FaceColor = 'flat';
b.CData     = bar_colors;
b.EdgeColor = 'none';
yline(chance, 'k:', 'Chance', 'LineWidth', 1.2);
set(gca, 'XTick', 1:6, 'XTickLabel', ...
    {'Within obj','Within loc','Obj→Loc','Loc→Obj','CT obj','CT cross'}, ...
    'XTickLabelRotation', 30);
ylabel('Peak accuracy (0–400ms)');
title('(F) Peak decoding summary');
ylim([0.4, 1.0]);
box off;

sgtitle(fig_title, 'FontSize', 13, 'FontWeight', 'bold');

%% Print summary
t_win2 = ta >= 0.1 & ta <= 0.5;
fprintf('\n══ Decoding Summary (%s) ══\n', opts.area_label);
fprintf('  (A) Within obj  peak: %.3f @ %.0fms\n', ...
    max(acc_wo(t_win2)), ta(argmax_in(acc_wo, t_win2))*1000);
fprintf('  (A) Within loc  peak: %.3f @ %.0fms\n', ...
    max(acc_wl(t_win2)), ta(argmax_in(acc_wl, t_win2))*1000);
fprintf('  (B) Obj→Loc     peak: %.3f @ %.0fms\n', ...
    max(acc_o2l(t_win2)), ta(argmax_in(acc_o2l, t_win2))*1000);
fprintf('  (B) Loc→Obj     peak: %.3f @ %.0fms\n', ...
    max(acc_l2o(t_win2)), ta(argmax_in(acc_l2o, t_win2))*1000);
fprintf('  (E) CT diag obj cross: %.3f\n', max(d_xctx(t_win2)));
fprintf('  Chance: %.3f\n', chance);
fprintf('══════════════════════════════\n');

end  % PlotDecoderResults


%% ── Helper ───────────────────────────────────────────────────────────────────
function idx = argmax_in(v, mask)
v_masked = v;
v_masked(~mask) = -inf;
[~, idx] = max(v_masked);
end
