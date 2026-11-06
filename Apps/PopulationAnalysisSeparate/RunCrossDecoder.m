function Results = RunCrossDecoder(X_obj, y_obj, X_loc, y_loc, time_axis, varargin)
% RunCrossDecoder — MATLAB wrapper for CrossDecoder_ValueContext.py
%
% Saves data to a temp .mat file, calls the Python script, loads results.
%
% Inputs:
%   X_obj      [n_obj_trials × nNeurons × nTimeBins]  object context trials
%   y_obj      [n_obj_trials × 1]   value labels (1=high, 0=low)
%   X_loc      [n_loc_trials × nNeurons × nTimeBins]  location context trials
%   y_loc      [n_loc_trials × 1]
%   time_axis  [1 × nTimeBins]  time in seconds relative to stimulus onset
%
% Optional name-value pairs:
%   'n_iter'            250     resampling iterations
%   'n_jobs'            -1      parallel workers (-1 = all cores)
%   'test_size'         0.2     fraction of trials for test set
%   'skip_crosstemporal' false  set true to skip T×T matrices (quick test)
%   'python_cmd'        'python3'  python executable
%   'script_dir'        fileparts(mfilename('fullpath'))  path to .py file
%   'tmp_dir'           tempdir   where to write intermediate .mat files
%
% Output:
%   Results struct with fields:
%     .acc_within_obj          [nTimeBins]      within-context accuracy
%     .acc_within_loc          [nTimeBins]
%     .acc_cross_obj2loc       [nTimeBins]      train obj → test loc
%     .acc_cross_loc2obj       [nTimeBins]      train loc → test obj
%     .acc_crosstemporal_obj   [nTimeBins × nTimeBins]
%     .acc_crosstemporal_loc   [nTimeBins × nTimeBins]
%     .acc_crosstemporal_cross [nTimeBins × nTimeBins]  train obj t → test loc t'
%     .time_axis               [nTimeBins]
%     .chance                  0.5
%
% Example:
%   Results = RunCrossDecoder(X_obj, y_obj, X_loc, y_loc, time_axis);
%   RunCrossDecoder_PlotResults(Results);  % see companion plot function

%% ── Parse optional arguments ────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'n_iter',             250,               @isnumeric);
addParameter(p, 'n_jobs',             -1,                @isnumeric);
addParameter(p, 'test_size',          0.2,               @isnumeric);
addParameter(p, 'skip_crosstemporal', false,             @islogical);
addParameter(p, 'python_cmd',         'python3',         @ischar);
addParameter(p, 'script_dir',         fileparts(mfilename('fullpath')), @ischar);
addParameter(p, 'tmp_dir',            tempdir,           @ischar);
parse(p, varargin{:});
opts = p.Results;

%% ── Validate inputs ─────────────────────────────────────────────────────────
assert(ndims(X_obj) == 3, 'X_obj must be [nTrials × nNeurons × nTimeBins]');
assert(ndims(X_loc) == 3, 'X_loc must be [nTrials × nNeurons × nTimeBins]');
assert(size(X_obj, 3) == size(X_loc, 3), ...
    'X_obj and X_loc must have the same number of time bins');
assert(numel(time_axis) == size(X_obj, 3), ...
    'time_axis length must equal nTimeBins');

y_obj = y_obj(:);   % ensure column
y_loc = y_loc(:);

fprintf('\n=== RunCrossDecoder ===\n');
fprintf('  Object ctx  : %d trials\n', size(X_obj,1));
fprintf('  Location ctx: %d trials\n', size(X_loc,1));
fprintf('  Neurons     : %d\n',        size(X_obj,2));
fprintf('  Time bins   : %d  (%.0f ms step)\n', size(X_obj,3), ...
        (time_axis(end)-time_axis(1))/(numel(time_axis)-1)*1000);

%% ── Check Python is available ───────────────────────────────────────────────
[status, ~] = system(sprintf('%s --version', opts.python_cmd));
if status ~= 0
    error(['Python not found at "%s".\n' ...
           'Set the ''python_cmd'' parameter, e.g.:\n' ...
           '  RunCrossDecoder(..., ''python_cmd'', ''/usr/bin/python3'')'], ...
          opts.python_cmd);
end

script_path = fullfile(opts.script_dir, 'CrossDecoder_ValueContext.py');
if ~isfile(script_path)
    error('CrossDecoder_ValueContext.py not found at:\n  %s', script_path);
end

%% ── Write input .mat file ───────────────────────────────────────────────────
input_mat  = fullfile(opts.tmp_dir, 'decoder_input.mat');
output_mat = fullfile(opts.tmp_dir, 'decoder_output.mat');

fprintf('\nSaving input data to: %s\n', input_mat);

% Cast to single to reduce file size and match Python float32 expectation
save(input_mat, 'X_obj', 'X_loc', 'y_obj', 'y_loc', 'time_axis', '-v7');

%% ── Build and run Python command ────────────────────────────────────────────
skip_flag = '';
if opts.skip_crosstemporal
    skip_flag = '--skip_crosstemporal';
end

cmd = sprintf(['%s "%s" '...
               '--input "%s" '...
               '--output "%s" '...
               '--n_iter %d '...
               '--n_jobs %d '...
               '--test_size %.3f '...
               '%s'], ...
    opts.python_cmd, script_path, ...
    input_mat, output_mat, ...
    opts.n_iter, opts.n_jobs, opts.test_size, ...
    skip_flag);

fprintf('\nRunning Python decoder...\n');
fprintf('Command: %s\n\n', cmd);

t_start = tic;
[status, cmdout] = system(cmd, '-echo');   % -echo streams output live
elapsed = toc(t_start);

if status ~= 0
    error('Python script failed (exit code %d).\nOutput:\n%s', status, cmdout);
end
fprintf('\nPython finished in %.1f s\n', elapsed);

%% ── Load results ────────────────────────────────────────────────────────────
if ~isfile(output_mat)
    error('Output file not found: %s', output_mat);
end

fprintf('Loading results from: %s\n', output_mat);
raw = load(output_mat);

Results = struct();
Results.acc_within_obj          = raw.acc_within_obj(:)';
Results.acc_within_loc          = raw.acc_within_loc(:)';
Results.acc_cross_obj2loc       = raw.acc_cross_obj2loc(:)';
Results.acc_cross_loc2obj       = raw.acc_cross_loc2obj(:)';
Results.acc_crosstemporal_obj   = raw.acc_crosstemporal_obj;
Results.acc_crosstemporal_loc   = raw.acc_crosstemporal_loc;
Results.acc_crosstemporal_cross = raw.acc_crosstemporal_cross;
Results.time_axis               = raw.time_axis(:)';
Results.chance                  = 0.5;
Results.n_iter                  = opts.n_iter;

%% ── Print summary ────────────────────────────────────────────────────────────
t_win = Results.time_axis >= 0.1 & Results.time_axis <= 0.5;
fprintf('\n--- Decoding summary (peak accuracy, 100–500 ms) ---\n');
fprintf('  Within obj ctx :  %.3f\n',  max(Results.acc_within_obj(t_win)));
fprintf('  Within loc ctx :  %.3f\n',  max(Results.acc_within_loc(t_win)));
fprintf('  Cross obj→loc  :  %.3f\n',  max(Results.acc_cross_obj2loc(t_win)));
fprintf('  Cross loc→obj  :  %.3f\n',  max(Results.acc_cross_loc2obj(t_win)));
fprintf('  Chance level   :  0.500\n');
fprintf('============================\n\n');

%% ── Clean up temp files ─────────────────────────────────────────────────────
delete(input_mat);
delete(output_mat);

end  % function RunCrossDecoder


%% ════════════════════════════════════════════════════════════════════════════
%  COMPANION PLOT FUNCTION  — call after RunCrossDecoder
% ════════════════════════════════════════════════════════════════════════════
function RunCrossDecoder_PlotResults(Results)
% RunCrossDecoder_PlotResults  Visualise cross-decoder output.
%
% Usage:
%   Results = RunCrossDecoder(X_obj, y_obj, X_loc, y_loc, time_axis);
%   RunCrossDecoder_PlotResults(Results);

ta   = Results.time_axis;
lw   = 2;
clim_acc = [0.4 1.0];

c_within_obj  = [0.15 0.65 0.30];  % green
c_within_loc  = [0.80 0.45 0.10];  % orange
c_cross_o2l   = [0.15 0.35 0.85];  % blue
c_cross_l2o   = [0.85 0.15 0.35];  % red

figure('Name', 'Cross-Decoder Results', ...
       'Position', [50 50 1400 850], 'Color', 'w');

% ── Row 1: Same-time accuracy timecourses ─────────────────────────────────
subplot(2, 3, 1);
hold on;
plot(ta, Results.acc_within_obj,    '-',  'Color', c_within_obj, 'LineWidth', lw, ...
     'DisplayName', 'Within obj');
plot(ta, Results.acc_within_loc,    '-',  'Color', c_within_loc, 'LineWidth', lw, ...
     'DisplayName', 'Within loc');
plot(ta, Results.acc_cross_obj2loc, '--', 'Color', c_cross_o2l,  'LineWidth', lw, ...
     'DisplayName', 'Train obj → Test loc');
plot(ta, Results.acc_cross_loc2obj, '--', 'Color', c_cross_l2o,  'LineWidth', lw, ...
     'DisplayName', 'Train loc → Test obj');
yline(0.5, 'k:', 'Chance', 'LineWidth', 1, 'LabelHorizontalAlignment', 'left');
xline(0,   'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time from stim onset (s)');
ylabel('Decoding accuracy');
title('Same-time accuracy');
legend('Location', 'northwest', 'FontSize', 8);
ylim([0.35 1.05]);
box off;

% ── Row 2: Cross-temporal within object context ───────────────────────────
subplot(2, 3, 2);
imagesc(ta, ta, Results.acc_crosstemporal_obj, clim_acc);
colorbar; axis xy;
xline(0, 'w--', 'LineWidth', 1); yline(0, 'w--', 'LineWidth', 1);
xlabel('Test time (s)'); ylabel('Train time (s)');
title('Cross-temporal: within object ctx');
colormap(gca, hot(256));
axis square;

% ── Row 3: Cross-temporal within location context ─────────────────────────
subplot(2, 3, 3);
imagesc(ta, ta, Results.acc_crosstemporal_loc, clim_acc);
colorbar; axis xy;
xline(0, 'w--', 'LineWidth', 1); yline(0, 'w--', 'LineWidth', 1);
xlabel('Test time (s)'); ylabel('Train time (s)');
title('Cross-temporal: within location ctx');
colormap(gca, hot(256));
axis square;

% ── Row 4: Cross-temporal cross-context (train obj → test loc) ────────────
subplot(2, 3, 4);
imagesc(ta, ta, Results.acc_crosstemporal_cross, clim_acc);
colorbar; axis xy;
xline(0, 'w--', 'LineWidth', 1); yline(0, 'w--', 'LineWidth', 1);
xlabel('Test time, location ctx (s)');
ylabel('Train time, object ctx (s)');
title('Cross-temporal cross-context: train obj → test loc');
colormap(gca, hot(256));
axis square;

% ── Row 5: Diagonal comparison (within vs cross-context) ──────────────────
subplot(2, 3, 5);
hold on;
d_obj   = diag(Results.acc_crosstemporal_obj);
d_loc   = diag(Results.acc_crosstemporal_loc);
d_cross = diag(Results.acc_crosstemporal_cross);
plot(ta, d_obj,   '-',  'Color', c_within_obj, 'LineWidth', lw, 'DisplayName', 'CT within obj');
plot(ta, d_loc,   '-',  'Color', c_within_loc, 'LineWidth', lw, 'DisplayName', 'CT within loc');
plot(ta, d_cross, '--', 'Color', c_cross_o2l,  'LineWidth', lw, 'DisplayName', 'CT cross (obj→loc)');
yline(0.5, 'k:', 'Chance', 'LineWidth', 1);
xline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time (s)'); ylabel('Accuracy (diagonal)');
title('Diagonal: cross-temporal accuracy');
legend('Location', 'northwest', 'FontSize', 8);
ylim([0.35 1.05]);
box off;

% ── Row 6: Symmetry check — cross-context obj→loc vs loc→obj diagonal ─────
subplot(2, 3, 6);
hold on;
plot(ta, Results.acc_cross_obj2loc, '-',  'Color', c_cross_o2l, 'LineWidth', lw, ...
     'DisplayName', 'Train obj → Test loc');
plot(ta, Results.acc_cross_loc2obj, '--', 'Color', c_cross_l2o, 'LineWidth', lw, ...
     'DisplayName', 'Train loc → Test obj');
yline(0.5, 'k:', 'Chance', 'LineWidth', 1);
xline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time (s)'); ylabel('Accuracy');
title('Cross-context: generalisation symmetry');
legend('Location', 'northwest', 'FontSize', 8);
ylim([0.35 1.05]);
box off;

sgtitle(sprintf('Cross-Decoder: Value × Context  (n_{iter}=%d)', Results.n_iter), ...
        'FontSize', 13, 'FontWeight', 'bold');

end  % function RunCrossDecoder_PlotResults

%% ════════════════════════════════════════════════════════════════════════════
%  ORTHOGONALITY ENTRYPOINT — Value axis angle (Object vs Location)
% ════════════════════════════════════════════════════════════════════════════
function OrthResults = RunValueAxisOrthogonality( ...
    psth_high_value_object_RF, psth_low_value_object_RF, ...
    psth_high_value_location_RF, psth_low_value_location_RF, time_axis, varargin)
% RunValueAxisOrthogonality
% Convenience wrapper for ComputeValueAxisOrthogonality.m using the same
% PopulationAnalysisSeparate workflow style.
fprintf('\n=== RunValueAxisOrthogonality ===\n');
OrthResults = ComputeValueAxisOrthogonality( ...
    psth_high_value_object_RF, psth_low_value_object_RF, ...
    psth_high_value_location_RF, psth_low_value_location_RF, ...
    time_axis, varargin{:});
fprintf('=================================\n\n');
end

%% ════════════════════════════════════════════════════════════════════════════
%  USAGE EXAMPLE BLOCK — RunCrossDecoder_Example
% ════════════════════════════════════════════════════════════════════════════
%{
%% Example: build X/y from your PSTH arrays and run the decoder
%
%  Assume you have (from Batch_ObjLocValue_Organize or similar):
%    psth_obj: [nTrials_obj × nNeurons × nTimeBins]
%    labels_obj: [nTrials_obj × 1]  (1=high value, 0=low value)
%    psth_loc: [nTrials_loc × nNeurons × nTimeBins]
%    labels_loc: [nTrials_loc × 1]
%    time_axis: [1 × nTimeBins]  seconds

% Quick test (skip T×T matrices, fast):
Results = RunCrossDecoder(psth_obj, labels_obj, psth_loc, labels_loc, time_axis, ...
    'n_iter', 50, 'skip_crosstemporal', true);

% Full run:
Results = RunCrossDecoder(psth_obj, labels_obj, psth_loc, labels_loc, time_axis, ...
    'n_iter', 250);

% Plot:
RunCrossDecoder_PlotResults(Results);

% Save:
save('CrossDecoder_Results.mat', 'Results');
%
% Orthogonality analysis (mean-difference value axis):
% OrthResults = RunValueAxisOrthogonality( ...
%     psth_high_value_object_RF, psth_low_value_object_RF, ...
%     psth_high_value_location_RF, psth_low_value_location_RF, ...
%     time_axis, 'window_ms', [100 200], 'n_repeats', 200);
%}
