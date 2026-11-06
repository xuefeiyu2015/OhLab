function Results = CrossDecoder_ValueContext(X_obj, y_obj, X_loc, y_loc, time_axis, varargin)
% CrossDecoder_ValueContext  Cross-context & cross-temporal value decoder.
%
% Implements four decodings of high vs low value using logistic regression,
% following Tafazoli et al. (2026) Nature 650:164-172.
%
% Inputs:
%   X_obj      [nTrials_obj x nNeurons x nTimeBins]  object context, RF trials
%   y_obj      [nTrials_obj x 1]  labels: 1=high value, 0=low value
%   X_loc      [nTrials_loc x nNeurons x nTimeBins]  location context, RF trials
%   y_loc      [nTrials_loc x 1]
%   time_axis  [1 x nTimeBins]  seconds re stimulus onset
%
% Optional name-value pairs:
%   'n_folds'        10       number of CV folds
%   'C'              60       LR regularisation (lambda = 1/C)
%   'quick_test'     false    run 2 folds x 20 time bins only
%   'save_path'      ''       full path for results .mat ('' = no save)
%   'ckpt_dir'       ''       checkpoint directory ('' = auto in save_path dir)
%   'use_parallel'   true     use parfor if PCT available
%   'seed'           42       random seed for reproducibility
%
% Output:
%   Results struct (also saved to save_path if provided):
%     .acc_within_obj        [nT x 1]     decoding A, object context
%     .acc_within_loc        [nT x 1]     decoding A, location context
%     .acc_cross_obj2loc     [nT x 1]     decoding B, train obj test loc
%     .acc_cross_loc2obj     [nT x 1]     decoding B, train loc test obj
%     .acc_ct_obj            [nT x nT]    decoding C, cross-temporal object
%     .acc_ct_loc            [nT x nT]    decoding C, cross-temporal location
%     .acc_ct_cross          [nT x nT]    decoding D, cross-temporal cross-ctx
%     .weights_A_obj         [nT x nF x nN]  fold weights for A, object
%     .weights_A_loc         [nT x nF x nN]
%     .weights_B_obj2loc     [nT x nF x nN]
%     .weights_B_loc2obj     [nT x nF x nN]
%     .weights_CT_obj        [nT x nF x nN]  weights at each t_train, C obj
%     .weights_CT_loc        [nT x nF x nN]
%     .weights_CT_cross      [nT x nF x nN]
%     .time_axis, .chance, .n_folds, .n_neurons, .n_time, .quick_test
%     .run_timestamp
%
% Usage:
%   Results = CrossDecoder_ValueContext(X_obj, y_obj, X_loc, y_loc, time_axis,...
%       'save_path', 'Results/decoder_fef.mat', 'quick_test', true);

%% ── Parse options ────────────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'n_folds',      10,    @isnumeric);
addParameter(p, 'C',            60,    @isnumeric);
addParameter(p, 'quick_test',   false, @islogical);
addParameter(p, 'save_path',    '',    @ischar);
addParameter(p, 'ckpt_dir',     '',    @ischar);
addParameter(p, 'use_parallel', true,  @islogical);
addParameter(p, 'seed',         42,    @isnumeric);
parse(p, varargin{:});
opts = p.Results;

rng(opts.seed);
lambda = 1 / opts.C;

%% ── Validate inputs ──────────────────────────────────────────────────────────
assert(ndims(X_obj)==3 && ndims(X_loc)==3, ...
    'X_obj and X_loc must be [trials x neurons x timebins]');
assert(size(X_obj,3)==size(X_loc,3), 'X_obj and X_loc must have same nTimeBins');
assert(size(X_obj,2)==size(X_loc,2), 'X_obj and X_loc must have same nNeurons');
assert(numel(time_axis)==size(X_obj,3), 'time_axis length must equal nTimeBins');

y_obj = y_obj(:);
y_loc = y_loc(:);
assert(all(ismember(y_obj,[0 1])) && all(ismember(y_loc,[0 1])), ...
    'Labels must be 0 (low value) or 1 (high value)');

[~, nNeurons, nTimeBins] = size(X_obj);
nFolds = opts.n_folds;

%% ── Quick-test mode ──────────────────────────────────────────────────────────
if opts.quick_test
    nFolds  = 2;
    % Use first 20 and last 20 bins to cover pre-stim and response period
    qt_bins = [1:20, nTimeBins-19:nTimeBins];
    qt_bins = unique(qt_bins);
    X_obj   = X_obj(:, :, qt_bins);
    X_loc   = X_loc(:, :, qt_bins);
    time_axis = time_axis(qt_bins);
    nTimeBins = numel(qt_bins);
    log_msg('⚠  QUICK TEST MODE: 2 folds, %d time bins. Results not valid for analysis.', nTimeBins);
end

%% ── Setup: parallel, checkpoint dir ─────────────────────────────────────────
% Check Parallel Computing Toolbox
pct_available = license('test','Distrib_Computing_Toolbox') && ...
                ~isempty(ver('parallel'));
use_parfor    = opts.use_parallel && pct_available;

if opts.use_parallel && ~pct_available
    log_msg('Parallel Computing Toolbox not found — using serial loops.');
end
if use_parfor
    log_msg('Parallel Computing Toolbox detected — using parfor.');
    % Open pool if not already open
    if isempty(gcp('nocreate'))
        parpool('local');
    end
end

% Checkpoint directory
if isempty(opts.ckpt_dir)
    if ~isempty(opts.save_path)
        opts.ckpt_dir = fullfile(fileparts(opts.save_path), 'checkpoints');
    else
        opts.ckpt_dir = fullfile(tempdir, 'decoder_checkpoints');
    end
end
if ~exist(opts.ckpt_dir, 'dir'), mkdir(opts.ckpt_dir); end

log_msg('Checkpoint dir: %s', opts.ckpt_dir);
log_msg('Neurons=%d | TimeBins=%d | Folds=%d | Lambda=%.4f', ...
    nNeurons, nTimeBins, nFolds, lambda);

%% ── Build stratified CV folds ────────────────────────────────────────────────
log_msg('Building %d-fold stratified CV partitions...', nFolds);

cv_obj = cvpartition(y_obj, 'KFold', nFolds, 'Stratify', true);
cv_loc = cvpartition(y_loc, 'KFold', nFolds, 'Stratify', true);

% Pre-compute fold indices as cell arrays (needed for parfor broadcast)
train_idx_obj = cell(nFolds, 1);
test_idx_obj  = cell(nFolds, 1);
train_idx_loc = cell(nFolds, 1);
test_idx_loc  = cell(nFolds, 1);
for f = 1:nFolds
    train_idx_obj{f} = find(cv_obj.training(f));
    test_idx_obj{f}  = find(cv_obj.test(f));
    train_idx_loc{f} = find(cv_loc.training(f));
    test_idx_loc{f}  = find(cv_loc.test(f));
end

%% ── Build balanced cross-context test set ────────────────────────────────────
% Used for decodings B and D (fixed test set across all folds)
log_msg('Building balanced cross-context test sets...');

[X_test_cross_obj, y_test_cross_obj] = balance_testset(X_obj, y_obj, opts.seed);
[X_test_cross_loc, y_test_cross_loc] = balance_testset(X_loc, y_loc, opts.seed);

n_test_cross_loc = size(X_test_cross_loc, 1);
n_test_cross_obj = size(X_test_cross_obj, 1);

log_msg('Cross-ctx test sets: loc=%d trials, obj=%d trials', ...
    n_test_cross_loc, n_test_cross_obj);

%% ══════════════════════════════════════════════════════════════════════════════
%  DECODING A — Within-context same-time
%  Stratified 10-fold CV within each context separately
% ══════════════════════════════════════════════════════════════════════════════
ckpt_A = fullfile(opts.ckpt_dir, 'ckpt_A.mat');

if exist(ckpt_A, 'file')
    log_msg('Checkpoint found: loading Decoding A results...');
    c = load(ckpt_A);
    acc_within_obj  = c.acc_within_obj;
    acc_within_loc  = c.acc_within_loc;
    weights_A_obj   = c.weights_A_obj;
    weights_A_loc   = c.weights_A_loc;
else
    log_msg('══ Decoding A: Within-context same-time ══');
    t_A = tic;

    [acc_within_obj, weights_A_obj] = decode_same_time(X_obj, y_obj, ...
        X_obj, y_obj, train_idx_obj, test_idx_obj, ...
        nFolds, nNeurons, nTimeBins, lambda, use_parfor, 'A_obj');

    [acc_within_loc, weights_A_loc] = decode_same_time(X_loc, y_loc, ...
        X_loc, y_loc, train_idx_loc, test_idx_loc, ...
        nFolds, nNeurons, nTimeBins, lambda, use_parfor, 'A_loc');

    log_msg('Decoding A done (%.1fs) | peak obj=%.3f @ %.0fms | peak loc=%.3f @ %.0fms', ...
        toc(t_A), max(acc_within_obj), time_axis(argmax(acc_within_obj))*1000, ...
        max(acc_within_loc), time_axis(argmax(acc_within_loc))*1000);

    save(ckpt_A, 'acc_within_obj','acc_within_loc','weights_A_obj','weights_A_loc');
end

%% ══════════════════════════════════════════════════════════════════════════════
%  DECODING B — Cross-context same-time
%  Train on context A (10-fold), test each fold on balanced context B trials
% ══════════════════════════════════════════════════════════════════════════════
ckpt_B = fullfile(opts.ckpt_dir, 'ckpt_B.mat');

if exist(ckpt_B, 'file')
    log_msg('Checkpoint found: loading Decoding B results...');
    c = load(ckpt_B);
    acc_cross_obj2loc  = c.acc_cross_obj2loc;
    acc_cross_loc2obj  = c.acc_cross_loc2obj;
    weights_B_obj2loc  = c.weights_B_obj2loc;
    weights_B_loc2obj  = c.weights_B_loc2obj;
else
    log_msg('══ Decoding B: Cross-context same-time ══');
    t_B = tic;

    % Train obj, test loc
    [acc_cross_obj2loc, weights_B_obj2loc] = decode_same_time(X_obj, y_obj, ...
        X_test_cross_loc, y_test_cross_loc, train_idx_obj, [], ...
        nFolds, nNeurons, nTimeBins, lambda, use_parfor, 'B_obj2loc');

    % Train loc, test obj
    [acc_cross_loc2obj, weights_B_loc2obj] = decode_same_time(X_loc, y_loc, ...
        X_test_cross_obj, y_test_cross_obj, train_idx_loc, [], ...
        nFolds, nNeurons, nTimeBins, lambda, use_parfor, 'B_loc2obj');

    log_msg('Decoding B done (%.1fs) | peak obj->loc=%.3f @ %.0fms | peak loc->obj=%.3f @ %.0fms', ...
        toc(t_B), max(acc_cross_obj2loc), time_axis(argmax(acc_cross_obj2loc))*1000, ...
        max(acc_cross_loc2obj), time_axis(argmax(acc_cross_loc2obj))*1000);

    save(ckpt_B, 'acc_cross_obj2loc','acc_cross_loc2obj', ...
        'weights_B_obj2loc','weights_B_loc2obj');
end

%% ══════════════════════════════════════════════════════════════════════════════
%  DECODING C — Within-context cross-temporal (object context)
%  For each t_train: train at t_train, predict at ALL t_test (vectorised)
% ══════════════════════════════════════════════════════════════════════════════
ckpt_C_obj = fullfile(opts.ckpt_dir, 'ckpt_C_obj.mat');

if exist(ckpt_C_obj, 'file')
    log_msg('Checkpoint found: resuming Decoding C (object)...');
    c = load(ckpt_C_obj);
    acc_ct_obj     = c.acc_ct_obj;
    weights_CT_obj = c.weights_CT_obj;
    t_resume_C_obj = c.last_t + 1;
else
    acc_ct_obj     = zeros(nTimeBins, nTimeBins, 'single');
    weights_CT_obj = zeros(nTimeBins, nFolds, nNeurons, 'single');
    t_resume_C_obj = 1;
end

if t_resume_C_obj <= nTimeBins
    log_msg('══ Decoding C: Within-context cross-temporal (object) ══');
    log_msg('Starting from t_train=%d/%d', t_resume_C_obj, nTimeBins);
    t_C = tic;

    for t = t_resume_C_obj : nTimeBins
        [acc_row, w_row] = decode_crosstemporal_one_tbin(...
            X_obj, y_obj, X_obj, y_obj, ...
            train_idx_obj, test_idx_obj, ...
            t, nFolds, nNeurons, nTimeBins, lambda, use_parfor, false);

        acc_ct_obj(t, :)        = acc_row;
        weights_CT_obj(t, :, :) = w_row;   % [nFolds x nNeurons]

        log_msg('C_obj  t_train=%3d/%d | fold-avg diag acc=%.3f | elapsed %s', ...
            t, nTimeBins, acc_ct_obj(t,t), fmt_time(t_C));

        last_t = t;
        save(ckpt_C_obj, 'acc_ct_obj', 'weights_CT_obj', 'last_t');
    end

    log_msg('Decoding C (object) done (%s)', fmt_time(t_C));
end

%% ── Decoding C — Within-context cross-temporal (location context) ────────────
ckpt_C_loc = fullfile(opts.ckpt_dir, 'ckpt_C_loc.mat');

if exist(ckpt_C_loc, 'file')
    log_msg('Checkpoint found: resuming Decoding C (location)...');
    c = load(ckpt_C_loc);
    acc_ct_loc     = c.acc_ct_loc;
    weights_CT_loc = c.weights_CT_loc;
    t_resume_C_loc = c.last_t + 1;
else
    acc_ct_loc     = zeros(nTimeBins, nTimeBins, 'single');
    weights_CT_loc = zeros(nTimeBins, nFolds, nNeurons, 'single');
    t_resume_C_loc = 1;
end

if t_resume_C_loc <= nTimeBins
    log_msg('══ Decoding C: Within-context cross-temporal (location) ══');
    log_msg('Starting from t_train=%d/%d', t_resume_C_loc, nTimeBins);
    t_C2 = tic;

    for t = t_resume_C_loc : nTimeBins
        [acc_row, w_row] = decode_crosstemporal_one_tbin(...
            X_loc, y_loc, X_loc, y_loc, ...
            train_idx_loc, test_idx_loc, ...
            t, nFolds, nNeurons, nTimeBins, lambda, use_parfor, false);

        acc_ct_loc(t, :)        = acc_row;
        weights_CT_loc(t, :, :) = w_row;

        log_msg('C_loc  t_train=%3d/%d | fold-avg diag acc=%.3f | elapsed %s', ...
            t, nTimeBins, acc_ct_loc(t,t), fmt_time(t_C2));

        last_t = t;
        save(ckpt_C_loc, 'acc_ct_loc', 'weights_CT_loc', 'last_t');
    end

    log_msg('Decoding C (location) done (%s)', fmt_time(t_C2));
end

%% ══════════════════════════════════════════════════════════════════════════════
%  DECODING D — Cross-context cross-temporal
%  Train on object at t_train, predict at ALL t_test on balanced location trials
% ══════════════════════════════════════════════════════════════════════════════
ckpt_D = fullfile(opts.ckpt_dir, 'ckpt_D.mat');

if exist(ckpt_D, 'file')
    log_msg('Checkpoint found: resuming Decoding D...');
    c = load(ckpt_D);
    acc_ct_cross     = c.acc_ct_cross;
    weights_CT_cross = c.weights_CT_cross;
    t_resume_D       = c.last_t + 1;
else
    acc_ct_cross     = zeros(nTimeBins, nTimeBins, 'single');
    weights_CT_cross = zeros(nTimeBins, nFolds, nNeurons, 'single');
    t_resume_D       = 1;
end

if t_resume_D <= nTimeBins
    log_msg('══ Decoding D: Cross-context cross-temporal (train obj, test loc) ══');
    log_msg('Starting from t_train=%d/%d', t_resume_D, nTimeBins);
    t_D = tic;

    for t = t_resume_D : nTimeBins
        [acc_row, w_row] = decode_crosstemporal_one_tbin(...
            X_obj, y_obj, X_test_cross_loc, y_test_cross_loc, ...
            train_idx_obj, [], ...
            t, nFolds, nNeurons, nTimeBins, lambda, use_parfor, true);

        acc_ct_cross(t, :)        = acc_row;
        weights_CT_cross(t, :, :) = w_row;

        log_msg('D_cross t_train=%3d/%d | fold-avg acc at t=0: %.3f | elapsed %s', ...
            t, nTimeBins, mean(acc_ct_cross(t,:)), fmt_time(t_D));

        last_t = t;
        save(ckpt_D, 'acc_ct_cross', 'weights_CT_cross', 'last_t');
    end

    log_msg('Decoding D done (%s)', fmt_time(t_D));
end

%% ── Pack Results ─────────────────────────────────────────────────────────────
log_msg('Packing results...');

Results = struct();
Results.acc_within_obj      = acc_within_obj;
Results.acc_within_loc      = acc_within_loc;
Results.acc_cross_obj2loc   = acc_cross_obj2loc;
Results.acc_cross_loc2obj   = acc_cross_loc2obj;
Results.acc_ct_obj          = acc_ct_obj;
Results.acc_ct_loc          = acc_ct_loc;
Results.acc_ct_cross        = acc_ct_cross;
Results.weights_A_obj       = weights_A_obj;
Results.weights_A_loc       = weights_A_loc;
Results.weights_B_obj2loc   = weights_B_obj2loc;
Results.weights_B_loc2obj   = weights_B_loc2obj;
Results.weights_CT_obj      = weights_CT_obj;
Results.weights_CT_loc      = weights_CT_loc;
Results.weights_CT_cross    = weights_CT_cross;
Results.time_axis           = time_axis;
Results.chance              = 0.5;
Results.n_folds             = nFolds;
Results.n_neurons           = nNeurons;
Results.n_time              = nTimeBins;
Results.quick_test          = opts.quick_test;
Results.run_timestamp       = datestr(now);

%% ── Save ─────────────────────────────────────────────────────────────────────
if ~isempty(opts.save_path)
    save_dir = fileparts(opts.save_path);
    if ~isempty(save_dir) && ~exist(save_dir,'dir')
        mkdir(save_dir);
    end
    save(opts.save_path, 'Results', '-v7.3');
    log_msg('Results saved to: %s', opts.save_path);
end

%% ── Clean up checkpoints ─────────────────────────────────────────────────────
ckpt_files = {ckpt_A, ckpt_B, ckpt_C_obj, ckpt_C_loc, ckpt_D};
for k = 1:numel(ckpt_files)
    if exist(ckpt_files{k}, 'file'), delete(ckpt_files{k}); end
end
log_msg('Checkpoint files cleaned up.');
log_msg('══ CrossDecoder_ValueContext complete ══');

end  % main function


%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTION: decode_same_time
%  Decodings A and B — fit at each time bin, parfor over time bins
% ════════════════════════════════════════════════════════════════════════════
function [acc_vec, weights_out] = decode_same_time(...
    X_train_ctx, y_train_ctx, X_test_ctx, y_test_ctx, ...
    train_idx, test_idx, ...
    nFolds, nNeurons, nTimeBins, lambda, use_parfor, label)
% For each time bin: 10-fold CV (within-context) or fixed test set (cross-context)
% test_idx = [] signals cross-context mode (use full X_test_ctx as test set)

cross_ctx_mode = isempty(test_idx);

acc_vec    = zeros(nTimeBins, 1, 'single');
weights_out = zeros(nTimeBins, nFolds, nNeurons, 'single');

if use_parfor
    % parfor over time bins (embarrassingly parallel at this level)
    parfor t = 1:nTimeBins
        [a, w] = decode_one_timebin(X_train_ctx, y_train_ctx, ...
            X_test_ctx, y_test_ctx, ...
            train_idx, test_idx, ...
            t, nFolds, nNeurons, lambda, cross_ctx_mode);
        acc_vec(t)           = a;
        weights_out(t, :, :) = w;   % [nFolds x nNeurons]
    end
else
    for t = 1:nTimeBins
        [a, w] = decode_one_timebin(X_train_ctx, y_train_ctx, ...
            X_test_ctx, y_test_ctx, ...
            train_idx, test_idx, ...
            t, nFolds, nNeurons, lambda, cross_ctx_mode);
        acc_vec(t)           = a;
        weights_out(t, :, :) = w;

        if mod(t, 10) == 0
            log_msg('  [%s] t=%d/%d done', label, t, nTimeBins);
        end
    end
end
end  % decode_same_time


%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTION: decode_one_timebin
%  Fit + evaluate at a single time bin across all folds
% ════════════════════════════════════════════════════════════════════════════
function [acc_mean, weights_folds] = decode_one_timebin(...
    X_train_ctx, y_train_ctx, X_test_ctx, y_test_ctx, ...
    train_idx, test_idx, t, nFolds, nNeurons, lambda, cross_ctx_mode)

fold_acc     = zeros(nFolds, 1);
fold_weights = zeros(nNeurons, nFolds);

for f = 1:nFolds
    tr_idx = train_idx{f};

    X_tr = squeeze(X_train_ctx(tr_idx, :, t));   % [n_train x nNeurons]
    y_tr = y_train_ctx(tr_idx);

    if cross_ctx_mode
        % Test on full balanced cross-context set (fixed, not fold-split)
        X_te = squeeze(X_test_ctx(:, :, t));      % [n_test x nNeurons]
        y_te = y_test_ctx;
    else
        te_idx = test_idx{f};
        X_te   = squeeze(X_train_ctx(te_idx, :, t));  % within-ctx test fold
        y_te   = y_train_ctx(te_idx);
    end

    % Z-score using training statistics
    [X_tr_z, mu, sd] = zscore_fit(X_tr);
    X_te_z           = zscore_apply(X_te, mu, sd);

    % Fit logistic regression
    Mdl = fit_lr(X_tr_z, y_tr, lambda);

    % Predict
    scores = X_te_z * Mdl.Beta + Mdl.Bias;
    preds  = double(scores > 0);
    fold_acc(f)      = mean(preds == y_te);
    fold_weights(:,f) = Mdl.Beta;
end

acc_mean     = mean(fold_acc);
weights_folds = fold_weights';   % [nFolds x nNeurons]
end  % decode_one_timebin


%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTION: decode_crosstemporal_one_tbin
%  For one t_train: parfor over folds, vectorised predict over ALL t_test
% ════════════════════════════════════════════════════════════════════════════
function [acc_row, weights_row] = decode_crosstemporal_one_tbin(...
    X_train_ctx, y_train_ctx, X_test_ctx, y_test_ctx, ...
    train_idx, test_idx, ...
    t_train, nFolds, nNeurons, nTimeBins, lambda, use_parfor, cross_ctx_mode)
% Returns:
%   acc_row     [1 x nTimeBins]  accuracy at every t_test
%   weights_row [nFolds x nNeurons]

fold_results = cell(nFolds, 1);

if use_parfor
    parfor f = 1:nFolds
        fold_results{f} = ct_one_fold(X_train_ctx, y_train_ctx, ...
            X_test_ctx, y_test_ctx, ...
            train_idx, test_idx, ...
            t_train, f, nNeurons, nTimeBins, lambda, cross_ctx_mode);
    end
else
    for f = 1:nFolds
        fold_results{f} = ct_one_fold(X_train_ctx, y_train_ctx, ...
            X_test_ctx, y_test_ctx, ...
            train_idx, test_idx, ...
            t_train, f, nNeurons, nTimeBins, lambda, cross_ctx_mode);
    end
end

% Aggregate across folds
acc_row     = zeros(1, nTimeBins, 'single');
weights_row = zeros(nFolds, nNeurons, 'single');
for f = 1:nFolds
    acc_row     = acc_row + fold_results{f}.acc / nFolds;
    weights_row(f, :) = fold_results{f}.w;
end
end  % decode_crosstemporal_one_tbin


%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTION: ct_one_fold
%  One fold of cross-temporal decoding — the innermost unit
% ════════════════════════════════════════════════════════════════════════════
function result = ct_one_fold(...
    X_train_ctx, y_train_ctx, X_test_ctx, y_test_ctx, ...
    train_idx, test_idx, ...
    t_train, fold, nNeurons, nTimeBins, lambda, cross_ctx_mode)

tr_idx = train_idx{fold};
X_tr   = squeeze(X_train_ctx(tr_idx, :, t_train));   % [n_train x nNeurons]
y_tr   = y_train_ctx(tr_idx);

% Z-score stats from training data at t_train
[X_tr_z, mu, sd] = zscore_fit(X_tr);

% Fit classifier at t_train
Mdl = fit_lr(X_tr_z, y_tr, lambda);
w   = Mdl.Beta;    % [nNeurons x 1]
b   = Mdl.Bias;    % scalar

% Build test set
if cross_ctx_mode
    X_te_all = X_test_ctx;             % [n_test x nNeurons x nTimeBins]
    y_te     = y_test_ctx;
else
    te_idx   = test_idx{fold};
    X_te_all = X_train_ctx(te_idx,:,:);  % [n_test x nNeurons x nTimeBins]
    y_te     = y_train_ctx(te_idx);
end
n_test = size(X_te_all, 1);

% Apply training z-score stats to ALL test time bins
% mu, sd: [1 x nNeurons] -> broadcast over [n_test x nNeurons x nTimeBins]
mu_3d = reshape(mu, 1, nNeurons, 1);
sd_3d = reshape(sd, 1, nNeurons, 1);
X_te_z = (X_te_all - mu_3d) ./ sd_3d;   % [n_test x nNeurons x nTimeBins]

% ── KEY SPEED TRICK: vectorised prediction over all t_test ────────────────
% Reshape to [n_test*nTimeBins x nNeurons], multiply by w, reshape back
X_2d       = reshape(permute(X_te_z, [1 3 2]), [], nNeurons); % [n_test*nT x nN]
scores_all = X_2d * w + b;                                      % [n_test*nT x 1]
preds_all  = reshape(scores_all > 0, n_test, nTimeBins);       % [n_test x nT]
acc_vec    = mean(preds_all == y_te, 1);                       % [1 x nTimeBins]

result.acc = single(acc_vec);
result.w   = single(w');          % [1 x nNeurons]
end  % ct_one_fold


%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTION: balance_testset
%  Draw balanced (equal high/low) cross-context test set
% ════════════════════════════════════════════════════════════════════════════
function [X_bal, y_bal] = balance_testset(X, y, seed)
rng(seed + 1);
idx_high = find(y == 1);
idx_low  = find(y == 0);
n_each   = min(numel(idx_high), numel(idx_low));
idx_high = idx_high(randperm(numel(idx_high), n_each));
idx_low  = idx_low(randperm(numel(idx_low),   n_each));
idx_bal  = [idx_high; idx_low];
X_bal    = X(idx_bal, :, :);
y_bal    = y(idx_bal);
end  % balance_testset


%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTION: fit_lr
%  Fit logistic regression with L2 regularisation
% ════════════════════════════════════════════════════════════════════════════
function Mdl = fit_lr(X, y, lambda)
% fitclinear: objective = (1/N) sum(loss) + lambda * ||w||^2
% Equivalent to L2 logistic regression with regularisation strength lambda
Mdl = fitclinear(X, y, ...
    'Learner',        'logistic', ...
    'Regularization', 'ridge', ...
    'Lambda',          lambda, ...
    'Solver',         'lbfgs', ...
    'Verbose',         0);
end  % fit_lr


%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTION: zscore_fit
%  Compute mean/std from X and return normalised copy
% ════════════════════════════════════════════════════════════════════════════
function [X_z, mu, sd] = zscore_fit(X)
mu   = mean(X, 1);                % [1 x nNeurons]
sd   = std(X, 0, 1) + 1e-8;      % [1 x nNeurons], epsilon avoids div/0
X_z  = (X - mu) ./ sd;
end


%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTION: zscore_apply
%  Apply pre-computed mean/std to new data
% ════════════════════════════════════════════════════════════════════════════
function X_z = zscore_apply(X, mu, sd)
X_z = (X - mu) ./ sd;
end


%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTION: log_msg
% ════════════════════════════════════════════════════════════════════════════
function log_msg(fmt, varargin)
msg = sprintf(fmt, varargin{:});
fprintf('[CrossDecoder %s] %s\n', datestr(now,'HH:MM:SS'), msg);
end


%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTION: fmt_time
% ════════════════════════════════════════════════════════════════════════════
function s = fmt_time(t_start)
e = toc(t_start);
h = floor(e/3600);
m = floor(mod(e,3600)/60);
s_val = floor(mod(e,60));
s = sprintf('%02d:%02d:%02d', h, m, s_val);
end


%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTION: argmax
% ════════════════════════════════════════════════════════════════════════════
function idx = argmax(v)
[~, idx] = max(v(:));
end
