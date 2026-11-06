"""
CrossDecoder_ValueContext.py
============================
Cross-context and cross-temporal logistic regression decoding of value
(high vs low) for object and location task contexts.

Based on Tafazoli et al. (2026) Nature 650:164-172.
  - Classifier: logistic regression, L2 regularisation (C=60, i.e. λ=1/60)
  - 250 resampling iterations with balanced train/test splits
  - KEY SPEED TRICK: after fitting at t_train, predict at ALL t_test
    simultaneously via a single matrix multiply → O(T) fits not O(T²)
  - Parallelised over t_train using joblib

Input  (.mat file written by RunCrossDecoder.m):
    X_obj  [n_obj_trials  × nNeurons × nTimeBins]  float32
    X_loc  [n_loc_trials  × nNeurons × nTimeBins]  float32
    y_obj  [n_obj_trials]   int   (1=high value, 0=low value)
    y_loc  [n_loc_trials]   int
    time_axis  [nTimeBins]  float  (seconds)

Output (.mat file loaded back by RunCrossDecoder.m):
    acc_within_obj   [nTimeBins]        within-context accuracy, object
    acc_within_loc   [nTimeBins]        within-context accuracy, location
    acc_cross_obj2loc [nTimeBins]       train obj → test loc
    acc_cross_loc2obj [nTimeBins]       train loc → test obj
    acc_crosstemporal_obj  [nTimeBins × nTimeBins]  train t → test t', object ctx
    acc_crosstemporal_loc  [nTimeBins × nTimeBins]  train t → test t', location ctx
    acc_crosstemporal_cross [nTimeBins × nTimeBins] train obj at t → test loc at t'
    time_axis [nTimeBins]
    chance    scalar (0.5)

Usage:
    python CrossDecoder_ValueContext.py --input decoder_input.mat \
                                        --output decoder_output.mat \
                                        --n_iter 250 \
                                        --n_jobs -1

Google Colab usage:
    !python CrossDecoder_ValueContext.py --input /content/decoder_input.mat \
                                         --output /content/decoder_output.mat
"""

import argparse
import time
import numpy as np
import scipy.io as sio
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedShuffleSplit
from joblib import Parallel, delayed

# ── Classifier factory ────────────────────────────────────────────────────────
def make_clf():
    """Logistic regression matching Tafazoli et al.: L2, λ=1/60 → C=60."""
    return LogisticRegression(
        C=60.0,
        penalty='l2',
        solver='lbfgs',
        max_iter=1000,
        tol=1e-4,
    )


# ── Core: train at t_train, predict at ALL t_test (vectorised) ───────────────
def _fit_and_predict_all_times(X_train_t, y_train, X_test_all, y_test):
    """
    Train a classifier on X_train_t [n_train × nNeurons] at one time bin,
    then predict at every time bin in X_test_all [n_test × nNeurons × nTimeBins].

    Returns accuracy vector of shape [nTimeBins].

    Speed trick: after fitting, decision score for all test times is
        scores = X_test_all @ w + b   →  [n_test × nTimeBins]
    which is a single BLAS dgemm call instead of T separate predict() calls.
    """
    clf = make_clf()
    clf.fit(X_train_t, y_train)

    w = clf.coef_[0]        # [nNeurons]
    b = clf.intercept_[0]   # scalar

    # X_test_all: [n_test × nNeurons × nTimeBins]
    # scores_all: [n_test × nTimeBins]
    scores_all = np.einsum('ijk,j->ik', X_test_all, w) + b
    preds_all  = (scores_all > 0).astype(np.int8)  # [n_test × nTimeBins]

    acc = np.mean(preds_all == y_test[:, np.newaxis], axis=0)  # [nTimeBins]
    return acc


# ── One resampling iteration ──────────────────────────────────────────────────
def _one_iteration(X_train_ctx, y_train_ctx,
                   X_test_ctx,  y_test_ctx,
                   train_idx, test_idx,
                   n_time):
    """
    For a single train/test split, fit classifiers at every t_train and
    collect accuracy at all t_test. Returns [nTimeBins × nTimeBins].
    """
    X_tr = X_train_ctx[train_idx]   # [n_tr × nNeurons × nTimeBins]
    y_tr = y_train_ctx[train_idx]

    X_te = X_test_ctx[test_idx]     # [n_te × nNeurons × nTimeBins]
    y_te = y_test_ctx[test_idx]

    acc_mat = np.zeros((n_time, n_time), dtype=np.float32)

    for t in range(n_time):
        acc_mat[t, :] = _fit_and_predict_all_times(
            X_tr[:, :, t], y_tr, X_te, y_te
        )
    return acc_mat


# ── Cross-temporal decoding (within one context or cross-context) ─────────────
def cross_temporal_decode(X_train_ctx, y_train_ctx,
                          X_test_ctx,  y_test_ctx,
                          n_iter=250, test_size=0.2, n_jobs=-1, seed=0):
    """
    Train at every t_train, test at every t_test.
    Returns mean accuracy matrix [nTimeBins × nTimeBins] across iterations.

    For within-context: pass same X/y for train and test ctx (split internally).
    For cross-context:  pass different X/y; test_size applies to test ctx only.
    """
    n_time = X_train_ctx.shape[2]
    rng    = np.random.default_rng(seed)

    # Pre-generate all train/test splits for reproducibility
    sss_train = StratifiedShuffleSplit(n_splits=n_iter,
                                       test_size=test_size,
                                       random_state=seed)
    sss_test  = StratifiedShuffleSplit(n_splits=n_iter,
                                       test_size=test_size,
                                       random_state=seed + 1)

    splits_train = list(sss_train.split(
        np.zeros(len(y_train_ctx)), y_train_ctx))
    splits_test  = list(sss_test.split(
        np.zeros(len(y_test_ctx)),  y_test_ctx))

    results = Parallel(n_jobs=n_jobs, prefer='threads')(
        delayed(_one_iteration)(
            X_train_ctx, y_train_ctx,
            X_test_ctx,  y_test_ctx,
            splits_train[i][0],  # train index into train_ctx
            splits_test[i][1],   # test  index into test_ctx
            n_time
        )
        for i in range(n_iter)
    )
    return np.mean(results, axis=0)  # [nTimeBins × nTimeBins]


# ── Diagonal (same-time) decoding ─────────────────────────────────────────────
def same_time_decode(X_train_ctx, y_train_ctx,
                     X_test_ctx,  y_test_ctx,
                     n_iter=250, test_size=0.2, n_jobs=-1, seed=0):
    """
    Train and test at the same timepoint (diagonal of the cross-temporal matrix).
    Returns accuracy vector [nTimeBins].
    Faster than extracting the diagonal of cross_temporal_decode().
    """
    n_time = X_train_ctx.shape[2]

    sss_train = StratifiedShuffleSplit(n_splits=n_iter, test_size=test_size,
                                       random_state=seed)
    sss_test  = StratifiedShuffleSplit(n_splits=n_iter, test_size=test_size,
                                       random_state=seed + 1)

    splits_train = list(sss_train.split(np.zeros(len(y_train_ctx)), y_train_ctx))
    splits_test  = list(sss_test.split( np.zeros(len(y_test_ctx)),  y_test_ctx))

    def _one_time(t):
        accs = []
        for i in range(n_iter):
            tr_idx = splits_train[i][0]
            te_idx = splits_test[i][1]
            clf = make_clf()
            clf.fit(X_train_ctx[tr_idx, :, t], y_train_ctx[tr_idx])
            acc = clf.score(X_test_ctx[te_idx, :, t], y_test_ctx[te_idx])
            accs.append(acc)
        return np.mean(accs)

    accs = Parallel(n_jobs=n_jobs, prefer='threads')(
        delayed(_one_time)(t) for t in range(n_time)
    )
    return np.array(accs)


# ── I/O helpers ───────────────────────────────────────────────────────────────
def load_mat(path):
    """Load .mat file; handles both v5 (scipy) and v7.3 (h5py/mat73)."""
    try:
        data = sio.loadmat(path, squeeze_me=True)
        return data
    except Exception:
        try:
            import mat73
            return mat73.loadmat(path)
        except ImportError:
            raise ImportError(
                "mat73 not installed. For MATLAB v7.3 files run:\n"
                "  pip install mat73"
            )


def ensure_float32_chw(X):
    """Ensure array is float32 and shape [nTrials × nNeurons × nTimeBins]."""
    X = np.array(X, dtype=np.float32)
    if X.ndim != 3:
        raise ValueError(f"Expected 3D array [trials × neurons × time], got shape {X.shape}")
    return X


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description='Cross-context & cross-temporal value decoder')
    parser.add_argument('--input',   required=True,  help='Path to input .mat file')
    parser.add_argument('--output',  required=True,  help='Path to output .mat file')
    parser.add_argument('--n_iter',  type=int, default=250,
                        help='Number of resampling iterations (default 250)')
    parser.add_argument('--n_jobs',  type=int, default=-1,
                        help='Parallel jobs: -1=all cores (default)')
    parser.add_argument('--test_size', type=float, default=0.2,
                        help='Fraction of trials held out for testing (default 0.2)')
    parser.add_argument('--skip_crosstemporal', action='store_true',
                        help='Skip cross-temporal matrices (faster, for quick tests)')
    args = parser.parse_args()

    t0 = time.time()
    print(f"\n{'='*60}")
    print("  Cross-Decoder: Value Context")
    print(f"{'='*60}")
    print(f"  Input  : {args.input}")
    print(f"  Output : {args.output}")
    print(f"  n_iter : {args.n_iter}  |  n_jobs : {args.n_jobs}  |  test_size : {args.test_size}")
    print()

    # ── Load data ─────────────────────────────────────────────────────────────
    print("Loading data...")
    mat = load_mat(args.input)

    X_obj   = ensure_float32_chw(mat['X_obj'])    # [n_obj × nNeurons × nTimeBins]
    X_loc   = ensure_float32_chw(mat['X_loc'])    # [n_loc × nNeurons × nTimeBins]
    y_obj   = np.array(mat['y_obj'], dtype=int).ravel()
    y_loc   = np.array(mat['y_loc'], dtype=int).ravel()
    time_ax = np.array(mat['time_axis']).ravel()

    n_obj, nNeurons, nTimeBins = X_obj.shape
    n_loc = X_loc.shape[0]

    print(f"  Object ctx : {n_obj} trials  |  Location ctx : {n_loc} trials")
    print(f"  Neurons    : {nNeurons}       |  Time bins     : {nTimeBins}")
    print(f"  Time range : {time_ax[0]:.3f} – {time_ax[-1]:.3f} s\n")

    # ── Within-context same-time decoding ─────────────────────────────────────
    print("(1/5) Within-context decoding (same time, diagonal)...")
    t1 = time.time()
    acc_within_obj = same_time_decode(
        X_obj, y_obj, X_obj, y_obj,
        n_iter=args.n_iter, test_size=args.test_size, n_jobs=args.n_jobs)
    acc_within_loc = same_time_decode(
        X_loc, y_loc, X_loc, y_loc,
        n_iter=args.n_iter, test_size=args.test_size, n_jobs=args.n_jobs)
    print(f"  Done ({time.time()-t1:.1f}s)  "
          f"peak obj={acc_within_obj.max():.3f}  peak loc={acc_within_loc.max():.3f}")

    # ── Cross-context same-time decoding ──────────────────────────────────────
    print("(2/5) Cross-context same-time decoding...")
    t1 = time.time()
    acc_cross_obj2loc = same_time_decode(
        X_obj, y_obj, X_loc, y_loc,
        n_iter=args.n_iter, test_size=args.test_size, n_jobs=args.n_jobs)
    acc_cross_loc2obj = same_time_decode(
        X_loc, y_loc, X_obj, y_obj,
        n_iter=args.n_iter, test_size=args.test_size, n_jobs=args.n_jobs)
    print(f"  Done ({time.time()-t1:.1f}s)  "
          f"peak obj→loc={acc_cross_obj2loc.max():.3f}  "
          f"peak loc→obj={acc_cross_loc2obj.max():.3f}")

    if args.skip_crosstemporal:
        print("\n  [skipping cross-temporal matrices as requested]")
        acc_ct_obj   = np.full((nTimeBins, nTimeBins), np.nan, dtype=np.float32)
        acc_ct_loc   = np.full((nTimeBins, nTimeBins), np.nan, dtype=np.float32)
        acc_ct_cross = np.full((nTimeBins, nTimeBins), np.nan, dtype=np.float32)
    else:
        # ── Cross-temporal within object context ──────────────────────────────
        print("(3/5) Cross-temporal: within object context  [T×T matrix]...")
        t1 = time.time()
        acc_ct_obj = cross_temporal_decode(
            X_obj, y_obj, X_obj, y_obj,
            n_iter=args.n_iter, test_size=args.test_size, n_jobs=args.n_jobs)
        print(f"  Done ({time.time()-t1:.1f}s)")

        # ── Cross-temporal within location context ────────────────────────────
        print("(4/5) Cross-temporal: within location context  [T×T matrix]...")
        t1 = time.time()
        acc_ct_loc = cross_temporal_decode(
            X_loc, y_loc, X_loc, y_loc,
            n_iter=args.n_iter, test_size=args.test_size, n_jobs=args.n_jobs)
        print(f"  Done ({time.time()-t1:.1f}s)")

        # ── Cross-temporal cross-context (train obj → test loc) ───────────────
        print("(5/5) Cross-temporal cross-context: train obj → test loc  [T×T]...")
        t1 = time.time()
        acc_ct_cross = cross_temporal_decode(
            X_obj, y_obj, X_loc, y_loc,
            n_iter=args.n_iter, test_size=args.test_size, n_jobs=args.n_jobs)
        print(f"  Done ({time.time()-t1:.1f}s)")

    # ── Save ──────────────────────────────────────────────────────────────────
    print(f"\nSaving to {args.output}...")
    sio.savemat(args.output, {
        'acc_within_obj':         acc_within_obj,
        'acc_within_loc':         acc_within_loc,
        'acc_cross_obj2loc':      acc_cross_obj2loc,
        'acc_cross_loc2obj':      acc_cross_loc2obj,
        'acc_crosstemporal_obj':  acc_ct_obj,
        'acc_crosstemporal_loc':  acc_ct_loc,
        'acc_crosstemporal_cross': acc_ct_cross,
        'time_axis':              time_ax,
        'chance':                 0.5,
        'n_iter':                 args.n_iter,
        'n_neurons':              nNeurons,
        'n_time':                 nTimeBins,
    })

    total = time.time() - t0
    print(f"\nTotal time: {total:.1f}s  ({total/60:.1f} min)")
    print("Done.\n")


if __name__ == '__main__':
    main()
