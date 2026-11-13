"""Chebyshev interpolation in w: first-kind nodes, barycentric basis, its derivative, tail estimate."""
import numpy as np


def nodes(K, lo, hi):
    k = np.arange(K)
    x = np.cos((2 * k + 1) * np.pi / (2 * K))
    lam = (-1.0) ** k * np.sin((2 * k + 1) * np.pi / (2 * K))
    return lo + (hi - lo) * (x + 1) / 2, x, lam


def basis(wq, K, lo, hi):
    """(len(wq), K) barycentric Lagrange basis at query points wq."""
    _, x, lam = nodes(K, lo, hi)
    if hi == lo:
        return np.full((len(np.atleast_1d(wq)), K), 1.0 / K)
    xq = 2 * (np.asarray(wq, float) - lo) / (hi - lo) - 1
    D = xq[:, None] - x[None, :]
    hit = np.abs(D) < 1e-15
    with np.errstate(divide='ignore', invalid='ignore'):
        Wt = lam[None, :] / D
    Wt = np.where(hit, 0.0, Wt)
    B = Wt / Wt.sum(1, keepdims=True)
    rows = hit.any(1)
    B[rows] = hit[rows].astype(float)
    return B


def diff_matrix(K, lo, hi):
    """D with (D f)_i = p'(w_i) for the interpolant p through (w_k, f_k), in w units."""
    _, x, lam = nodes(K, lo, hi)
    X = x[:, None] - x[None, :]
    np.fill_diagonal(X, 1.0)
    D = (lam[None, :] / lam[:, None]) / X
    np.fill_diagonal(D, 0.0)
    np.fill_diagonal(D, -D.sum(1))
    return D * (2.0 / (hi - lo)) if hi > lo else D * 0.0


def basis_deriv(wq, K, lo, hi):
    """(len(wq), K) derivative basis: p'(w) = basis_deriv(w) @ f."""
    return basis(wq, K, lo, hi) @ diff_matrix(K, lo, hi)


def coeffs(F):
    """Chebyshev coefficients of the interpolant from values at first-kind nodes, F: (..., K)."""
    K = F.shape[-1]
    k = np.arange(K)
    T = np.cos(np.pi * np.outer(np.arange(K), (2 * k + 1)) / (2 * K))     # T_j(x_k)
    c = F @ T.T * (2.0 / K)
    c[..., 0] /= 2
    return c


def tail(F):
    """Relative size of the last two Chebyshev coefficients, max over rows: an error estimate."""
    c = np.abs(coeffs(F))
    return float((c[..., -2:].sum(-1) / np.abs(F).max(-1)).max())


def tail_min(F, tol=1e-11, noise_c=64.0, err_cap=1e-10):
    """Tail test value for node values F (rows = patterns, columns = K first-kind nodes): max over rows of
    tail_u / threshold_u, where tail_u = (|c_K-2| + |c_K-1|) / min_k F_uk and threshold_u =
    max(tol, min(noise_c eps max_k F_uk / min_k F_uk, err_cap)) (see psigc.gpu.tail_threshold). Passes when
    <= 1; inf if a row has a non-positive node value."""
    F = np.asarray(F, float)
    c = np.abs(coeffs(F))
    mn = F.min(-1)
    if (mn <= 0).any():
        return float('inf')
    thr = np.maximum(tol, np.minimum(noise_c * np.finfo(float).eps * F.max(-1) / mn, err_cap))
    return float((c[..., -2:].sum(-1) / mn / thr).max())
